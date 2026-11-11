from __future__ import annotations

import csv
import itertools
import json
import math
import os
import shutil
import subprocess
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass
from pathlib import Path
from statistics import fmean, pstdev
from tempfile import TemporaryDirectory
from time import perf_counter
from typing import Any

import gemmi

from .scoring import SimRNAResult, _engine_provenance, _native_score, _prepare_input, _sha256

BASES = "ACGU"
Contact = tuple[str, str, str, str]
State = tuple[tuple[int, str], ...]


@dataclass(frozen=True)
class SimRNAFrustrationResult:
    frustration_file: Path
    variants_file: Path
    run_file: Path
    contact_count: int
    engine_calls: int


def _base_atoms(residue: Any) -> tuple[Any, Any, Any]:
    names = ("N9", "C2", "C6") if residue.name in "AG" else ("N1", "C2", "C4")
    return tuple(next(a for a in residue if a.name == name) for name in names)


def _frame(atoms: tuple[Any, Any, Any]) -> tuple[Any, tuple[Any, Any, Any]]:
    anchor, c2, edge = (atom.pos for atom in atoms)
    origin = (anchor + edge) * 0.5
    x = (edge - anchor).normalized()
    normal = x.cross(c2 - anchor)
    if not math.isfinite(normal.length()) or normal.length() < 1e-6:
        raise ValueError("Degenerate base triangle.")
    z = normal.normalized()
    return origin, (x, z.cross(x), z)


def _replacement(residue: Any, template: Any) -> Any:
    replacement = residue.clone()
    if residue.name == template.name:
        return replacement
    origin, axes = _frame(_base_atoms(residue))
    template_atoms = _base_atoms(template)
    template_origin, template_axes = _frame(template_atoms)
    for index in reversed(range(len(replacement))):
        if replacement[index].name not in {"P", "C4'"}:
            del replacement[index]
    replacement.name = template.name
    for template_atom in template_atoms:
        displacement = template_atom.pos - template_origin
        xyz = gemmi.Vec3(origin.x, origin.y, origin.z)
        for axis, template_axis in zip(axes, template_axes, strict=True):
            xyz = xyz + axis * displacement.dot(template_axis)
        atom = template_atom.clone()
        atom.pos = gemmi.Position(xyz.x, xyz.y, xyz.z)
        replacement.add_atom(atom)
    return replacement


def _write_variant(
    structure: Any,
    locations: list[tuple[int, int]],
    replacements: dict[tuple[int, str], Any],
    state: State,
    path: Path,
) -> None:
    mutant = structure.clone()
    for index, base in state:
        chain_index, residue_index = locations[index]
        chain = mutant[0][chain_index]
        del chain[residue_index]
        chain.add_residue(replacements[index, base].clone(), pos=residue_index)
    mutant.write_pdb(str(path))


def _state(i: int, j: int, a: str, b: str, native: list[str]) -> State:
    return tuple(sorted((k, base) for k, base in ((i, a), (j, b)) if base != native[k]))


def _scan(
    pdb: Path, templates_pdb: Path, output: Path, vendor: Path, options: dict[str, Any]
) -> None:
    start = perf_counter()
    with TemporaryDirectory(prefix="simrna-scan-") as temporary:
        work = Path(temporary)
        local_vendor = work / "engine"
        local_vendor.mkdir()
        for name in ("SimRNA", "config_n0.dat"):
            shutil.copy2(vendor / name, local_vendor / name)
        shutil.copytree(vendor / "data", local_vendor / "data")
        provenance = _engine_provenance(local_vendor)

        reference_path = work / "reference.pdb"
        _prepare_input(pdb, reference_path)
        structure = gemmi.read_structure(str(reference_path))
        template_path = work / "templates.pdb"
        _prepare_input(templates_pdb, template_path)
        templates: dict[str, Any] = {}
        template_ids = {}
        for chain in gemmi.read_structure(str(template_path))[0]:
            for residue in chain:
                if residue.name not in templates:
                    templates[residue.name] = residue.clone()
                    template_ids[residue.name] = [chain.name, str(residue.seqid)]
        if set(templates) != set(BASES):
            raise ValueError("templates_pdb must contain all four canonical bases.")

        locations = [
            (ci, ri) for ci, chain in enumerate(structure[0]) for ri, _ in enumerate(chain)
        ]
        residues = [structure[0][ci][ri] for ci, ri in locations]
        identifiers = [
            (structure[0][ci].name, str(residue.seqid))
            for (ci, _), residue in zip(locations, residues, strict=True)
        ]
        lookup = {key: index for index, key in enumerate(identifiers)}
        native = [residue.name for residue in residues]
        explicit = options["contacts"]
        if explicit is None:
            centers = [_frame(_base_atoms(residue))[0] for residue in residues]
            pairs = [
                (i, j)
                for i, j in itertools.combinations(range(len(residues)), 2)
                if centers[i].dist(centers[j]) <= options["cutoff"]
            ]
        else:
            pairs = [(lookup[(ci, ri)], lookup[(cj, rj)]) for ci, ri, cj, rj in explicit]
        if not pairs or any(i == j for i, j in pairs):
            raise ValueError("Select at least one pair of distinct residues.")
        if len({tuple(sorted(pair)) for pair in pairs}) != len(pairs):
            raise ValueError("Duplicate contact.")

        replacements = {
            (index, base): _replacement(residues[index], templates[base])
            for index in {index for pair in pairs for index in pair}
            for base in BASES
            if base != native[index]
        }
        combinations = list(itertools.product(BASES, repeat=2))
        states = sorted({_state(i, j, a, b, native) for i, j in pairs for a, b in combinations})
        preparation_seconds = perf_counter() - start

        reference_dir = output / "reference"
        reference_dir.mkdir()
        reference = _native_score(
            reference_path, reference_dir, local_vendor, 1, options["timeout"], provenance
        )
        scores: dict[State, tuple[float, dict[str, float]]] = {
            (): (reference.total_energy, reference.components)
        }
        workers = min(options["workers"], os.cpu_count() or 1)
        scoring_start = perf_counter()

        def score(state: State) -> tuple[State, tuple[float, dict[str, float]]]:
            with TemporaryDirectory(dir=work, prefix="variant-") as folder:
                variant_dir = Path(folder)
                variant_path = variant_dir / "variant.pdb"
                _write_variant(structure, locations, replacements, state, variant_path)
                score_dir = variant_dir / "score"
                score_dir.mkdir()
                try:
                    result: SimRNAResult = _native_score(
                        variant_path, score_dir, local_vendor, 1, options["timeout"], provenance
                    )
                except Exception as error:
                    failure = output / ("failed_" + "_".join(f"{i}{b}" for i, b in state))
                    shutil.copytree(variant_dir, failure)
                    raise RuntimeError(
                        f"SimRNA failed for state {state}. See {failure.name}."
                    ) from error
                return state, (result.total_energy, result.components)

        with ThreadPoolExecutor(max_workers=workers) as pool:
            scores.update(pool.map(score, (state for state in states if state)))

    component_names = list(reference.components)
    with (
        (output / "frustration.csv").open("w", newline="", encoding="utf-8") as f,
        (output / "contact_variants.csv").open("w", newline="", encoding="utf-8") as v,
    ):
        summary, variants = csv.writer(f), csv.writer(v)
        ids = ["chain_i", "residue_i", "chain_j", "residue_j"]
        summary.writerow([*ids, "E_native", "E_mean", "E_sd", "F"])
        variants.writerow([*ids, "base_i", "base_j", "E_total", "delta_E_total", *component_names])
        for i, j in pairs:
            contact = (*identifiers[i], *identifiers[j])
            deltas = []
            for a, b in combinations:
                energy, parts = scores[_state(i, j, a, b, native)]
                delta = energy - reference.total_energy
                deltas.append(delta)
                variants.writerow(
                    [*contact, a, b, energy, delta, *(parts[k] for k in component_names)]
                )
            mean, sigma = fmean(deltas), pstdev(deltas)
            summary.writerow(
                [
                    *contact,
                    reference.total_energy,
                    reference.total_energy + mean,
                    sigma,
                    mean / sigma if sigma else "",
                ]
            )

    (output / "run.json").write_text(
        json.dumps(
            {
                **provenance,
                "engine": "SimRNA 3.20",
                "iterations": 0,
                "seed": 1,
                "contact_selection": "explicit"
                if explicit is not None
                else "base_midpoint_distance",
                "contact_cutoff_angstrom": options["cutoff"] if explicit is None else None,
                "mutation_geometry": (
                    "fixed N-to-edge midpoint and oriented base plane; template triangle"
                ),
                "reference_sha256": _sha256(pdb),
                "template_sha256": _sha256(templates_pdb),
                "template_residues": template_ids,
                "contact_count": len(pairs),
                "engine_calls": len(states),
                "workers": workers,
                "preparation_seconds": preparation_seconds,
                "variant_scoring_seconds": perf_counter() - scoring_start,
                "worker_seconds": perf_counter() - start,
            },
            indent=2,
        ),
        encoding="utf-8",
    )


def _container_main() -> None:
    options = json.loads(Path("/options.json").read_text(encoding="utf-8"))
    _scan(Path("/input.pdb"), Path("/templates.pdb"), Path("/output"), Path("/simrna"), options)


def run_simrna_frustration(
    pdb_path: str | Path,
    output_dir: str | Path,
    *,
    simrna_dir: str | Path,
    templates_pdb: str | Path,
    contacts: list[Contact] | None = None,
    cutoff: float = 8.0,
    workers: int = 4,
    timeout_per_variant: float = 120,
    docker_image: str | None = "simrna-adapter:local",
) -> SimRNAFrustrationResult:
    """Run an independent SimRNA contact scan with cached sequence states."""
    pdb = Path(pdb_path).resolve(strict=True)
    templates = Path(templates_pdb).resolve(strict=True)
    vendor = Path(simrna_dir).resolve(strict=True)
    output = Path(output_dir).resolve()
    if workers < 1 or not all(
        math.isfinite(value) and value > 0 for value in (cutoff, timeout_per_variant)
    ):
        raise ValueError("Workers, cutoff and timeout must be positive.")
    if pdb.suffix.lower() != ".pdb" or templates.suffix.lower() != ".pdb":
        raise ValueError("Input and templates must be PDB files.")
    output.mkdir(parents=True, exist_ok=True)
    if any(output.iterdir()):
        raise ValueError(f"Output directory must be empty: {output}")
    options = dict(contacts=contacts, cutoff=cutoff, workers=workers, timeout=timeout_per_variant)
    if docker_image is None:
        _scan(pdb, templates, output, vendor, options)
    else:
        with TemporaryDirectory(prefix="simrna-options-") as temporary:
            manifest = Path(temporary) / "options.json"
            manifest.write_text(json.dumps(options), encoding="utf-8")
            command = ["docker", "run", "--rm", "--network", "none", "--entrypoint", "python"]
            for path, target in [
                (pdb, "/input.pdb"),
                (templates, "/templates.pdb"),
                (vendor, "/simrna"),
                (manifest, "/options.json"),
            ]:
                command.extend(["--mount", f"type=bind,source={path},target={target},readonly"])
            command.extend(
                [
                    "--mount",
                    f"type=bind,source={output},target=/output",
                    docker_image,
                    "-c",
                    "from simrna_adapter.frustration import _container_main; _container_main()",
                ]
            )
            process = subprocess.run(command, capture_output=True, text=True, check=False)
            (output / "docker.log").write_text(process.stdout + process.stderr, encoding="utf-8")
            if process.returncode:
                raise RuntimeError(f"SimRNA scan failed. See {output / 'docker.log'}.")
    run_file = output / "run.json"
    result = json.loads(run_file.read_text(encoding="utf-8"))
    return SimRNAFrustrationResult(
        output / "frustration.csv",
        output / "contact_variants.csv",
        run_file,
        result["contact_count"],
        result["engine_calls"],
    )
