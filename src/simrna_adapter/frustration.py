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

from .scoring import (
    SimRNAResult,
    _beads,
    _engine_provenance,
    _native_score,
    _prepare_input,
    _sha256,
)

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

    @property
    def scores_file(self) -> Path:
        return self.run_file.parent / "contact_scores.csv"


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


def _distance_contacts(structure: Any, cutoff: float) -> list[Contact]:
    if not math.isfinite(cutoff) or cutoff <= 0:
        raise ValueError("Contact cutoff must be finite and positive.")
    residues = [
        (chain.name, str(residue.seqid), residue) for chain in structure[0] for residue in chain
    ]
    centers = [_frame(_base_atoms(r))[0] for _, _, r in residues]
    return [
        (ci, ri, cj, rj)
        for i, (ci, ri, _) in enumerate(residues)
        for j, (cj, rj, _) in enumerate(residues)
        if j > i and centers[i].dist(centers[j]) <= cutoff
    ]


def select_simrna_contacts(pdb_path: str | Path, *, cutoff: float = 8.0) -> list[Contact]:
    """Return the adapter's native base-midpoint contacts without running SimRNA.

    The point is (N9+C6)/2 for A/G and (N1+C4)/2 for C/U. This threshold
    selects mutation targets, not the licensed engine's interaction range.
    """
    path = Path(pdb_path).resolve(strict=True)
    if path.suffix.lower() != ".pdb":
        raise ValueError("Input must be a PDB file.")
    structure = gemmi.read_structure(str(path))
    _beads(structure)
    return _distance_contacts(structure, cutoff)


SUGAR_ATOMS = {
    "P",
    "OP1",
    "OP2",
    "OP3",
    "O5'",
    "C5'",
    "C4'",
    "O4'",
    "C3'",
    "O3'",
    "C2'",
    "O2'",
    "C1'",
}
BASE_ATOMS = {
    "A": {"N9", "C8", "N7", "C5", "C6", "N6", "N1", "C2", "N3", "C4"},
    "G": {"N9", "C8", "N7", "C5", "C6", "O6", "N1", "C2", "N2", "N3", "C4"},
    "C": {"N1", "C2", "O2", "N3", "C4", "N4", "C5", "C6"},
    "U": {"N1", "C2", "O2", "N3", "C4", "O4", "C5", "C6"},
}


def _glycosidic_frame(residue: Any) -> tuple[Any, tuple[Any, Any, Any]]:
    """Anchor at N and preserve N-to-C1' direction and the chi-defining base plane."""
    names = ("N9", "C1'", "C4") if residue.name in "AG" else ("N1", "C1'", "C2")
    points = []
    for name in names:
        atoms = [atom for atom in residue if atom.name == name]
        if len(atoms) != 1 or atoms[0].has_altloc() or atoms[0].occ <= 0:
            raise ValueError(
                f"Mutation geometry requires one occupied {name} at {residue.seqid}. "
                "Use a full-atom reference/template."
            )
        points.append(atoms[0].pos)
    anchor, sugar, base = points
    bond = sugar - anchor
    if bond.length() < 1e-6:
        raise ValueError("Degenerate glycosidic bond.")
    x = bond.normalized()
    normal = x.cross(base - anchor)
    if normal.length() < 1e-6:
        raise ValueError("Degenerate glycosidic frame.")
    z = normal.normalized()
    return anchor, (x, z.cross(x), z)


def _replacement(residue: Any, template: Any, *, preserve_same_class_ring: bool = False) -> Any:
    replacement = residue.clone()
    if residue.name == template.name:
        return replacement
    origin, axes = _glycosidic_frame(residue)
    template_origin, template_axes = _glycosidic_frame(template)
    for index in reversed(range(len(replacement))):
        if replacement[index].name not in SUGAR_ATOMS:
            del replacement[index]
    replacement.name = template.name
    for template_atom in template:
        if template_atom.name not in BASE_ATOMS[template.name]:
            continue
        displacement = template_atom.pos - template_origin
        xyz = gemmi.Vec3(origin.x, origin.y, origin.z)
        for axis, template_axis in zip(axes, template_axes, strict=True):
            xyz = xyz + axis * displacement.dot(template_axis)
        atom = template_atom.clone()
        atom.pos = gemmi.Position(xyz.x, xyz.y, xyz.z)
        replacement.add_atom(atom)
    if preserve_same_class_ring and (residue.name in "AG") == (template.name in "AG"):
        # Same-class mutations need no change to the ring's observed coordinates.
        # Exocyclic groups still follow the target identity; this is not a clash fix.
        ring = (
            {"N9", "C8", "N7", "C5", "C6", "N1", "C2", "N3", "C4"}
            if residue.name in "AG"
            else {"N1", "C2", "N3", "C4", "C5", "C6"}
        )
        parents = {"N6": "C6", "O6": "C6", "N2": "C2", "O2": "C2", "N4": "C4", "O4": "C4"}
        aligned = {atom.name: atom.pos for atom in replacement}
        native = {atom.name: atom.pos for atom in residue}
        for atom in replacement:
            if atom.name in ring and atom.name in native:
                atom.pos = native[atom.name]
            elif atom.name in parents:
                parent = parents[atom.name]
                atom.pos = native[parent] + (aligned[atom.name] - aligned[parent])
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


def _isolate_pair(structure: Any, pair: tuple[int, int]) -> Any:
    """Keep full chain topology; separate every other residue beyond interaction range.

    The selected pair keeps its exact coordinates. Only nonlocal energy is read
    from this diagnostic structure; its artificial backbone strain is discarded.
    """
    masked = structure.clone()
    residues = [r for c in masked[0] for r in c]
    offset = max(1000.0, max(a.pos.x for r in residues for a in r) + 100.0)
    for index, residue in enumerate(residues):
        if index in pair:
            continue
        anchor = residue[0].pos
        shift = gemmi.Vec3(
            offset + 50 * (index % 40) - anchor.x,
            1000 + 50 * ((index // 40) % 40) - anchor.y,
            1000 + 50 * (index // 1600) - anchor.z,
        )
        for atom in residue:
            xyz = gemmi.Position(atom.pos.x + shift.x, atom.pos.y + shift.y, atom.pos.z + shift.z)
            if any(not -999 < v < 9999 for v in (xyz.x, xyz.y, xyz.z)):
                raise ValueError("Pair isolation exceeds PDB coordinate limits.")
            atom.pos = xyz
    return masked


def _pair_energy(result: SimRNAResult) -> float:
    # Both values are printed to six decimals; component subtotals have only three.
    energy = result.total_energy - result.components["E_local_geometry"]
    subtotal = result.components["E_base_base"] + result.components["E_base_backbone"]
    if not math.isclose(energy, subtotal, rel_tol=0.0, abs_tol=0.002):
        raise ValueError("Isolated nonlocal energy does not match component subtotals.")
    return energy


def _energy_statistics(
    native: float, values: list[float], *, resolution: float = 0.0
) -> list[float | str]:
    if len(values) != 16 or not all(math.isfinite(v) for v in [native, *values]):
        raise ValueError("Expected 16 finite variant energies.")
    deltas = [value - native for value in values]
    mean, sigma = fmean(deltas), pstdev(deltas)
    return [native, native + mean, sigma, mean / sigma if sigma > resolution else ""]


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
        chemical_structure = gemmi.read_structure(str(pdb))
        template_path = work / "templates.pdb"
        _prepare_input(templates_pdb, template_path)
        templates: dict[str, Any] = {}
        template_ids = {}
        for chain in gemmi.read_structure(str(templates_pdb))[0]:
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
            pairs = [
                (lookup[(ci, ri)], lookup[(cj, rj)])
                for ci, ri, cj, rj in _distance_contacts(structure, options["cutoff"])
            ]
        else:
            pairs = [(lookup[(ci, ri)], lookup[(cj, rj)]) for ci, ri, cj, rj in explicit]
        if not pairs or any(i == j for i, j in pairs):
            raise ValueError("Select at least one pair of distinct residues.")
        if len({tuple(sorted(pair)) for pair in pairs}) != len(pairs):
            raise ValueError("Duplicate contact.")

        replacements = {
            (index, base): _replacement(
                chemical_structure[0][locations[index][0]][locations[index][1]],
                templates[base],
                preserve_same_class_ring=options.get("preserve_same_class_ring", False),
            )
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

        pair_scores: dict[tuple[int, int, str, str], float] = {}
        pair_calls = 0
        if options.get("pair_energies", True):
            tasks = [(i, j, a, b) for i, j in pairs for a, b in combinations]

            def score_pair(
                task: tuple[int, int, str, str],
            ) -> tuple[tuple[int, int, str, str], float]:
                i, j, a, b = task
                with TemporaryDirectory(dir=work, prefix="pair-") as folder:
                    variant_dir = Path(folder)
                    path = variant_dir / "pair.pdb"
                    masked = _isolate_pair(structure, (i, j))
                    _write_variant(
                        masked, locations, replacements, _state(i, j, a, b, native), path
                    )
                    score_dir = variant_dir / "score"
                    score_dir.mkdir()
                    try:
                        result = _native_score(
                            path, score_dir, local_vendor, 1, options["timeout"], provenance
                        )
                        energy = _pair_energy(result)
                    except Exception as error:
                        failure = output / f"failed_pair_{i}_{j}_{a}{b}"
                        shutil.copytree(variant_dir, failure)
                        raise RuntimeError(f"Pair scoring failed; see {failure.name}.") from error
                return task, energy

            with ThreadPoolExecutor(max_workers=workers) as pool:
                pair_scores.update(pool.map(score_pair, tasks))
            pair_calls = len(tasks)

    component_names = list(reference.components)
    ids = ["chain_i", "residue_i", "chain_j", "residue_j"]
    with (
        (output / "frustration.csv").open("w", newline="", encoding="utf-8") as f,
        (output / "contact_variants.csv").open("w", newline="", encoding="utf-8") as v,
        (output / "contact_scores.csv").open("w", newline="", encoding="utf-8") as s,
    ):
        summary, variants, dual = csv.writer(f), csv.writer(v), csv.writer(s)
        summary.writerow([*ids, "E_native", "E_mean", "E_sd", "F", "F_total", "F_pair"])
        dual.writerow(
            [
                *ids,
                "native_E_pair",
                "mean_E_pair",
                "sd_E_pair",
                "F_pair",
                "native_E_total",
                "mean_E_total",
                "sd_E_total",
                "F_total",
            ]
        )
        variants.writerow(
            [
                *ids,
                "base_i",
                "base_j",
                "E_total",
                "delta_E_total",
                *component_names,
                "native_E_total",
                "native_E_pair",
                "E_pair",
                "delta_E_pair",
            ]
        )
        for i, j in pairs:
            contact = (*identifiers[i], *identifiers[j])
            pair_native = pair_scores.get((i, j, native[i], native[j]))
            totals, pair_values = [], []
            for a, b in combinations:
                energy, parts = scores[_state(i, j, a, b, native)]
                totals.append(energy)
                pair_value = pair_scores.get((i, j, a, b))
                if pair_value is not None:
                    pair_values.append(pair_value)
                variants.writerow(
                    [
                        *contact,
                        a,
                        b,
                        energy,
                        energy - reference.total_energy,
                        *(parts[k] for k in component_names),
                        reference.total_energy,
                        pair_native if pair_native is not None else "",
                        pair_value if pair_value is not None else "",
                        pair_value - pair_native
                        if pair_value is not None and pair_native is not None
                        else "",
                    ]
                )
            total_stats = _energy_statistics(reference.total_energy, totals)
            pair_stats = (
                _energy_statistics(pair_native, pair_values, resolution=2e-6)
                if pair_native is not None
                else [""] * 4
            )
            dual.writerow([*contact, *pair_stats, *total_stats])
            summary.writerow([*contact, *total_stats, total_stats[-1], pair_stats[-1]])

    (output / "run.json").write_text(
        json.dumps(
            {
                **provenance,
                "engine": "SimRNA 3.20",
                "iterations": 0,
                "mutation_protocol_version": 4
                if options.get("preserve_same_class_ring", False)
                else 3,
                "preserve_same_class_ring": options.get("preserve_same_class_ring", False),
                "adapter_source_sha256": _sha256(Path(__file__)),
                "seed": 1,
                "contact_selection": "explicit"
                if explicit is not None
                else "base_midpoint_distance",
                "contact_cutoff_angstrom": options["cutoff"] if explicit is None else None,
                "mutation_geometry": (
                    "fixed native sugar and N1/N9; N-C1 direction and glycosidic chi preserved; "
                    + (
                        "native ring retained for A/G and C/U exchanges"
                        if options.get("preserve_same_class_ring", False)
                        else "all changed bases placed from templates"
                    )
                ),
                "reference_sha256": _sha256(pdb),
                "template_sha256": _sha256(templates_pdb),
                "template_residues": template_ids,
                "contact_count": len(pairs),
                "engine_calls": len(states) + pair_calls,
                "total_engine_calls": len(states),
                "pair_engine_calls": pair_calls,
                "energy_strategies": ["total", "pair"] if pair_scores else ["total"],
                "pair_energy_definition": (
                    "Nonlocal interaction of i and j: base-base plus both directions of "
                    "base-backbone; full chain topology retained, all other residues separated. "
                    "Local backbone strain excluded; stacking subset not added twice."
                )
                if pair_scores
                else None,
                "pair_F_minimum_sd": 2e-6,
                "pair_energy_numerical_precision": "approximately 2e-6 (engine log rounding)",
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
    pair_energies: bool = True,
    preserve_same_class_ring: bool = False,
    docker_image: str | None = "simrna-adapter:local",
) -> SimRNAFrustrationResult:
    """Score both E_total and E_pair by default. Disable pair_energies for total-only speed.

    Opt in to native ring coordinates for A<->G and C<->U exchanges with
    preserve_same_class_ring=True. Purine/pyrimidine exchanges keep the existing
    glycosidic template mapping. This does not minimize energy or guarantee no clashes.

    Pair energy includes nonlocal base-base and base-backbone interactions,
    excluding local backbone geometry. It needs 16 extra engine calls per pair.
    """
    if not isinstance(preserve_same_class_ring, bool):
        raise ValueError("preserve_same_class_ring must be a bool.")
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
    options = dict(
        contacts=contacts,
        cutoff=cutoff,
        workers=workers,
        timeout=timeout_per_variant,
        pair_energies=pair_energies,
        preserve_same_class_ring=preserve_same_class_ring,
    )
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
    if preserve_same_class_ring and (
        result.get("preserve_same_class_ring") is not True
        or result.get("mutation_protocol_version") != 4
    ):
        raise RuntimeError(
            "SimRNA image did not apply preserve_same_class_ring. Rebuild the adapter image."
        )
    return SimRNAFrustrationResult(
        output / "frustration.csv",
        output / "contact_variants.csv",
        run_file,
        result["contact_count"],
        result["engine_calls"],
    )
