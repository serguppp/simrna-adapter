from __future__ import annotations

import hashlib
import json
import math
import os
import re
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path
from tempfile import TemporaryDirectory
from time import perf_counter
from typing import Any

import gemmi


@dataclass(frozen=True)
class SimRNAResult:
    total_energy: float
    components: dict[str, float]
    output_dir: Path
    run_file: Path


def _sha256(path: Path) -> str:
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def _beads(structure: Any) -> dict[tuple[str, str, str], tuple[float, float, float]]:
    if len(structure) != 1:
        raise ValueError("Provide exactly one PDB model.")
    positions = {}
    for chain in structure[0]:
        for residue in chain:
            if residue.name not in ("A", "C", "G", "U"):
                raise ValueError(
                    f"Noncanonical/non-RNA residue: {chain.name}:{residue.seqid} "
                    f"{residue.name}. Prepare the input explicitly first."
                )
            names = (
                ("P", "C4'", "N9", "C2", "C6")
                if residue.name in ("A", "G")
                else ("P", "C4'", "N1", "C2", "C4")
            )
            for name in names:
                atoms = [atom for atom in residue if atom.name == name]
                if len(atoms) != 1 or atoms[0].has_altloc() or atoms[0].occ <= 0:
                    raise ValueError(
                        f"Expected one occupied {name} without altloc: "
                        f"{chain.name}:{residue.seqid} {residue.name}"
                    )
                pos = atoms[0].pos
                coords = (pos.x, pos.y, pos.z)
                if not all(math.isfinite(value) for value in coords):
                    raise ValueError("Nonfinite coordinates.")
                key = (chain.name, str(residue.seqid), name)
                if key in positions:
                    raise ValueError(f"Duplicate residue/atom identifier: {key}")
                positions[key] = coords
    if not positions:
        raise ValueError("No RNA nucleotides found.")
    return positions


def _prepare_input(source: Path, destination: Path) -> dict[str, Any]:
    if source.suffix.lower() != ".pdb":
        raise ValueError("Provide a PDB file.")
    structure = gemmi.read_structure(str(source))
    positions = _beads(structure)
    for chain in structure[0]:
        for residue in chain:
            required = (
                {"P", "C4'", "N9", "C2", "C6"}
                if residue.name in ("A", "G")
                else {"P", "C4'", "N1", "C2", "C4"}
            )
            for index in range(len(residue) - 1, -1, -1):
                if residue[index].name not in required:
                    del residue[index]
    structure.write_pdb(str(destination))
    return {"nucleotides": len(positions) // 5, "beads": len(positions)}


def _parse_output(log: str, trajectory: str) -> tuple[float, dict[str, float]]:
    lines = [line for line in trajectory.splitlines() if line.strip()]
    if len(lines) != 2:
        raise ValueError("Expected exactly one trajectory frame (header and coordinates).")
    header = lines[0].split()
    if len(header) != 5:
        raise ValueError("Expected five columns in the trajectory header.")
    restrained, energy = float(header[2]), float(header[3])
    if not all(math.isfinite(value) for value in (restrained, energy)):
        raise ValueError("Nonfinite SimRNA energy.")
    if not math.isclose(restrained, energy, abs_tol=1e-5):
        raise ValueError("Unexpected restraint contribution in an unrestrained score.")
    matches = re.findall(r"Total energy of system\s+([-+0-9.eE]+)", log)
    if len(matches) != 1 or not math.isclose(float(matches[0]), energy, abs_tol=2e-5):
        raise ValueError("Log and trajectory energies disagree or are missing.")
    patterns = {
        "E_base_base": r"Base-Base interactions energy\s+([-+0-9.eE]+)",
        "E_base_backbone": r"Base-Backbone interact\. energy\s+([-+0-9.eE]+)",
        "E_local_geometry": r"Local geometry energy\s+([-+0-9.eE]+)",
        "E_short_stacking_subset": r"short stacking energy\s+([-+0-9.eE]+)",
    }
    components = {}
    for name, pattern in patterns.items():
        found = re.findall(pattern, log)
        if len(found) != 1 or not math.isfinite(float(found[0])):
            raise ValueError(f"Missing/invalid energy component: {name}")
        components[name] = float(found[0])
    main_sum = sum(
        components[name] for name in ("E_base_base", "E_base_backbone", "E_local_geometry")
    )
    if not math.isclose(main_sum, energy, abs_tol=0.002):
        raise ValueError("Main energy components do not sum to E_total.")
    return energy, components


def _check_geometry(input_pdb: Path, reduced_pdb: Path, trajectory: Path) -> float:
    before = _beads(gemmi.read_structure(str(input_pdb)))
    after = _beads(gemmi.read_structure(str(reduced_pdb)))
    if before.keys() != after.keys():
        raise ValueError("SimRNA changed residue/atom identifiers.")
    displacement = max(math.dist(before[key], after[key]) for key in before)
    if displacement > 0.002:
        raise ValueError(f"Zero-step run moved an input bead by {displacement} angstrom.")
    coordinates = [float(value) for value in trajectory.read_text().splitlines()[1].split()]
    if len(coordinates) != 3 * len(before) or not all(map(math.isfinite, coordinates)):
        raise ValueError("Incomplete/nonfinite trajectory coordinates.")
    # TRAFL and the reduced PDB use P, C4', N9/N1, C2, C6/C4 order.
    expected = [component for xyz in after.values() for component in xyz]
    if any(abs(a - b) > 0.002 for a, b in zip(coordinates, expected, strict=True)):
        raise ValueError("Trajectory and reduced-PDB coordinates disagree.")
    return displacement


def _engine_provenance(vendor: Path) -> dict[str, Any]:
    return {
        "data_files_sha256": {
            str(path.relative_to(vendor / "data")): _sha256(path)
            for path in sorted((vendor / "data").rglob("*"))
            if path.is_file()
        },
        "binary_sha256": _sha256(vendor / "SimRNA"),
        "config_sha256": _sha256(vendor / "config_n0.dat"),
    }


def _native_score(
    source: Path,
    output: Path,
    vendor: Path,
    seed: int,
    timeout: float,
    provenance: dict[str, Any] | None = None,
) -> SimRNAResult:
    executable = vendor / "SimRNA"
    config = vendor / "config_n0.dat"
    if not executable.is_file() or not config.is_file() or not (vendor / "data").is_dir():
        raise FileNotFoundError("simrna_dir must contain SimRNA, config_n0.dat and data/.")
    with TemporaryDirectory(prefix="simrna-") as directory:
        work = Path(directory)
        (work / "data").symlink_to(vendor / "data", target_is_directory=True)
        shutil.copyfile(config, work / "config.dat")
        info = _prepare_input(source, work / "input.pdb")
        command = [
            str(executable),
            "-p",
            "input.pdb",
            "-c",
            "config.dat",
            "-n",
            "0",
            "-R",
            str(seed),
            "-o",
            "score",
        ]
        environment = dict(os.environ, OMP_NUM_THREADS="1")
        start = perf_counter()
        try:
            process = subprocess.run(
                command,
                cwd=work,
                capture_output=True,
                text=True,
                timeout=timeout,
                env=environment,
                check=False,
            )
        except subprocess.TimeoutExpired as error:
            (output / "failure.json").write_text(json.dumps({"error": "timeout"}))
            raise RuntimeError(f"SimRNA exceeded {timeout}s. See {output}") from error
        elapsed = perf_counter() - start
        (output / "stdout.log").write_text(process.stdout, encoding="utf-8")
        (output / "stderr.log").write_text(process.stderr, encoding="utf-8")
        for path in work.iterdir():
            if path.is_file():
                shutil.copyfile(path, output / path.name)
        if process.returncode:
            raise RuntimeError(
                f"SimRNA exited with {process.returncode}; see {output / 'stderr.log'}"
            )
        total, components = _parse_output(
            process.stdout + process.stderr, (output / "score.trafl").read_text()
        )
        displacement = _check_geometry(
            output / "input.pdb", output / "score-000001.pdb", output / "score.trafl"
        )
        metadata = {
            "max_bead_displacement_angstrom": displacement,
            **(provenance if provenance is not None else _engine_provenance(vendor)),
            "engine": "SimRNA",
            "version": "3.20",
            "iterations": 0,
            "seed": seed,
            "units": "SimRNA statistical score (not kcal/mol)",
            "input_sha256": _sha256(source),
            "command": command,
            "seconds_engine": elapsed,
            "total_energy": total,
            "components": components,
            **info,
        }
        run_file = output / "run.json"
        run_file.write_text(json.dumps(metadata, indent=2), encoding="utf-8")
    return SimRNAResult(total, components, output, run_file)


def _load_result(output: Path) -> SimRNAResult:
    run_file = output / "run.json"
    metadata = json.loads(run_file.read_text(encoding="utf-8"))
    return SimRNAResult(float(metadata["total_energy"]), metadata["components"], output, run_file)


def run_simrna(
    pdb_path: str | Path,
    output_dir: str | Path,
    *,
    simrna_dir: str | Path,
    docker_image: str | None = "simrna-adapter:local",
    seed: int = 1,
    timeout: float = 120,
) -> SimRNAResult:
    """Score a prepared canonical RNA PDB; never mutate or minimize it.

    Full-atom and complete five-bead inputs are accepted. Extra atoms are removed
    without moving the five retained beads. Missing beads are never reconstructed.
    Use docker_image=None only on Linux with the official dynamic executable.
    """
    source = Path(pdb_path).resolve(strict=True)
    vendor = Path(simrna_dir).resolve(strict=True)
    output = Path(output_dir).resolve()
    if seed < 1 or not math.isfinite(timeout) or timeout <= 0:
        raise ValueError("Use a positive seed and timeout.")
    if output.exists() and any(output.iterdir()):
        raise ValueError(f"Output directory must be empty: {output}")
    output.mkdir(parents=True, exist_ok=True)
    if docker_image is None:
        return _native_score(source, output, vendor, seed, timeout)
    command = [
        "docker",
        "run",
        "--rm",
        "--network",
        "none",
        "--mount",
        f"type=bind,source={source},target=/input.pdb,readonly",
        "--mount",
        f"type=bind,source={vendor},target=/simrna,readonly",
        "--mount",
        f"type=bind,source={output},target=/output",
        docker_image,
        "/input.pdb",
        "/output",
        "--simrna-dir",
        "/simrna",
        "--native",
        "--seed",
        str(seed),
        "--timeout",
        str(timeout),
    ]
    process = subprocess.run(
        command, capture_output=True, text=True, timeout=timeout + 30, check=False
    )
    (output / "docker.log").write_text(process.stdout + process.stderr, encoding="utf-8")
    if process.returncode:
        raise RuntimeError(f"SimRNA container exited with {process.returncode}; see {output}")
    return _load_result(output)


def run_simrna_batch(
    structures: list[tuple[str, Path]],
    output_dir: str | Path,
    *,
    simrna_dir: str | Path,
    docker_image: str | None = "simrna-adapter:local",
    timeout_per_structure: float = 120,
) -> dict[str, SimRNAResult]:
    """Score already prepared PDBs sequentially. This is not incremental scoring."""
    names = [name for name, _ in structures]
    if len(set(names)) != len(names) or any(
        not re.fullmatch(r"[A-Za-z0-9_-]+", name) for name in names
    ):
        raise ValueError("Use unique simple names (letters, digits, underscore, hyphen).")
    root = Path(output_dir)
    results = {}
    for name, source in structures:
        results[name] = run_simrna(
            source,
            root / name,
            simrna_dir=simrna_dir,
            docker_image=docker_image,
            timeout=timeout_per_structure,
        )
    return results
