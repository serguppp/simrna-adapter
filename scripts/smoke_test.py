"""Run manually with licensed SimRNA installed locally. Does not download the engine."""

import argparse
import json
from pathlib import Path
from time import perf_counter

from simrna_adapter import run_simrna


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("pdb", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--simrna-dir", required=True, type=Path)
    args = parser.parse_args()
    rows = []
    for number, seed in enumerate((1, 1, 2), start=1):
        start = perf_counter()
        result = run_simrna(
            args.pdb, args.output / f"repeat_{number}", simrna_dir=args.simrna_dir, seed=seed
        )
        metadata = json.loads(result.run_file.read_text())
        rows.append(
            {
                "seed": seed,
                "E_total": result.total_energy,
                "wall_seconds": perf_counter() - start,
                "engine_seconds": metadata["seconds_engine"],
                "max_bead_displacement_angstrom": metadata["max_bead_displacement_angstrom"],
            }
        )
    if len({row["E_total"] for row in rows}) != 1:
        raise RuntimeError("Zero-step energy unexpectedly depends on repetition/seed.")
    (args.output / "smoke_test.json").write_text(json.dumps(rows, indent=2))
    print(json.dumps(rows, indent=2))


if __name__ == "__main__":
    main()
