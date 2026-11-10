import argparse
from pathlib import Path

from .scoring import run_simrna


def main() -> None:
    parser = argparse.ArgumentParser(description="SimRNA energy without simulation (n=0)")
    parser.add_argument("pdb", type=Path)
    parser.add_argument("output_dir", type=Path)
    parser.add_argument("--simrna-dir", required=True, type=Path)
    parser.add_argument("--native", action="store_true", help="Linux only; skip Docker")
    parser.add_argument("--seed", type=int, default=1)
    parser.add_argument("--timeout", type=float, default=120)
    args = parser.parse_args()
    result = run_simrna(
        args.pdb,
        args.output_dir,
        simrna_dir=args.simrna_dir,
        docker_image=None if args.native else "simrna-adapter:local",
        seed=args.seed,
        timeout=args.timeout,
    )
    print(f"SimRNA energy: {result.total_energy:.6f}")
    print(f"Details: {result.run_file}")


if __name__ == "__main__":
    main()
