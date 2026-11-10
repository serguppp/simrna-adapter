# Copy into your application only if you want a small wrapper around the adapter.
from pathlib import Path

from simrna_adapter import SimRNAResult, run_simrna

SIMRNA_DIR = Path.home() / "Desktop/tools/SimRNA-3.20/installation_dynamic/SimRNA_64bitIntel_Linux"


def score_structure(pdb_path: Path, output_dir: Path) -> SimRNAResult:
    return run_simrna(pdb_path, output_dir, simrna_dir=SIMRNA_DIR)
