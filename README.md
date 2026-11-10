# SimRNA adapter

Python wrapper for **SimRNA 3.20 zero-iteration energy scoring** (`-n 0`).
This repository contains our adapter, not a fork or copy of the SimRNA engine.
No SimRNA executable, energy tables, or licensed documents are included.

## Windows setup

Install the official **64-bit Linux dynamic build**, after accepting its academic
license. Keep it outside this repository. In this installation it is located at:

`C:/Users/sergi/Desktop/tools/SimRNA-3.20/installation_dynamic/SimRNA_64bitIntel_Linux`

Build the adapter image (it installs libgomp1; the SimRNA engine is mounted at runtime):

```powershell
docker build -t simrna-adapter:local .
uv sync --extra dev
```

The image contains only our adapter and its open-source dependencies. The official
static build crashed during our smoke test; use the dynamic build listed above.

## Add to the RNA-frustratometer package

Run this yourself from `C:/Users/sergi/Desktop/rna-frustratometer`:

```powershell
uv add --editable "C:/Users/sergi/Desktop/tools/simrna-adapter"
```

The package's source and dependency files have not been changed by this setup.

## Python

```python
from pathlib import Path
from simrna_adapter import run_simrna

result = run_simrna(
    Path("reference.pdb"),
    Path("results/simrna/reference"),
    simrna_dir=Path.home()
    / "Desktop/tools/SimRNA-3.20/installation_dynamic/SimRNA_64bitIntel_Linux",
)
print(result.total_energy)
print(result.components)
```

Use a new or empty output directory for every call. `run.json` contains the score,
three main energy contributions, a stacking subset, timings and file checksums.
Output geometry and the single-frame trajectory are checked against the input.
The stacking subset belongs to base-base energy: **do not add it again**.
Logs, the one-frame trajectory and the engine's reduced PDB are also retained.
On Linux, pass `docker_image=None` to invoke the executable directly.

CLI: `uv run simrna-score reference.pdb results/simrna/reference --simrna-dir PATH`.

## Inputs and interpretation

Provide one model with canonical A/C/G/U residues, without waters/ligands/proteins.
Each nucleotide must have P, C4' and three base atoms: N9/C2/C6 for A/G,
N1/C2/C4 for C/U. Full-atom inputs are reduced to those five atoms without moving
coordinates. Missing phosphates, modified bases and alternate conformations must
be handled explicitly beforehand. cgRNASP three-bead inputs are insufficient.

This API scores supplied structures. It does not generate mutations or relax them,
and provides no per-contact energy decomposition. `run_simrna_batch` runs prepared
PDBs sequentially; it still starts the engine for each PDB. It is a correctness
baseline, not the optimized single-process cgRNASP implementation.

The three main contributions are base-base, base-backbone and local geometry.
Scores are statistical model units, not kcal/mol or experimental free energies.
Do not mix SimRNA energies with BRiQ or cgRNASP values. Define and validate a separate
five-bead mutation protocol before interpreting a 16-variant frustration index.

## License

Our code is MIT-licensed. SimRNA has a **separate academic license**. Acceptance does
not make its engine or tables open source. Keep those files local, cite the SimRNA
paper, and check the vendor's terms before deploying a public/commercial service.
The license restricts redistribution outside the licensed institution and requires
separate authorization for commercial uses/services. See the official terms:

- [Download and license](https://iimcb.genesilico.pl/index.php/software/simrna)
- [Academic license](https://iimcb.genesilico.pl/DOWNLOAD/simrna/version_3.20/SimRNA_academic_license.pdf)
- [SimRNA paper](https://doi.org/10.1093/nar/gkv1479)

## Checks

`uv run pytest`, `uv run ruff check .`, `uv run ruff format --check .`,
`uv run mypy src`. Licensed-engine smoke tests are recorded separately in the
RNA-frustratometer package's results directory, outside this repository.
