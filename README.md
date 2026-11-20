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

## Mutational frustration geometry

Protocol version 3 fixes the native N1/N9 and sugar coordinates. A rigid template
rotation aligns the N-to-C1' direction and the glycosidic chi-defining base plane
(C4 for purines, C2 for pyrimidines). This preserves native chi using the native
O4'-C1'-N reference. All base distances within the template are retained.
Same-identity substitutions keep the input coordinates exactly.

Scans now require full-atom reference and template PDBs containing C1' and the
chi-defining base atom. Conversion to five-bead input happens after mutation.
Contact selection still uses native base midpoints (default 8 angstrom).
Run metadata identifies protocol 3. Earlier scans used different mutation geometry
and should be rerun for comparisons under this protocol.

## Pair and total frustration

`run_simrna_frustration` now computes both strategies by default. Set
`pair_energies=False` to retain the previous total-only runtime.
`contact_scores.csv` has `F_pair`, `F_total`, their native energies, means and
population SDs. `contact_variants.csv` retains total components and adds
`native_E_pair`, `E_pair`, `delta_E_pair`, `native_E_total`.
The legacy `F` column in `frustration.csv` still means `F_total`.

For SimRNA, E_pair means nonlocal base-base plus base-backbone interactions
between the selected residues (both directions). Short stacking is already in
base-base. Local backbone strain is excluded. To obtain it from the licensed
closed engine, all other residues are rigidly translated apart, preserving the
full residue sequence, chain membership and adjacency. The chosen pair stays
at its variant coordinates. E_total minus local geometry from this diagnostic
run gives the nonlocal pair interaction. Do not interpret the artificial masked
structure's total energy. No atoms in the actual scored mutants are relaxed.

This needs 16 additional engine calls per contact, inside the same container.
Pair energies have about 2e-6 numerical uncertainty from the six-decimal log
output. F_pair is undefined (empty) when population SD is at most 2e-6.
Native reconstruction of all pair energies is regression-tested against the
full structure, allowing for printed rounding. Input topology and mutation
geometry remain the same for both strategies. This is an adapter decomposition,
not a native SimRNA per-contact export. Rebuild the Docker image after updating.

## Shared mutation targets

`select_simrna_contacts(pdb_path, cutoff=8.0)` returns the adapter's native
contact list without running the licensed engine. It uses the same selector
as a full scan: midpoint (N9+C6)/2 for A/G or (N1+C4)/2 for C/U.
The cutoff selects pairs to mutate; it does not truncate SimRNA total energy.
To compare models on exactly the same targets, pass the chosen list via
`contacts=...` and the same four-ID CSV to cgRNASP `--contacts`.


## Optional native-ring mutation geometry

Pass `preserve_same_class_ring=True` to `run_simrna_frustration` to keep the
observed base-ring coordinates for A/G and C/U exchanges. All three base beads
then stay exactly fixed for these exchanges. Cross-class exchanges still use
the N1/N9-anchored, glycosidic-chi-preserving template placement.

The default is `False`, preserving the existing protocol. The opt-in policy is
recorded as mutation protocol 4 in `run.json`; the existing policy remains 3.
This is a geometry-sensitivity option, not energy minimization or a guarantee
of clash-free mutants. Target-specific exocyclic groups and SimRNA's base shape
still change. It adds no SimRNA engine calls.

Rebuild `simrna-adapter:local` after updating the source. The Python API rejects
an old image that silently ignores the requested native-ring policy.
