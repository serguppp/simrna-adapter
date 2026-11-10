# Local validation ? 2026-11-10

Engine: official SimRNA 3.20, dynamic 64-bit Linux build, one OpenMP thread.
Runtime: Docker Desktop Linux, adapter image with libgomp1, no simulation (`-n 0`).

| Input | RNA size | E_total | Engine time | Entire Python/Docker call |
|---|---:|---:|---:|---:|
| CASP R1209 reference | 72 nt | -743.558261 | 0.26?0.29 s | 1.53?2.51 s |
| CASP R1205 reference | 55 nt | -565.061737 | 0.32 s | 1.61 s |
| R1209, pre-reduced to five beads | 72 nt | -743.558261 | 0.31 s | 1.56 s |

R1209 was scored three times with seeds 1, 1 and 2; energies were identical.
Maximum displacement of retained beads was 0.000 angstrom for all successful runs.
Full-atom and pre-reduced R1209 gave the same energy. TRAFL, log and component sum
checks passed. The stacking subset was not added to the three main contributions.
These are smoke-test timings on this workstation, not full-frustration timings.

The first NMR model of 1ANR was intentionally tried as another input. It was
rejected because the first nucleotide A:17 G lacks the required occupied P bead.
It was not scored, repaired or silently shortened. This is an input limitation.

The official static build segfaulted (-11) after loading tables for complete R1209
(full-atom and reduced input). The dynamic build succeeded on the same structure.
This establishes the working build for this environment, not a general diagnosis
of the static executable.

16 unit tests pass; Ruff lint, Ruff format and mypy pass. Unit tests cover missing
beads, ligands, multiple models, geometry preservation, wrong/restrained/nonfinite
energies, engine crashes, existing-output protection and batch-name validation.

Raw logs, input/output PDBs, trajectory, file hashes and timings remain outside
this repository at:
`C:/Users/sergi/Desktop/rna-frustratometer/results/simrna_integration_20261110/`.

No engine binaries, energy tables or licensed documents belong in this repository.
