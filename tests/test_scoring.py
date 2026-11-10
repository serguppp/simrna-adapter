import subprocess
from pathlib import Path

import gemmi
import pytest

from simrna_adapter import run_simrna, run_simrna_batch
from simrna_adapter.scoring import (
    _beads,
    _check_geometry,
    _native_score,
    _parse_output,
    _prepare_input,
)

LOG = """Total energy of system -12.000000
Base-Base interactions energy -8.000
short stacking energy -2.000
Base-Backbone interact. energy -1.000
Local geometry energy -3.000000
"""
TRAJ = "1 0 -12.000000 -12.000000 1.35\n0 1 2\n"


def make_pdb(path: Path) -> Path:
    structure = gemmi.Structure()
    model = gemmi.Model("1")
    chain = gemmi.Chain("A")
    residue = gemmi.Residue()
    residue.name = "A"
    residue.seqid = gemmi.SeqId(1, " ")
    for index, name in enumerate(("P", "C4'", "N9", "C2", "C6", "O5'")):
        atom = gemmi.Atom()
        atom.name = name
        atom.pos = gemmi.Position(index, index + 1, index + 2)
        atom.occ = 1
        residue.add_atom(atom)
    chain.add_residue(residue)
    model.add_chain(chain)
    structure.add_model(model)
    structure.write_pdb(str(path))
    return path


def test_total_and_subset_are_not_double_counted() -> None:
    total, parts = _parse_output(LOG, TRAJ)
    assert total == -12
    assert parts["E_short_stacking_subset"] == -2
    assert sum(parts[k] for k in ("E_base_base", "E_base_backbone", "E_local_geometry")) == total


@pytest.mark.parametrize(
    "trajectory",
    [
        "",
        TRAJ + TRAJ,
        TRAJ.replace("-12.000000 -12.000000", "-10 -12"),
        TRAJ.replace("-12.000000 -12.000000", "nan nan"),
        TRAJ.replace("-12.000000 -12.000000", "-14 -14"),
    ],
)
def test_rejects_incomplete_restrained_or_wrong_energies(trajectory: str) -> None:
    with pytest.raises(ValueError):
        _parse_output(LOG, trajectory)


def test_reduction_preserves_all_bead_coordinates(tmp_path: Path) -> None:
    source = make_pdb(tmp_path / "full input.pdb")
    before = _beads(gemmi.read_structure(str(source)))
    dest = tmp_path / "reduced.pdb"
    metadata = _prepare_input(source, dest)
    after_structure = gemmi.read_structure(str(dest))
    assert before == _beads(after_structure)
    assert metadata == {"beads": 5, "nucleotides": 1}
    assert len(after_structure[0][0][0]) == 5


@pytest.mark.parametrize("atom", ["P", "C6", "N9"])
def test_rejects_missing_required_beads(tmp_path: Path, atom: str) -> None:
    source = make_pdb(tmp_path / "input.pdb")
    structure = gemmi.read_structure(str(source))
    residue = structure[0][0][0]
    residue.remove_atom(atom, "\0")
    structure.write_pdb(str(source))
    with pytest.raises(ValueError, match="Expected one occupied"):
        _prepare_input(source, tmp_path / "output.pdb")


def test_rejects_ligand_instead_of_silently_dropping_it(tmp_path: Path) -> None:
    source = make_pdb(tmp_path / "input.pdb")
    structure = gemmi.read_structure(str(source))
    structure[0][0][0].name = "SAM"
    structure.write_pdb(str(source))
    with pytest.raises(ValueError, match="Noncanonical/non-RNA"):
        _prepare_input(source, tmp_path / "output.pdb")


def test_rejects_multiple_models(tmp_path: Path) -> None:
    source = make_pdb(tmp_path / "input.pdb")
    structure = gemmi.read_structure(str(source))
    structure.add_model(structure[0].clone())
    with pytest.raises(ValueError, match="exactly one"):
        _beads(structure)


def test_crash_is_not_a_valid_score(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    vendor = tmp_path / "vendor"
    vendor.mkdir()
    (vendor / "SimRNA").write_text("fake")
    (vendor / "config_n0.dat").write_text("fake")
    (vendor / "data").mkdir()
    source = make_pdb(tmp_path / "input.pdb")
    output = tmp_path / "out"
    output.mkdir()
    monkeypatch.setattr(Path, "symlink_to", lambda *args, **kwargs: None)
    monkeypatch.setattr(
        subprocess, "run", lambda *args, **kwargs: subprocess.CompletedProcess([], -11, "", "crash")
    )
    with pytest.raises(RuntimeError, match="-11"):
        _native_score(source, output, vendor, 1, 10)
    assert not (output / "run.json").exists()
    assert (output / "stderr.log").read_text() == "crash"


def test_existing_output_is_never_overwritten(tmp_path: Path) -> None:
    source = make_pdb(tmp_path / "input.pdb")
    output = tmp_path / "out"
    output.mkdir()
    sentinel = output / "important.txt"
    sentinel.write_text("keep")
    with pytest.raises(ValueError, match="empty"):
        run_simrna(source, output, simrna_dir=tmp_path)
    assert sentinel.read_text() == "keep"


def test_batch_rejects_path_escape_before_execution(tmp_path: Path) -> None:
    with pytest.raises(ValueError, match="unique simple names"):
        run_simrna_batch([("../escape", tmp_path)], tmp_path, simrna_dir=tmp_path)


def test_detects_zero_step_geometry_change(tmp_path: Path) -> None:
    source = make_pdb(tmp_path / "input.pdb")
    output = tmp_path / "reduced.pdb"
    _prepare_input(source, output)
    trajectory = tmp_path / "score.trafl"
    positions = _beads(gemmi.read_structure(str(output)))
    coords = " ".join(str(value) for xyz in positions.values() for value in xyz)
    trajectory.write_text("1 1 0 0 1.35\n" + coords + "\n")
    assert _check_geometry(source, output, trajectory) == 0
    structure = gemmi.read_structure(str(output))
    structure[0][0][0][0].pos.x += 1
    structure.write_pdb(str(output))
    with pytest.raises(ValueError, match="moved"):
        _check_geometry(source, output, trajectory)
