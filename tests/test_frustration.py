import gemmi
import pytest

from simrna_adapter.frustration import _base_atoms, _frame, _glycosidic_frame, _replacement, _state


def residue(name, points):
    result = gemmi.Residue()
    result.name = name
    result.seqid = gemmi.SeqId(1, " ")
    names = ("N9", "C2", "C6") if name in "AG" else ("N1", "C2", "C4")
    for label, xyz in zip(("P", "C4'", *names), points, strict=True):
        atom = gemmi.Atom()
        atom.name = label
        atom.pos = gemmi.Position(*xyz)
        result.add_atom(atom)
    anchor = _base_atoms(result)[0].pos
    for label, vector in (
        ("C1'", gemmi.Vec3(-0.7, -0.3, 0.9)),
        ("O4'", gemmi.Vec3(-0.3, 0.3, 0.4)),
    ):
        atom = gemmi.Atom()
        atom.name = label
        xyz = gemmi.Vec3(anchor.x, anchor.y, anchor.z) + vector
        atom.pos = gemmi.Position(xyz.x, xyz.y, xyz.z)
        result.add_atom(atom)
    if name in "AG":
        atoms = _base_atoms(result)
        xyz = anchor + (atoms[1].pos - anchor) * 0.8 + (atoms[2].pos - anchor) * 0.2
        atom = gemmi.Atom()
        atom.name = "C4"
        atom.pos = xyz
        result.add_atom(atom)
    return result


@pytest.mark.parametrize("source", "ACGU")
@pytest.mark.parametrize("target", "ACGU")
def test_template_preserves_backbone_and_nitrogen(source, target):
    native = residue(source, [(0, 0, 0), (1, 0, 0), (1, 1, 1), (2, 1, 0), (1, 3, 1)])
    template = residue(target, [(0, 0, 0), (1, 0, 0), (0, 0, 0), (1, 1, 0), (3, 0, 0)])
    mutant = _replacement(native, template)
    assert mutant.name == target
    for index in (0, 1):
        assert mutant[index].pos.dist(native[index].pos) == 0
    before = _base_atoms(native)[0].pos
    after = _base_atoms(mutant)[0].pos
    assert before.dist(after) == 0
    if source == target:
        assert [tuple(a.pos) for a in mutant] == [tuple(a.pos) for a in native]
        return
    _, before_axes = _glycosidic_frame(native)
    _, after_axes = _glycosidic_frame(mutant)
    for a, b in zip(before_axes, after_axes, strict=True):
        assert a.dot(b) == pytest.approx(1)
    for a, b in zip(_base_atoms(mutant), _base_atoms(template), strict=True):
        origin = _base_atoms(template)[0].pos
        assert a.pos.dist(after) == pytest.approx(b.pos.dist(origin))


def test_reuses_native_and_single_mutations_across_contacts():
    native = ["A", "G", "U"]
    assert _state(0, 1, "A", "G", native) == ()
    assert _state(0, 1, "C", "G", native) == _state(0, 2, "C", "U", native) == ((0, "C"),)
    assert _state(0, 1, "C", "A", native) == _state(1, 0, "A", "C", native)


def test_contact_midpoint_definition_is_independent_of_mutation_anchor():
    reference = residue("G", [(0, 0, 0), (1, 0, 0), (0, 1, 0), (2, 2, 0), (4, 1, 0)])
    midpoint, _ = _frame(_base_atoms(reference))
    assert tuple(midpoint) == (2, 1, 0)


@pytest.mark.parametrize("source", "ACGU")
@pytest.mark.parametrize("target", "ACGU")
def test_glycosidic_chi_is_preserved(source, target):
    native = residue(source, [(0, 0, 0), (1, 0, 0), (1, 1, 1), (2, 1, 0), (1, 3, 1)])
    template = residue(target, [(0, 0, 0), (1, 0, 0), (0, 0, 0), (1, 1, 0), (3, 0, 0)])
    mutant = _replacement(native, template)

    def chi(r):
        labels = ("O4'", "C1'", "N9", "C4") if r.name in "AG" else ("O4'", "C1'", "N1", "C2")
        return gemmi.calculate_dihedral(*(next(a.pos for a in r if a.name == n) for n in labels))

    assert chi(mutant) == pytest.approx(chi(native))
    for name in ("P", "C4'", "C1'", "O4'"):
        assert (
            next(a.pos for a in mutant if a.name == name).dist(
                next(a.pos for a in native if a.name == name)
            )
            == 0
        )


def test_pair_isolation_preserves_sequence_topology_and_selected_coordinates():
    from simrna_adapter.frustration import _isolate_pair

    structure = gemmi.Structure()
    model = gemmi.Model("1")
    chain = gemmi.Chain("A")
    for index in range(4):
        r = residue(
            "G", [(index, 0, 0), (index, 1, 0), (index, 2, 0), (index, 3, 0), (index, 4, 1)]
        )
        r.seqid = gemmi.SeqId(index + 1, " ")
        chain.add_residue(r)
    model.add_chain(chain)
    structure.add_model(model)
    isolated = _isolate_pair(structure, (0, 2))
    assert [str(r.seqid) for r in isolated[0][0]] == ["1", "2", "3", "4"]
    for index in (0, 2):
        assert [tuple(a.pos) for a in isolated[0][0][index]] == [
            tuple(a.pos) for a in structure[0][0][index]
        ]
    for index in (1, 3):
        assert isolated[0][0][index][0].pos.dist(structure[0][0][index][0].pos) > 500
        assert isolated[0][0][index][0].pos.dist(isolated[0][0][index][1].pos) == pytest.approx(1)
    assert isolated[0][0][1][0].pos.dist(isolated[0][0][3][0].pos) > 40


def test_pair_energy_removes_local_geometry_without_double_counting_stacking():
    from pathlib import Path

    from simrna_adapter.frustration import _pair_energy
    from simrna_adapter.scoring import SimRNAResult

    result = SimRNAResult(
        995.123456,
        {
            "E_local_geometry": 1000,
            "E_base_base": -4.5,
            "E_base_backbone": -0.377,
            "E_short_stacking_subset": -2,
        },
        Path("."),
        Path("run.json"),
    )
    assert _pair_energy(result) == pytest.approx(-4.876544)


def test_pair_statistics_define_zero_variance_and_engine_rounding_as_undefined():
    from simrna_adapter.frustration import _energy_statistics

    assert _energy_statistics(2.0, [2.0] * 16)[-1] == ""
    assert _energy_statistics(2.0, [2 + 1e-7 * k for k in range(16)], resolution=2e-6)[-1] == ""
    values = list(range(16))
    from statistics import fmean, pstdev

    assert _energy_statistics(3.0, values)[-1] == pytest.approx(
        (fmean(values) - 3) / pstdev(values)
    )


def test_contact_selection_uses_base_midpoints_including_cutoff_boundary(tmp_path):
    from simrna_adapter import select_simrna_contacts

    structure = gemmi.Structure()
    model = gemmi.Model("1")
    chain = gemmi.Chain("A")
    for index, x in enumerate([0, 8, 16.001]):
        r = residue("G", [(x, 0, 0), (x, 1, 0), (x, 2, 0), (x + 1, 3, 0), (x + 2, 2, 0)])
        r.seqid = gemmi.SeqId(index + 1, " ")
        chain.add_residue(r)
    model.add_chain(chain)
    structure.add_model(model)
    path = tmp_path / "selection.pdb"
    structure.write_pdb(str(path))
    assert select_simrna_contacts(path, cutoff=8) == [("A", "1", "A", "2")]
    assert select_simrna_contacts(path, cutoff=8.01) == [("A", "1", "A", "2"), ("A", "2", "A", "3")]
    for invalid in (0, -1, float("nan"), float("inf")):
        with pytest.raises(ValueError):
            select_simrna_contacts(path, cutoff=invalid)


@pytest.mark.parametrize("source", "ACGU")
@pytest.mark.parametrize("target", "ACGU")
def test_optional_ring_policy_keeps_same_class_beads_and_cross_class_protocol(source, target):
    native = residue(source, [(0, 0, 0), (1, 0, 0), (1, 1, 1), (2, 1, 0), (1, 3, 1)])
    template = residue(target, [(0, 0, 0), (1, 0, 0), (0, 0, 0), (1, 1, 0), (3, 0, 0)])
    default = _replacement(native, template)
    fixed = _replacement(native, template, preserve_same_class_ring=True)
    if (source in "AG") == (target in "AG"):
        assert [tuple(a.pos) for a in _base_atoms(fixed)] == [
            tuple(a.pos) for a in _base_atoms(native)
        ]
    else:
        assert [(a.name, tuple(a.pos)) for a in fixed] == [(a.name, tuple(a.pos)) for a in default]
    for name in ("P", "C4'", "C1'", "O4'"):
        assert (
            next(a.pos for a in fixed if a.name == name).dist(
                next(a.pos for a in native if a.name == name)
            )
            == 0
        )
    if source == target:
        assert [(a.name, tuple(a.pos)) for a in fixed] == [(a.name, tuple(a.pos)) for a in native]


def test_requested_ring_policy_rejects_stale_docker_image(tmp_path, monkeypatch):
    import json
    import subprocess

    from simrna_adapter import run_simrna_frustration

    pdb = tmp_path / "input.pdb"
    pdb.write_text("")
    vendor = tmp_path / "vendor"
    vendor.mkdir()
    output = tmp_path / "output"

    def old_image(*args, **kwargs):
        (output / "run.json").write_text(
            json.dumps(
                {
                    "mutation_protocol_version": 3,
                    "contact_count": 1,
                    "engine_calls": 1,
                }
            )
        )
        return subprocess.CompletedProcess(args[0], 0, "", "")

    monkeypatch.setattr(subprocess, "run", old_image)
    with pytest.raises(RuntimeError, match="Rebuild"):
        run_simrna_frustration(
            pdb,
            output,
            simrna_dir=vendor,
            templates_pdb=pdb,
            preserve_same_class_ring=True,
        )
