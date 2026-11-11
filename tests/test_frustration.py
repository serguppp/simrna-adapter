import gemmi
import pytest

from simrna_adapter.frustration import _base_atoms, _frame, _replacement, _state


def residue(name, points):
    result = gemmi.Residue()
    result.name = name
    result.seqid = gemmi.SeqId(1, " ")
    names = ("N9", "C2", "C6") if name == "G" else ("N1", "C2", "C4")
    for label, xyz in zip(("P", "C4'", *names), points, strict=True):
        atom = gemmi.Atom()
        atom.name = label
        atom.pos = gemmi.Position(*xyz)
        result.add_atom(atom)
    return result


def test_template_preserves_backbone_and_base_frame():
    native = residue("G", [(0, 0, 0), (1, 0, 0), (1, 1, 1), (2, 1, 0), (1, 3, 1)])
    template = residue("C", [(0, 0, 0), (1, 0, 0), (0, 0, 0), (1, 1, 0), (3, 0, 0)])
    mutant = _replacement(native, template)
    assert mutant.name == "C"
    for index in (0, 1):
        assert mutant[index].pos.dist(native[index].pos) == 0
    before, before_axes = _frame(_base_atoms(native))
    after, after_axes = _frame(_base_atoms(mutant))
    assert before.dist(after) == pytest.approx(0)
    for a, b in zip(before_axes, after_axes, strict=True):
        assert a.dot(b) == pytest.approx(1)
    for a, b in zip(_base_atoms(mutant), _base_atoms(template), strict=True):
        origin, _ = _frame(_base_atoms(template))
        assert a.pos.dist(after) == pytest.approx(b.pos.dist(origin))


def test_reuses_native_and_single_mutations_across_contacts():
    native = ["A", "G", "U"]
    assert _state(0, 1, "A", "G", native) == ()
    assert _state(0, 1, "C", "G", native) == _state(0, 2, "C", "U", native) == ((0, "C"),)
    assert _state(0, 1, "C", "A", native) == _state(1, 0, "A", "C", native)
