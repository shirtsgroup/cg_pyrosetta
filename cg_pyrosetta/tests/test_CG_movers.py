"""
Tests for the movers in cg_pyrosetta.CG_movers.
"""

import numpy as np
import pytest

import cg_pyrosetta
from cg_pyrosetta import CG_movers


def coordinates(pose):
    return np.array([[pose.xyz(cg_pyrosetta.pyrosetta.AtomID(atom, residue))[i] for i in range(3)]
                     for residue in range(1, pose.size() + 1)
                     for atom in range(1, pose.residue(residue).natoms() + 1)])


def atom_key(atom):
    return (atom.rsd(), atom.atomno())


@pytest.mark.parametrize("mover_class", [
    CG_movers.CGSmallMover,
    CG_movers.CGSmallAngleMover,
    CG_movers.CGBondLengthMover,
    CG_movers.CGBBSmallMover,
    CG_movers.CGSCSmallMover,
    CG_movers.CGBBSmallAngleMover,
    CG_movers.CGSCSmallAngleMover,
])
def test_perturbation_movers_move_atoms(pose, mover_class):
    """
    Each perturbation mover must change the coordinates of the pose
    """
    old = coordinates(pose)
    mover = mover_class(pose)
    for _ in range(50):
        mover.apply(pose)
    assert not np.allclose(old, coordinates(pose))


def test_internal_coordinates_are_found(pose):
    assert len(CG_movers.CGSmallMover(pose).torsions) > 0
    assert len(CG_movers.CGSmallAngleMover(pose).bond_angles) > 0
    assert len(CG_movers.CGBondLengthMover(pose).bond_lengths) > 0


def test_bond_length_mover_changes_one_bond_length(pose):
    mover = CG_movers.CGBondLengthMover(pose, d_bond=0.5)
    conf = pose.conformation()
    old = [conf.bond_length(*bond) for bond in mover.bond_lengths]
    mover.apply(pose)
    new = [conf.bond_length(*bond) for bond in mover.bond_lengths]
    changed = [not np.isclose(a, b) for a, b in zip(old, new)]
    assert sum(changed) == 1
    assert max(abs(a - b) for a, b in zip(old, new)) <= 0.25


@pytest.mark.parametrize("all_class, bb_class, sc_class, attribute", [
    (CG_movers.CGSmallMover, CG_movers.CGBBSmallMover, CG_movers.CGSCSmallMover, "torsions"),
    (CG_movers.CGSmallAngleMover, CG_movers.CGBBSmallAngleMover, CG_movers.CGSCSmallAngleMover, "bond_angles"),
])
def test_backbone_and_sidechain_movers_partition_the_move_set(pose, all_class, bb_class, sc_class, attribute):
    """
    Backbone movers hold only backbone atoms, sidechain movers hold at least
    one sidechain atom, and together they cover the unrestricted move set
    """
    conf = pose.conformation()
    everything = getattr(all_class(pose), attribute)
    backbone = getattr(bb_class(pose), attribute)
    sidechain = getattr(sc_class(pose), attribute)

    assert len(backbone) > 0
    assert len(sidechain) > 0
    assert len(backbone) + len(sidechain) == len(everything)
    for atoms in backbone:
        assert all(conf.atom_is_backbone_norefold(atom.rsd(), atom.atomno()) for atom in atoms)
    for atoms in sidechain:
        assert not all(conf.atom_is_backbone_norefold(atom.rsd(), atom.atomno()) for atom in atoms)


def test_randomizeTorsions(pose):
    conf = pose.conformation()
    randomize = CG_movers.randomizeTorsions(pose)
    old = [conf.torsion_angle(*atoms) for atoms in randomize.torsions]
    randomize.apply(pose)
    new = [conf.torsion_angle(*atoms) for atoms in randomize.torsions]
    assert not np.allclose(old, new)


def test_randomizeAngles(pose):
    conf = pose.conformation()
    randomize = CG_movers.randomizeAngles(pose)
    old = [conf.bond_angle(*atoms) for atoms in randomize.bond_angles]
    randomize.apply(pose)
    new = [conf.bond_angle(*atoms) for atoms in randomize.bond_angles]
    assert not np.allclose(old, new)


def test_setBondLengths(pose):
    conf = pose.conformation()
    set_lengths = CG_movers.setBondLengths(pose, {"BB1 SC1": 2.5, "BB2 BB3": 1.7})
    set_lengths.apply(pose)

    checked = {"BB1 SC1": 0, "BB2 BB3": 0}
    for atoms, name in zip(set_lengths.bond_lengths, set_lengths.bond_names):
        reverse_name = " ".join(reversed(name.split(" ")))
        if name == "BB1 SC1" or reverse_name == "BB1 SC1":
            assert conf.bond_length(*atoms) == pytest.approx(2.5)
            checked["BB1 SC1"] += 1
        elif name == "BB2 BB3" or reverse_name == "BB2 BB3":
            assert conf.bond_length(*atoms) == pytest.approx(1.7)
            checked["BB2 BB3"] += 1
        else:
            assert conf.bond_length(*atoms) == pytest.approx(1.0)
    assert all(count > 0 for count in checked.values())


def test_setBondAngle_degrees(pose):
    conf = pose.conformation()
    angles = {"BB1 BB2 BB3": 100}
    set_angle = CG_movers.setBondAngle(pose, angles)
    set_angle.apply(pose)

    # the caller's dictionary must not be converted in place
    assert angles == {"BB1 BB2 BB3": 100}

    n_checked = 0
    for atoms, name in zip(set_angle.bond_angles, set_angle.angle_names):
        if name in ("BB1 BB2 BB3", "BB3 BB2 BB1"):
            assert conf.bond_angle(*atoms) == pytest.approx(100 * np.pi / 180)
            n_checked += 1
    assert n_checked > 0


def test_setBondAngle_radians(pose):
    conf = pose.conformation()
    set_angle = CG_movers.setBondAngle(pose, {"BB1 BB2 BB3": 2.0}, radians=True)
    set_angle.apply(pose)
    for atoms, name in zip(set_angle.bond_angles, set_angle.angle_names):
        if name in ("BB1 BB2 BB3", "BB3 BB2 BB1"):
            assert conf.bond_angle(*atoms) == pytest.approx(2.0)


def test_setTorsions(pose):
    """
    Each named torsion is set to the requested value. Torsions that share a
    central bond are coupled, so only one torsion name is set here.
    """
    conf = pose.conformation()
    set_torsions = CG_movers.setTorsions(pose, {"BB1 BB2 BB3 BB1": 60})
    set_torsions.apply(pose)

    n_checked = 0
    for atoms, name in zip(set_torsions.torsions, set_torsions.torsion_names):
        if name == "BB1 BB2 BB3 BB1":
            assert conf.torsion_angle(*atoms) == pytest.approx(60 * np.pi / 180)
            n_checked += 1
    assert n_checked > 0


def test_setTorsions_reversed_name(pose):
    conf = pose.conformation()
    set_torsions = CG_movers.setTorsions(pose, {"BB1 BB3 BB2 BB1": -45})
    set_torsions.apply(pose)

    n_checked = 0
    for atoms, name in zip(set_torsions.torsions, set_torsions.torsion_names):
        if name == "BB1 BB2 BB3 BB1":
            assert conf.torsion_angle(*atoms) == pytest.approx(-45 * np.pi / 180)
            n_checked += 1
    assert n_checked > 0


def test_setTorsions_ignores_unknown_names(pose):
    old = coordinates(pose)
    CG_movers.setTorsions(pose, {"NOT A REAL TORSION": 60}).apply(pose)
    assert np.allclose(old, coordinates(pose))
