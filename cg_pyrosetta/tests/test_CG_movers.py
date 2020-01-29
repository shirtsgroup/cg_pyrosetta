"""
Unit and regression test for the cg_pyrosetta package.
"""

import cg_pyrosetta
import pyrosetta
import numpy as np
import pytest
import sys
import pytest

pyrosetta.init()

# @pytest.fixture()
# def pose():
#     return(pyrosetta.pose_from_sequence('X[CG11x3:CGLower]X[CG11x3]X[CG11x3][CG11x3:CGUpper]'))


def testCGSmallMover():
    """
    Testing CGSmallMover works
    """
    pose = pyrosetta.pose_from_sequence('X[CG11x3:CGLower]X[CG11x3]X[CG11x3][CG11x3:CGUpper]')
    x_old = pose.xyz(pyrosetta.AtomID(1, pose.size())).x
    y_old = pose.xyz(pyrosetta.AtomID(1, pose.size())).y
    z_old = pose.xyz(pyrosetta.AtomID(1, pose.size())).z

    small_mover = cg_pyrosetta.CG_movers.CGSmallMover(pose)

    for _ in range(200):
        small_mover.apply(pose)

    x_new = pose.xyz(pyrosetta.AtomID(1, pose.size())).x
    y_new = pose.xyz(pyrosetta.AtomID(1, pose.size())).y
    z_new = pose.xyz(pyrosetta.AtomID(1, pose.size())).z

    assert x_new != x_old
    assert y_new != y_old
    assert z_new != z_old


def testCGShearMover():
    """
    Testing CGShearMover works
    """
    pose = pyrosetta.pose_from_sequence('X[CG11x3:CGLower]X[CG11x3]X[CG11x3][CG11x3:CGUpper]')

    x_old = pose.xyz(pyrosetta.AtomID(1, pose.size())).x
    y_old = pose.xyz(pyrosetta.AtomID(1, pose.size())).y
    z_old = pose.xyz(pyrosetta.AtomID(1, pose.size())).z

    small_mover = cg_pyrosetta.CG_movers.CGShearMover(pose)

    for _ in range(200):
        small_mover.apply(pose)

    x_new = pose.xyz(pyrosetta.AtomID(1, pose.size())).x
    y_new = pose.xyz(pyrosetta.AtomID(1, pose.size())).y
    z_new = pose.xyz(pyrosetta.AtomID(1, pose.size())).z

    assert x_new != x_old
    assert y_new != y_old
    assert z_new != z_old


def test_cgsmallscmover():
    pass
    # assert False


def testCGSmallAngleMover():
    pose = pyrosetta.pose_from_sequence('X[CG11x3:CGLower]X[CG11x3]X[CG11x3][CG11x3:CGUpper]')

    x_old = pose.xyz(pyrosetta.AtomID(1, pose.size())).x
    y_old = pose.xyz(pyrosetta.AtomID(1, pose.size())).y
    z_old = pose.xyz(pyrosetta.AtomID(1, pose.size())).z

    small_angle_mover = cg_pyrosetta.CG_movers.CGSmallAngleMover(pose)

    for _ in range(200):
        small_angle_mover.apply(pose)

    x_new = pose.xyz(pyrosetta.AtomID(1, pose.size())).x
    y_new = pose.xyz(pyrosetta.AtomID(1, pose.size())).y
    z_new = pose.xyz(pyrosetta.AtomID(1, pose.size())).z

    assert x_new != x_old
    assert y_new != y_old
    assert z_new != z_old


def test_randomizeBB():
    pose = pyrosetta.pose_from_sequence('X[CG11x3:CGLower]X[CG11x3]X[CG11x3][CG11x3:CGUpper]')

    small = cg_pyrosetta.CG_movers.CGSmallMover(pose)
    old_angles = []
    for atoms in small.dihes:
        old_angles.append(pose.conformation().torsion_angle(atoms[0], atoms[1], atoms[2], atoms[3]))

    randomize = cg_pyrosetta.CG_movers.randomizeBackBone(pose)
    randomize.apply(pose)

    new_angles = []
    for atoms in small.dihes:
        new_angles.append(pose.conformation().torsion_angle(atoms[0], atoms[1], atoms[2], atoms[3]))

    bool_array = [old_angles[i] != new_angles[i] for i in range(len(small.dihes))]
    assert(all(bool_array))


def test_randomizeBBAngles():
    pose = pyrosetta.pose_from_sequence('X[CG11x3:CGLower]X[CG11x3]X[CG11x3][CG11x3:CGUpper]')

    small_angle = cg_pyrosetta.CG_movers.CGSmallAngleMover(pose)
    old_angles = []
    for atoms in small_angle.angles:
        old_angles.append(pose.conformation().bond_angle(atoms[0], atoms[1], atoms[2]))

    randomize = cg_pyrosetta.CG_movers.randomizeBackBoneAngles(pose)
    randomize.apply(pose)

    new_angles = []
    for atoms in small_angle.angles:
        new_angles.append(pose.conformation().bond_angle(atoms[0], atoms[1], atoms[2]))

    bool_array = [old_angles[i] != new_angles[i] for i in range(len(small_angle.angles))]
    assert(all(bool_array))


def test_setBBBL():
    pose = pyrosetta.pose_from_sequence('X[CG11x3:CGLower]X[CG11x3]X[CG11x3][CG11x3:CGUpper]')
    set_bbbl = cg_pyrosetta.CG_movers.setBackBoneBondLengths(pose,
                                                             {
                                                              "BB1 BB2": 20.3,
                                                              "BB2 BB3": 20.3,
                                                              "BB1 BB3": 20.3,
                                                             })
    set_bbbl.apply(pose)
    for atom1, atom2 in set_bbbl.bb_bonds[1:-1]:
        assert pose.conformation().bond_length(atom1, atom2) == 20.3


def test_setBBBL_variable_bb():
    pose = cg_pyrosetta.pyrosetta.pose_from_sequence('X[CG31]'*10)
    set_bbbl = cg_pyrosetta.CG_movers.setBackBoneBondLengths(pose, {"BB1 BB2": 2.3, "BB3 BB1": 3.3})
    set_bbbl.apply(pose)
    assert pose.conformation().bond_length(cg_pyrosetta.pyrosetta.AtomID(2, 1),
                                           cg_pyrosetta.pyrosetta.AtomID(3, 1),
                                           ) == pytest.approx(1.0)


def test_setBBAngles():
    pass

def test_newCGSmallAngleMover():
    pymol = pyrosetta.PyMOLMover()
    pose = pyrosetta.pose_from_sequence('X[CG11x3:CGLower]X[CG11x3]X[CG11x3][CG11x3:CGUpper]')
    conf = pose.conformation()
    small_mover = cg_pyrosetta.CG_movers.newCGSmallAngleMover(pose)
    # print(small_mover.bond_angles)
    assert len(small_mover.bond_angles) != 0
    pymol.apply(pose)
    for angle in small_mover.bond_angles:
        print("Changing Angle:", angle[0], angle[1], angle[2])
        old_angle = conf.bond_angle(angle[0], angle[1], angle[2])
        conf.set_bond_angle(
                            angle[0],
                            angle[1],
                            angle[2], 
                            old_angle + 10 * np.pi / 180
                            )
        pymol.apply(pose)
        assert conf.bond_angle(angle[0], angle[1], angle[2]) == pytest.approx(old_angle + 10 * np.pi / 180)
    # assert False