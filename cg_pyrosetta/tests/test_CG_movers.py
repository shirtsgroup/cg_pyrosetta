"""
Unit and regression test for the cg_pyrosetta package.
"""

import cg_pyrosetta
import pyrosetta
import pytest
import sys

@pytest.fixture()
def pose():
    return(pyrosetta.pose_from_sequence('X[CG11:CGLower]X[CG11]X[CG11]X[CG11]X[CG11]X[CG11]X[CG11:CGUpper]'))

def testCGSmallMover(pose):
    """
    Testing CGSmallMover works
    """
    x_old = pose.xyz(pyrosetta.AtomID(1,7)).x
    y_old = pose.xyz(pyrosetta.AtomID(1,7)).y
    z_old = pose.xyz(pyrosetta.AtomID(1,7)).z

    small_mover = cg_pyrosetta.CG_movers.CGSmallMover(pose)
    
    for _ in range(200):
        small_mover.apply(pose)
    
    x_new = pose.xyz(pyrosetta.AtomID(1,7)).x
    y_new = pose.xyz(pyrosetta.AtomID(1,7)).y
    z_new = pose.xyz(pyrosetta.AtomID(1,7)).z    

    assert x_new != x_old
    assert y_new != y_old
    assert z_new != z_old

def testCGShearMover(pose):
    """
    Testing CGShearMover works
    """
    x_old = pose.xyz(pyrosetta.AtomID(1,7)).x
    y_old = pose.xyz(pyrosetta.AtomID(1,7)).y
    z_old = pose.xyz(pyrosetta.AtomID(1,7)).z

    small_mover = cg_pyrosetta.CG_movers.CGShearMover(pose)
    
    for _ in range(200):
        small_mover.apply(pose)
    
    x_new = pose.xyz(pyrosetta.AtomID(1,7)).x
    y_new = pose.xyz(pyrosetta.AtomID(1,7)).y
    z_new = pose.xyz(pyrosetta.AtomID(1,7)).z    

    assert x_new != x_old
    assert y_new != y_old
    assert z_new != z_old


def testCGSmallAngleMover(pose):
    x_old = pose.xyz(pyrosetta.AtomID(1,7)).x
    y_old = pose.xyz(pyrosetta.AtomID(1,7)).y
    z_old = pose.xyz(pyrosetta.AtomID(1,7)).z

    small_angle_mover = cg_pyrosetta.CG_movers.CGSmallAngleMover(pose)

    for _ in range(200):
        small_angle_mover.apply(pose)
    
    x_new = pose.xyz(pyrosetta.AtomID(1,7)).x
    y_new = pose.xyz(pyrosetta.AtomID(1,7)).y
    z_new = pose.xyz(pyrosetta.AtomID(1,7)).z    

    assert x_new != x_old
    assert y_new != y_old
    assert z_new != z_old

def test_randomizeBB(pose):
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

def test_randomizeBBAngles(pose):
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

def test_setBB():
    pass

def test_setBBAngles():
    pass
