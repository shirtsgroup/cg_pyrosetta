"""
Unit and regression test for the cg_pyrosetta package.
"""

import cg_pyrosetta
import pyrosetta
import pytest
import sys

def testCGSmallMover():
    """
    Testing CGSmallMover works
    """
    pose = pyrosetta.pose_from_sequence('X[CG11]X[CG11]X[CG11]X[CG11]X[CG11]X[CG11]X[CG11]')
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

def testCGShearMover():
    """
    Testing CGShearMover works
    """
    pose = pyrosetta.pose_from_sequence('X[CG11]X[CG11]X[CG11]X[CG11]X[CG11]X[CG11]X[CG11]')
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


def testCGSmallAngleMover():
    pose = pyrosetta.pose_from_sequence('X[CG11]X[CG11]X[CG11]X[CG11]X[CG11]X[CG11]X[CG11]')
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

