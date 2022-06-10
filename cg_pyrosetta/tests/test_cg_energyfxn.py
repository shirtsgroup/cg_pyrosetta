import cg_pyrosetta
import math
import pytest

cg_pyrosetta.init()

def test_mm_twist():
    # Set CG1-CG1-CG1-CG1 torsion
    cg_pyrosetta.change_parameters.changeTorsionParameters({'CG2 CG1 CG1 CG2': [10, 6, 35]})
    # Build 1-1 CG model
    pose = cg_pyrosetta.pyrosetta.pose_from_sequence('X[CG11x3:CGLower]X[CG11x3]X[CG11x3]X[CG11x3]X[CG11x3]X[CG11x3]X[CG11x3]X[CG11x3:CGUpper]')

    # Randomize torsion in backbone
    randomize = cg_pyrosetta.CG_movers.randomizeBackBone(pose)
    randomize.apply(pose)

    # Scorefunction w/o mm_twist
    scorefunction = cg_pyrosetta.pyrosetta.ScoreFunction()
    scorefunction.set_weight(cg_pyrosetta.pyrosetta.rosetta.core.scoring.fa_atr, 1)
    scorefunction.set_weight(cg_pyrosetta.pyrosetta.rosetta.core.scoring.fa_rep, 1)

    # Scorefunction w/ mm_twist
    scorefunction_mm = cg_pyrosetta.pyrosetta.ScoreFunction()
    scorefunction_mm.set_weight(cg_pyrosetta.pyrosetta.rosetta.core.scoring.fa_atr, 1)
    scorefunction_mm.set_weight(cg_pyrosetta.pyrosetta.rosetta.core.scoring.fa_rep, 1)
    scorefunction_mm.set_weight(cg_pyrosetta.pyrosetta.rosetta.core.scoring.mm_twist, 1)
    print(scorefunction(pose))
    print(scorefunction_mm(pose))
    assert not math.isclose(scorefunction(pose), scorefunction_mm(pose))


def test_mm_bend():
    # Set CG1-CG1-CG1-CG1 torsion
    cg_pyrosetta.change_parameters.changeAngleParameters({'CG2 CG1 CG1': [10, 10]})
    # Build 1-1 CG model
    pose = cg_pyrosetta.pyrosetta.pose_from_sequence('X[CG11x3:CGLower]X[CG11x3]X[CG11x3]X[CG11x3]X[CG11x3]X[CG11x3]X[CG11x3]X[CG11x3:CGUpper]')

    # Randomize torsion in backbone
    randomize = cg_pyrosetta.CG_movers.randomizeBackBone(pose)
    randomize.apply(pose)

    # Scorefunction w/o mm_twist
    scorefunction = cg_pyrosetta.pyrosetta.ScoreFunction()
    scorefunction.set_weight(cg_pyrosetta.pyrosetta.rosetta.core.scoring.fa_atr, 1)
    scorefunction.set_weight(cg_pyrosetta.pyrosetta.rosetta.core.scoring.fa_rep, 1)

    # Scorefunction w/ mm_twist
    scorefunction_mm = cg_pyrosetta.pyrosetta.ScoreFunction()
    scorefunction_mm.set_weight(cg_pyrosetta.pyrosetta.rosetta.core.scoring.fa_atr, 1)
    scorefunction_mm.set_weight(cg_pyrosetta.pyrosetta.rosetta.core.scoring.fa_rep, 1)
    scorefunction_mm.set_weight(cg_pyrosetta.pyrosetta.rosetta.core.scoring.mm_bend, 1)
    print(scorefunction(pose))
    print(scorefunction_mm(pose))
    assert not math.isclose(scorefunction(pose), scorefunction_mm(pose))


def test_lj_terms():
    pass
