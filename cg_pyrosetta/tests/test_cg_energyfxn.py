import cg_pyrosetta
import pyrosetta
import math

def testMMTwist():
    #Set CG1-CG1-CG1-CG1 torsion
    cg_pyrosetta.change_parameters.changeTorsionParameters({'CG2 CG1 CG1 CG2':[10, 6, 35]})
    pyrosetta.init()
    # Build 1-1 CG model
    pose = pyrosetta.pose_from_sequence('X[CG11]X[CG11]X[CG11]X[CG11]X[CG11]X[CG11]X[CG11]X[CG11]')
    
    # Randomize torsion in backbone
    randomize = cg_pyrosetta.CG_movers.randomizeBackBone(pose)
    randomize.apply(pose)

    # Scorefunction w/o mm_twist
    scorefunction = pyrosetta.ScoreFunction()
    scorefunction.set_weight(pyrosetta.rosetta.core.scoring.fa_atr, 1)
    scorefunction.set_weight(pyrosetta.rosetta.core.scoring.fa_rep, 1)
    
    # Scorefunction w/ mm_twist
    scorefunction_mm = pyrosetta.ScoreFunction()
    scorefunction_mm.set_weight(pyrosetta.rosetta.core.scoring.fa_atr, 1)
    scorefunction_mm.set_weight(pyrosetta.rosetta.core.scoring.fa_rep, 1)
    scorefunction_mm.set_weight(pyrosetta.rosetta.core.scoring.mm_twist, 1)
    print(scorefunction(pose))
    print(scorefunction_mm(pose))
    assert not math.isclose(scorefunction(pose), scorefunction_mm(pose))

    


def testMMBend():
    #Set CG1-CG1-CG1-CG1 torsion
    cg_pyrosetta.change_parameters.changeAngleParameters({'CG2 CG1 CG1':[10, 50]})
    pyrosetta.init()
    # Build 1-1 CG model
    pose = pyrosetta.pose_from_sequence('X[CG11]X[CG11]X[CG11]X[CG11]X[CG11]X[CG11]X[CG11]X[CG11]')

    # Randomize torsion in backbone
    randomize = cg_pyrosetta.CG_movers.randomizeBackBone(pose)
    randomize.apply(pose)

    # Scorefunction w/o mm_twist
    scorefunction = pyrosetta.ScoreFunction()
    scorefunction.set_weight(pyrosetta.rosetta.core.scoring.fa_atr, 1)
    scorefunction.set_weight(pyrosetta.rosetta.core.scoring.fa_rep, 1)

    # Scorefunction w/ mm_twist
    scorefunction_mm = pyrosetta.ScoreFunction()
    scorefunction_mm.set_weight(pyrosetta.rosetta.core.scoring.fa_atr, 1)
    scorefunction_mm.set_weight(pyrosetta.rosetta.core.scoring.fa_rep, 1)
    scorefunction_mm.set_weight(pyrosetta.rosetta.core.scoring.mm_bend, 1)
    print(scorefunction(pose))
    print(scorefunction_mm(pose))
    assert not math.isclose(scorefunction(pose), scorefunction_mm(pose))
