"""
Tests confirming the molecular mechanics terms read the packaged parameters.
"""

import math

import pytest

import cg_pyrosetta


def lennard_jones_score_function(pyrosetta):
    score_function = pyrosetta.ScoreFunction()
    score_function.set_weight(pyrosetta.rosetta.core.scoring.fa_atr, 1)
    score_function.set_weight(pyrosetta.rosetta.core.scoring.fa_rep, 1)
    return score_function


@pytest.mark.parametrize("term", ["mm_twist", "mm_bend"])
def test_mm_term_contributes(pyrosetta, long_pose, term):
    """
    The packaged parameters have nonzero CG2-CG1-CG1-CG2 torsion and
    CG2-CG1-CG1 angle force constants, so each term must change the score.
    """
    randomize = cg_pyrosetta.CG_movers.randomizeTorsions(long_pose)
    randomize.apply(long_pose)

    score_function = lennard_jones_score_function(pyrosetta)
    score_function_mm = lennard_jones_score_function(pyrosetta)
    score_function_mm.set_weight(getattr(pyrosetta.rosetta.core.scoring, term), 1)

    assert not math.isclose(score_function(long_pose), score_function_mm(long_pose))
