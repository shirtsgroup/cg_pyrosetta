import numpy as np
import pyrosetta
import pytest
from cg_pyrosetta.CG_monte_carlo import CGMonteCarlo
import os
import sys
import warnings

current_path = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.abspath(current_path + '/../PyRosetta4.modified'))


@pytest.fixture
def pose():
    return(pyrosetta.pose_from_sequence('AAAAAA'))


@pytest.fixture
def score_function():
    score = pyrosetta.ScoreFunction()
    score.set_weight(pyrosetta.rosetta.core.scoring.fa_atr, 1)
    score.set_weight(pyrosetta.rosetta.core.scoring.fa_rep, 1)
    return(score)


@pytest.fixture
def seq_mover():
    small = pyrosetta.rosetta.protocols.simple_moves.SmallMover()
    rep_small = pyrosetta.RepeatMover(small, 20)
    return(rep_small)


@pytest.fixture
def cg_monte_carlo(pose, score_function, seq_mover):
    mc_obj = CGMonteCarlo(
        pose=pose,
        score=score_function,
        seq_mover=seq_mover,
        n_steps=30,
        output=False,)
    return(mc_obj)


def test_kT_property(cg_monte_carlo):
    cg_monte_carlo.kT = 30
    assert(cg_monte_carlo._kT == 30)
    assert(cg_monte_carlo.mc.temperature() == 30)
    assert(cg_monte_carlo.mc_trial.mc().temperature() == 30)


def test_get_energy(cg_monte_carlo):
    # fa_atr + fa_rep score of alanine hexamer
    assert(np.isclose(cg_monte_carlo.get_energy(), 5.416006218055426))


def test_run(cg_monte_carlo):
    cg_monte_carlo.run()
    assert(cg_monte_carlo.mc.total_trials() == 30)
