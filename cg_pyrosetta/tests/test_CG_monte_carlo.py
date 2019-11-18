import numpy as np
import pyrosetta
import pytest
from cg_pyrosetta.CG_monte_carlo import CGMonteCarlo, EnergyFunctionFactory
import os
import sys
import warnings

current_path = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.abspath(current_path + '/../PyRosetta4.modified'))


@pytest.fixture
def pose():
    return(pyrosetta.pose_from_sequence('AAAAAA'))

@pytest.fixture
def e_function_factory():
    return(EnergyFunctionFactory())

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
    assert np.isclose(cg_monte_carlo.get_energy(), 5.416006218055426)


def test_run(cg_monte_carlo):
    cg_monte_carlo.run()
    assert(cg_monte_carlo.mc.total_trials() == 30)

def test_run_output(cg_monte_carlo):
    cg_monte_carlo._output = True
    cg_monte_carlo.run()
    


def test_energy_function_factory(e_function_factory, pose):
    sf = e_function_factory.build_energy_function(["fa_atr", "fa_rep"], [1, 1])
    assert np.isclose(sf(pose), 5.416006218055426)

def test_energy_function_factory_invalid_e_term(e_function_factory):
    with pytest.warns(UserWarning):
        e_function_factory.build_energy_function(["not_score_term"], [1])