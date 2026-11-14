"""
Tests for the Monte Carlo objects in cg_pyrosetta.CG_monte_carlo.
"""

import os

import numpy as np
import pytest

import cg_pyrosetta
from cg_pyrosetta import CG_monte_carlo
from cg_pyrosetta.CG_monte_carlo import CGMonteCarlo, EnergyFunctionFactory, SequenceMoverFactory

SCORE_TERMS = {"fa_atr": 1, "fa_rep": 1, "fa_intra_atr": 1, "fa_intra_rep": 1, "mm_twist": 1, "mm_bend": 1}


@pytest.fixture
def score_function():
    return EnergyFunctionFactory().build_energy_function(SCORE_TERMS)


@pytest.fixture
def seq_mover(pose):
    return SequenceMoverFactory(pose).build_seq_mover({"small_dihe": 1, "small_angle": 1})


@pytest.fixture
def cg_monte_carlo(pose, score_function, seq_mover):
    return CGMonteCarlo(pose=pose, score=score_function, seq_mover=seq_mover, n_steps=30, kT=1, output=False)


def test_kT_property(cg_monte_carlo):
    cg_monte_carlo.kT = 30
    assert cg_monte_carlo._kT == 30
    assert cg_monte_carlo.mc.temperature() == 30
    assert cg_monte_carlo.mc_trial.mc().temperature() == 30


def test_get_energy(cg_monte_carlo, pose, score_function):
    assert np.isclose(cg_monte_carlo.get_energy(), score_function(pose))


def test_run(cg_monte_carlo):
    cg_monte_carlo.run()
    assert cg_monte_carlo.mc.total_trials() == 30
    assert cg_monte_carlo._total_trials == 30


def test_run_with_out_freq(pose, score_function, seq_mover):
    mc_sim = CGMonteCarlo(pose=pose, score=score_function, seq_mover=seq_mover, n_steps=30, output=False)
    mc_sim.out_freq = 10
    energy_observer = CG_monte_carlo.EnergyObserver(mc_sim, write_file=False)
    mc_sim.registerObserver(energy_observer)
    mc_sim.run()
    assert mc_sim.mc.total_trials() == 30
    assert len(energy_observer.energies) == 3


def test_minimum_energy_pose_is_not_higher_than_start(cg_monte_carlo, score_function):
    start_energy = cg_monte_carlo.get_energy()
    cg_monte_carlo.run()
    assert score_function(cg_monte_carlo.get_minimum_energy_pose()) <= start_energy + 1e-8


def test_energy_function_factory(pyrosetta, pose):
    built = EnergyFunctionFactory().build_energy_function({"fa_atr": 1, "fa_rep": 1})
    reference = pyrosetta.ScoreFunction()
    reference.set_weight(pyrosetta.rosetta.core.scoring.fa_atr, 1)
    reference.set_weight(pyrosetta.rosetta.core.scoring.fa_rep, 1)
    assert np.isclose(built(pose), reference(pose))


def test_energy_function_factory_invalid_e_term():
    with pytest.warns(UserWarning):
        EnergyFunctionFactory().build_energy_function({"not_score_term": 1})


def test_sequence_mover_factory(pose):
    seq_mover = SequenceMoverFactory(pose).build_seq_mover({"small_dihe": 10, "small_angle": 10})
    assert seq_mover.size() == 2
    assert seq_mover.apply(pose) is None


def test_sequence_mover_factory_extra_mover(pyrosetta, pose, score_function):
    mini = pyrosetta.rosetta.protocols.minimization_packing.MinMover()
    mini.min_type('lbfgs_armijo_nonmonotone')
    mini.score_function(score_function)
    seq_mover = SequenceMoverFactory(pose, {"mini": mini}).build_seq_mover({"small_dihe": 1, "mini": 1})
    assert seq_mover.size() == 2


@pytest.mark.parametrize("mover_name", ["unimplemented_mover", "shear_dihe"])
def test_sequence_mover_factory_warns_on_unimplemented_mover(pose, mover_name):
    with pytest.warns(UserWarning):
        seq_mover = SequenceMoverFactory(pose).build_seq_mover({mover_name: 10, "small_dihe": 10})
    assert seq_mover.size() == 1


def test_observers_write_output(cg_monte_carlo, tmp_path):
    energy_file = str(tmp_path / "energies.txt")
    structure_file = str(tmp_path / "structure.txt")
    pdb_file = str(tmp_path / "trajectory.pdb")
    cg_monte_carlo.out_freq = 10
    cg_monte_carlo.registerObserver(CG_monte_carlo.EnergyObserver(cg_monte_carlo, file_name=energy_file))
    cg_monte_carlo.registerObserver(
        CG_monte_carlo.StructureObserver(cg_monte_carlo, file_name=structure_file, pdb_file=pdb_file))
    cg_monte_carlo.run()

    with open(energy_file) as f:
        lines = f.readlines()
    assert [int(line.split(",")[0]) for line in lines] == [10, 20, 30]
    assert os.path.isfile(pdb_file)


def test_observers_remove_output_of_previous_run(cg_monte_carlo, tmp_path):
    """
    Output left behind by an interrupted simulation must not be appended to
    when the simulation is restarted (issue #34)
    """
    stale = {}
    for name in ("energies.txt", "structure.txt", "trajectory.pdb", "kt_schedule.txt", "acc_ratio.txt"):
        stale[name] = str(tmp_path / name)
        with open(stale[name], "w") as f:
            f.write("output of an interrupted run\n")

    CG_monte_carlo.EnergyObserver(cg_monte_carlo, file_name=stale["energies.txt"])
    CG_monte_carlo.StructureObserver(cg_monte_carlo, file_name=stale["structure.txt"],
                                     pdb_file=stale["trajectory.pdb"])
    CG_monte_carlo.KTObserver(cg_monte_carlo, file_name=stale["kt_schedule.txt"])
    CG_monte_carlo.AcceptanceRatioObserver(cg_monte_carlo, file_name=stale["acc_ratio.txt"])

    assert not any(os.path.isfile(filename) for filename in stale.values())


def test_annealer_runs_schedule(pose, score_function, seq_mover):
    """
    Small version of the annealing workflow used in the simulations directory
    """
    annealer_params = CG_monte_carlo.CGMonteCarloAnnealerParameters(
        n_inner=20,
        t_init=10,
        anneal_rate=0.5,
        n_anneals=3,
        annealer_criteron=CG_monte_carlo.Repeat1Convergence,
        mc_output=False,
        out_freq=20,
    )
    assert annealer_params.t_anneals == [10, 5, 2.5]

    start_energy = score_function(pose)
    annealer = CG_monte_carlo.CGMonteCarloAnnealer(
        seq_mover=seq_mover,
        score_function=score_function,
        pose=pose,
        param_file_object=annealer_params,
    )
    energy_observer = CG_monte_carlo.EnergyObserver(annealer.get_mc_sim(), write_file=False)
    annealer.registerObserver(energy_observer)
    annealer.run_schedule()

    mc_sim = annealer.get_mc_sim()
    assert mc_sim.kT == 2.5
    assert mc_sim.mc.total_trials() == 60
    assert len(energy_observer.energies) == 3
    assert annealer.kt_anneals == []
    assert score_function(mc_sim.get_minimum_energy_pose()) <= start_energy + 1e-8
