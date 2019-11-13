# This example file will be used to show how to read in a cif and sample the local backbone space around the starting conformation
import pytest
import cg_pyrosetta
import pyrosetta
import os
import numpy as np

# first generate a folded structure using CGFoldingAlgorithm

def test_read_in_pdb_functional():
    cg_pyrosetta.change_parameters.changeAtomParameters({"CG1":['X', 1.0, 0.2],
                                                         "CG2":['X', 1.0, 0.2]})
    cg_pyrosetta.change_parameters.changeAngleParameters({"CG1 CG1 CG1":[0, 0],
                                                          "CG2 CG1 CG1":[0, 0]})
    cg_pyrosetta.change_parameters.changeTorsionParameters({"CG1 CG1 CG1 CG1":[0, 0, 0],
                                                            "CG2 CG1 CG1 CG2":[0, 0, 0]})
    cg_pyrosetta.builder.buildCGPyRosetta()

    folding_object = cg_pyrosetta.CG_folding.CGFoldingAlgorithm('X[CG11x3:CGLower]X[CG11x3]X[CG11x3]X[CG11x3]X[CG11x3:CGUpper]')
    # Defining the various kt values used over course of sim.
    kt_i = 100
    kt_anneal = [kt_i*(0.9)**i for i in range(5)]

    # If running PyMOL this will ensure structure output during MC simulation
    folding_object.add_folding_move('default', folding_object.pymol)

    # Runs a folding MC simulation with 200 repeats of the 'default' folder at each kt
    folding_object.run_anneal_fold('default', 1, kt_anneal)

    # Dump the lowest energy structure from the MC simulation
    pose_from_seq = folding_object.mc.lowest_score_pose()
    pose_from_seq.dump_cif('1_1_example.cif')
    pose_from_seq.dump_pdb('1_1_example.cif')
    pose_from_cif = pyrosetta.pose_from_file('1_1_example.cif')
    folding_object.scorefxn(pose_from_cif)
    pose_from_cif.dump_cif('1_1_example_readin.cif')
    print(folding_object.scorefxn(pose_from_cif), folding_object.scorefxn(pose_from_seq))
    np.testing.assert_almost_equal(folding_object.scorefxn(pose_from_cif), folding_object.scorefxn(pose_from_seq), decimal=1)

    # Teardown
    os.remove("1_1_example.cif")
    os.remove("1_1_example_readin.cif")
