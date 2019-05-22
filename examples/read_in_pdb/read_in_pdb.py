# This example file will be used to show how to read in a PDB and sample the local backbone space around the starting conformation

import cg_pyrosetta
from cg_pyrosetta.CG_folding import pyrosetta


# first generate a folded structure using CGFoldingAlgorithm
folding_object = cg_pyrosetta.CG_folding.CGFoldingAlgorithm('X[CG11]X[CG11]X[CG11]X[CG11]X[CG11]')

# Defining the various kt values used over course of sim.
kt_i = 100
kt_anneal = [kt_i*(0.9)**i for i in range(50)]

# If running PyMOL this will ensure structure output during MC simulation
folding_object.add_folding_move('default', folding_object.pymol)

# Runs a folding MC simulation with 200 repeats of the 'default' folder at each kt
folding_object.run_anneal_fold('default', 200, kt_anneal)

# Dump the lowest energy structure from the MC simulation
folding_object.mc.lowest_score_pose().dump_pdb('1_1_example.pdb')

pose_from_pdb = pyrosetta.pose_from_pdb('1_1_example.pdb')

assert folding_object.scorefxn(pose_from_pdb) == folding_object.scorefxn(folding_object.pose)