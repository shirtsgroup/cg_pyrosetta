# This example file will be used to show how to generally use the CG_folding.FoldingAlgorithm object
import cg_pyrosetta
from cg_pyrosetta.CG_folding import pyrosetta as pyrosetta

# list of sequences of several CG models [CG11*5, CG21*5, CG31*5, mixed]
sequences = ['X[CG11]X[CG11]X[CG11]X[CG11]X[CG11]X[CG11]X[CG11]X[CG11]X[CG11]X[CG11]X[CG11]X[CG11]X[CG11]X[CG11]X[CG11]',
          'X[CG21]X[CG21]X[CG21]X[CG21]X[CG21]X[CG21]X[CG21]X[CG21]X[CG21]X[CG21]X[CG21]X[CG21]X[CG21]X[CG21]X[CG21]',
          'X[CG31]X[CG31]X[CG31]X[CG31]X[CG31]X[CG31]X[CG31]X[CG31]X[CG31]X[CG31]X[CG31]X[CG31]X[CG31]X[CG31]X[CG31]',
          'X[CG11]X[CG21]X[CG31]X[CG21]X[CG11]X[CG11]X[CG21]X[CG31]X[CG21]X[CG11]X[CG11]X[CG21]X[CG31]X[CG21]X[CG11]',
          'X[CG12]X[CG12]X[CG12]X[CG12]X[CG12]X[CG12]X[CG12]X[CG12]X[CG12]X[CG12]X[CG12]X[CG12]X[CG12]X[CG12]X[CG12]',
          'X[CG13]X[CG13]X[CG13]X[CG13]X[CG13]X[CG13]X[CG13]X[CG13]X[CG13]X[CG13]X[CG13]X[CG13]X[CG13]X[CG13]X[CG13]',
          'X[CG11]X[CG12]X[CG21]X[CG31]X[CG13]X[CG21]X[CG13]X[CG11]X[CG11]X[CG21]X[CG13]X[CG11]X[CG12]X[CG21]X[CG31]']

names = ['CG11', 'CG21', 'CG31', 'CG12', 'CG13', 'mixed']

# Defining the various kt values used over course of sim.
kt_i = 100
kt_anneal = [kt_i*(0.9)**i for i in range(50)]

for rep in range(5):

    for i in range(len(sequences)):
        # first generate a folded structure using CGFoldingAlgorithm
        folding_object = cg_pyrosetta.CG_folding.CGFoldingAlgorithm(sequences[i])

        # If running PyMOL this will ensure structure output during MC simulation
        # 'default' is the folding algorithm selected
        # this algorithn consists of:
        # 10x CGSmallMover
        # 10x CGShearMober
        # 10x MinMover
        # MC evaluation
        folding_object.scorefxn.set_weight(pyrosetta.rosetta.core.scoring.mm_twist, 0.1)
        folding_object.add_folding_move('default', folding_object.pymol)
        
        # Adding an angle mover to this folding algorithm
        small_sc_mover = cg_pyrosetta.CG_movers.CGSmallSCMover(folding_object.pose)
        small_sc_mover.angle = 180
        repeat_sc_mover = pyrosetta.RepeatMover(small_sc_mover, 10)
        folding_object.add_folding_move('default', repeat_sc_mover)

        # Runs a folding MC simulation with 200 repeats of the 'default' folder at each kt
        folding_object.run_anneal_fold('default', 5000, kt_anneal)

        # Dump the lowest energy structure from the MC simulation
        folding_object.mc.lowest_score_pose().dump_pdb('outputs/'+names[i]+'_example_angles_'+str(rep)+'.pdb')
