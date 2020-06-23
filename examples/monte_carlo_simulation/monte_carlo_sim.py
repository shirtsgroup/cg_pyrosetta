# Simple Monte Carlo script using the FoldingAlgorithm object to
# fold an example CG model

import cg_pyrosetta
import pyrosetta
pyrosetta.init()

# Set desired parameters
cg_pyrosetta.change_parameters.changeTorsionParameters(
    {
     'CG1 CG1 CG1 CG1':[0,0,0],
     'CG2 CG1 CG1 CG1':[0,0,0],
     'CG2 CG1 CG1 CG1':[0,0,0]
     }
)

cg_pyrosetta.change_parameters.changeAngleParameters(
    {
     'CG1 CG1 CG1':[0,0,0],
     'CG2 CG1 CG1':[0,0,0]
    }
)

cg_pyrosetta.change_parameters.changeAtomParameters(
    {
        'CG1':['X', 1.0, 0.2],
        'CG2':['X', 0.5, 0.2]
    }
)

# Setup CGFoldingAlgorithm object

sequence = "X[CG11x3:CGLower]X[CG11x3]X[CG11x3]X[CG11x3:CGUpper]"
kt = 1

cg_folding_object = cg_pyrosetta.CG_folding.CGFoldingAlgorithm(sequence)

# Create CG Movers
cg_small = cg_pyrosetta.CG_movers.CGSmallMover(cg_folding_object.pose)
cg_small.angle = 180

mini = pyrosetta.rosetta.protocols.minimization_packing.MinMover()
movemap = pyrosetta.MoveMap()
mini.score_function(cg_folding_object.scorefxn)
movemap.set_bb_true_range(1, cg_folding_object.pose.size())
mini.movemap(movemap)

pymol = pyrosetta.PyMOLMover()


# Add new folding moves
cg_folding_object.build_fold_alg("TorsionMC")
cg_folding_object.add_folding_move("TorsionMC", pyrosetta.RepeatMover(cg_small, 10))
cg_folding_object.add_folding_move("TorsionMC", pyrosetta.RepeatMover(mini, 20))
cg_folding_object.add_folding_move("TorsionMC", pymol)

# Run MC Simulation
cg_folding_object.run_folding_alg("TorsionMC", 10000)

# Write minimum energy structure
cg_folding_object.mc.lowest_score_pose().dump_pdb("outputs/min_energy.pdb")