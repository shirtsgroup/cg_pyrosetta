# This example file will be used to show how to generally use the CG_folding.FoldingAlgorithm object
import cg_pyrosetta
import pyrosetta
pyrosetta.init()


# list of sequences of several CG models [CG11*5, CG21*5, CG31*5, mixed]
sequences = [
        #  'X[CG13]X[CG13]X[CG13]X[CG13]X[CG13]X[CG13]X[CG13]X[CG13]X[CG13]X[CG13]X[CG13]X[CG13]X[CG13]X[CG13]X[CG13]',
        #  'X[CG21]X[CG21]X[CG21]X[CG21]X[CG21]X[CG21]X[CG21]X[CG21]X[CG21]X[CG21]X[CG21]X[CG21]X[CG21]X[CG21]X[CG21]',
           'X[CG31:CGLower]X[CG11x3]X[CG31]X[CG11x3]X[CG31]X[CG11x3]X[CG31]X[CG11x3]X[CG31]X[CG11x3]X[CG31]X[CG11x3]X[CG31]X[CG11x3]X[CG31]X[CG11x3]X[CG31]X[CG11x3]X[CG31]X[CG11x3]X[CG31]X[CG11x3]X[CG31]X[CG11x3]X[CG31]X[CG11x3]X[CG31]X[CG11x3:CGUpper]',
        #  'X[CG12]X[CG12]X[CG12]X[CG12]X[CG12]X[CG12]X[CG12]X[CG12]X[CG12]X[CG12]X[CG12]X[CG12]X[CG12]X[CG12]X[CG12]',
        #  'X[CG11]X[CG11]X[CG11]X[CG11]X[CG11]X[CG11]X[CG11]X[CG11]X[CG11]X[CG11]X[CG11]X[CG11]X[CG11]X[CG11]X[CG11]',
        #  'X[CG11]X[CG12]X[CG21]X[CG31]X[CG13]X[CG21]X[CG13]X[CG11]X[CG11]X[CG21]X[CG13]X[CG11]X[CG12]X[CG21]X[CG31]'
]

names = [# 'CG13',
         # 'CG21',
         'CG31',
         # 'CG12',
         # 'CG11', 
         # 'mixed'
        ]

# Defining the various kt values used over course of sim.
kt_i = 100
kt_anneal = [kt_i*(0.9)**i for i in range(50)]


cg_pyrosetta.change_parameters.changeTorsionParameters(
   {'CG1 CG1 CG1 CG1':[3,3,0],
    'CG2 CG1 CG1 CG2':[3,2,120],
    'CG2 CG1 CG1 CG1':[0,0,0],
    'X CG2 CG1 CG1':[0,0,0]},
    )

cg_pyrosetta.change_parameters.changeAngleParameters(
    {'CG1 CG1 CG1':[2,70],
     'CG2 CG1 CG1':[2,70],
     'CG1 CG1 CG2':[0,0],
      'X CG2 CG1':[0,0]}        
)

cg_pyrosetta.builder.buildCGPyRosetta()

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

        # folding_object.add_folding_move('default', folding_object.pymol)
        
        folding_object.build_fold_alg('AngleMC')

        folding_object.add_folding_move('AngleMC', pyrosetta.RepeatMover(folding_object.small, 1))
        # folding_object.add_folding_move('AngleMC', pyrosetta.RepeatMover(folding_object.shear, 5))

        # Adding an angle mover to this folding algorithm
        small_sc_mover = cg_pyrosetta.CG_movers.CGSmallSCMover(folding_object.pose)
        small_sc_mover.angle = 180

        # repeat_sc_mover = pyrosetta.RepeatMover(small_sc_mover, 5)
        # folding_object.add_folding_move('AngleMC', repeat_sc_mover)

        small_bb_angle_mover = cg_pyrosetta.CG_movers.CGSmallAngleMover(folding_object.pose)
        small_bb_angle_mover.angle = 10
        # repeat_bb_angle_mover = pyrosetta.RepeatMover(small_bb_angle_mover, 2)
        # folding_object.add_folding_move('AngleMC', repeat_bb_angle_mover)

        # Adding a PyMOL object to the folding sequence so we can see output
        
        folding_object.add_folding_move('AngleMC', pyrosetta.RepeatMover(folding_object.mini, 100))
                
        pymol = pyrosetta.PyMOLMover()
        folding_object.add_folding_move('AngleMC', pymol)


        # Runs a folding MC simulation with 200 repeats of the 'default' folder at each kt
        folding_object.run_anneal_fold('AngleMC', 5000, kt_anneal)

        # Dump the lowest energy structure from the MC simulation
        folding_object.mc.lowest_score_pose().dump_pdb('outputs/'+names[i]+'_example_angles_'+str(rep)+'.pdb')
        
        twist_score = pyrosetta.ScoreFunction()
        twist_score.set_weight(pyrosetta.rosetta.core.scoring.mm_twist, 1)
        # print(twist_score(folding_object.pose))

