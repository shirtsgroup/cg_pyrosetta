# This example file will be used to show how to generally use the CG_folding.FoldingAlgorithm object
from multiprocessing import Pool
import cg_pyrosetta
import pyrosetta

def run_anneal_folding(sequence, name,rep, kt_anneal):
    # first generate a folded structure using CGFoldingAlgorithm
    folding_object = cg_pyrosetta.CG_folding.CGFoldingAlgorithm(sequence)

    # If running PyMOL this will ensure structure output during MC simulation
    # 'default' is the folding algorithm selected
    # this algorithn consists of:
    # 10x CGSmallMover
    # 10x CGShearMober
    # 10x MinMover
    # MC evaluation

    # Adding a pymol object to the folding algorithm so we can confirm
    # folding is happening
    folding_object.build_fold_alg('no_min')
    
    # Editing the PDB_writer object in our folding object so we can get an output of
    # the annealing process
    folding_object.add_folding_move('no_min', pyrosetta.RepeatMover(folding_object.small, 10))
    folding_object.add_folding_move('no_min', pyrosetta.RepeatMover(folding_object.shear, 10))
    folding_object.add_folding_move('no_min', folding_object.pymol)

    # Runs a folding MC simulation with 200 repeats of the 'default' folder at each kt
    folding_object.run_anneal_fold('no_min', 5000, kt_anneal)

    # Dump the lowest energy structure from the MC simulation
    folding_object.mc.lowest_score_pose().dump_pdb('outputs/'+name+'_example_'+str(rep)+'.pdb')

def multiprocess_anneal_folding(sequence_name_rep_kts):
    run_anneal_folding(sequence_name_rep_kts[0], sequence_name_rep_kts[1], sequence_name_rep_kts[2], sequence_name_rep_kts[3])
    return(sequence_name_rep_kts[2])


def build_inputs(seqences, name, rep, kts):
    """
    object builds the various replicas of a given folding run we want
    """
    input_list = []
    for sequence, name in zip(sequences, names):
        for i in range(rep):
            input_list.append([sequence, name, i, kts])
    return(input_list)


if __name__ == '__main__':
    p = Pool(5)

    # list of sequences of several CG models [CG11*5, CG21*5, CG31*5, mixed]
    sequences = ['X[CG11]X[CG11]X[CG11]X[CG11]X[CG11]X[CG11]X[CG11]X[CG11]X[CG11]X[CG11]X[CG11]X[CG11]X[CG11]X[CG11]X[CG11]',
            'X[CG21]X[CG21]X[CG21]X[CG21]X[CG21]X[CG21]X[CG21]X[CG21]X[CG21]X[CG21]X[CG21]X[CG21]X[CG21]X[CG21]X[CG21]',
            'X[CG31]X[CG31]X[CG31]X[CG31]X[CG31]X[CG31]X[CG31]X[CG31]X[CG31]X[CG31]X[CG31]X[CG31]X[CG31]X[CG31]X[CG31]',
            'X[CG11]X[CG21]X[CG31]X[CG21]X[CG11]X[CG11]X[CG21]X[CG31]X[CG21]X[CG11]X[CG11]X[CG21]X[CG31]X[CG21]X[CG11]']

    names = ['CG11', 'CG21', 'CG31', 'mixed']
    # Defining the various kt values used over course of sim.
    kt_i = 100
    kt_anneal = [kt_i*(0.9)**i for i in range(50)]
    inputs = build_inputs(sequences, names, 30, kt_anneal)
    for i in p.imap_unordered(multiprocess_anneal_folding, inputs):
        print(i)
