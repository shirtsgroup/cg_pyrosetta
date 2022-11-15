from multiprocessing import Pool
import cg_pyrosetta
import pyrosetta
import yaml
import time
import timeit

# This code is similar to fold_sequence_multiprocess.py, except we time the specific runs

def run_anneal_folding(sequence, name,rep, kt_anneal):

    import cg_pyrosetta
    import pyrosetta

    # first generate a folded structure using CGFoldingAlgorithm
    folding_object = cg_pyrosetta.CG_folding.CGFoldingAlgorithm(sequence)

    # If running PyMOL this will ensure structure output during MC simulation
    # 'default' is the folding algorithm selected
    # this algorithn consists of:
    # 10x CGSmallMover
    # 10x CGShearMober
    # 10x MinMover
    # MC evaluation
    folding_object.build_fold_alg('no_min')
    folding_object.add_folding_move('no_min', pyrosetta.RepeatMover(folding_object.shear, 10))
    folding_object.add_folding_move('no_min', pyrosetta.RepeatMover(folding_object.small_angle, 5))
    folding_object.add_folding_move('no_min', folding_object.pymol)

    # Runs a folding MC simulation with 200 repeats of the 'default' folder at each kt

    folding_object.run_anneal_fold('no_min', 1000, kt_anneal)
    
    # Dump the lowest energy structure from the MC simulation
    return(folding_object.mc.total_trials())

def multiprocess_anneal_folding(sequence_name_rep_kts):
    time = run_anneal_folding(sequence_name_rep_kts[0], sequence_name_rep_kts[1], sequence_name_rep_kts[2], sequence_name_rep_kts[3])
    return(sequence_name_rep_kts[1], time)


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

    cg_pyrosetta.parameters.changeTorsionParameters(
        {'CG1 CG1 CG1 CG1':[5,3,0],
            'CG2 CG1 CG1 CG2':[0,0,0],
            'CG2 CG1 CG1 CG1':[0,0,0],
            'X CG2 CG1 CG1':[0,0,0]},
    )

    cg_pyrosetta.parameters.changeAngleParameters(
        {'CG1 CG1 CG1':[20,120],
        'CG2 CG1 CG1':[20,120],
        'CG1 CG1 CG2':[0,0],
        'X CG2 CG1':[0,0]}        
    )

    cg_pyrosetta.builder.buildCGPyRosetta()

    # list of sequences of several CG models [CG11*5, CG21*5, CG11*5, mixed]
    sequences = [
            'X[CG11]X[CG11]X[CG11]X[CG11]X[CG11]',
            'X[CG11]X[CG11]X[CG11]X[CG11]X[CG11]X[CG11]X[CG11]X[CG11]X[CG11]X[CG11]',
            'X[CG11]X[CG11]X[CG11]X[CG11]X[CG11]X[CG11]X[CG11]X[CG11]X[CG11]X[CG11]X[CG11]X[CG11]X[CG11]X[CG11]X[CG11]',
            'X[CG11]X[CG11]X[CG11]X[CG11]X[CG11]X[CG11]X[CG11]X[CG11]X[CG11]X[CG11]X[CG11]X[CG11]X[CG11]X[CG11]X[CG11]X[CG11]X[CG11]X[CG11]X[CG11]X[CG11]',
            'X[CG21]X[CG21]X[CG21]X[CG21]X[CG21]',
            'X[CG21]X[CG21]X[CG21]X[CG21]X[CG21]X[CG21]X[CG21]X[CG21]X[CG21]X[CG21]',
            'X[CG21]X[CG21]X[CG21]X[CG21]X[CG21]X[CG21]X[CG21]X[CG21]X[CG21]X[CG21]X[CG21]X[CG21]X[CG21]X[CG21]X[CG21]',
            'X[CG21]X[CG21]X[CG21]X[CG21]X[CG21]X[CG21]X[CG21]X[CG21]X[CG21]X[CG21]X[CG21]X[CG21]X[CG21]X[CG21]X[CG21]X[CG21]X[CG21]X[CG21]X[CG21]X[CG21]',
            'X[CG31]X[CG31]X[CG31]X[CG31]X[CG31]',
            'X[CG31]X[CG31]X[CG31]X[CG31]X[CG31]X[CG31]X[CG31]X[CG31]X[CG31]X[CG31]',
            'X[CG31]X[CG31]X[CG31]X[CG31]X[CG31]X[CG31]X[CG31]X[CG31]X[CG31]X[CG31]X[CG31]X[CG31]X[CG31]X[CG31]X[CG31]',
            'X[CG31]X[CG31]X[CG31]X[CG31]X[CG31]X[CG31]X[CG31]X[CG31]X[CG31]X[CG31]X[CG31]X[CG31]X[CG31]X[CG31]X[CG31]X[CG31]X[CG31]X[CG31]X[CG31]X[CG31]',
            ]

    names = ['CG11_5mer', 'CG11_10mer', 'CG11_15mer', 'CG11_20mer',
             'CG21_5mer', 'CG21_10mer', 'CG21_15mer', 'CG21_20mer',
             'CG31_5mer', 'CG31_10mer', 'CG31_15mer', 'CG31_20mer']
    # Defining the various kt values used over course of sim.
    kt_i = 100
    kt_anneal = [kt_i*(0.9)**i for i in range(100)]
    inputs = build_inputs(sequences, names, 1, kt_anneal)
    average_times = {}
    for out in p.imap_unordered(multiprocess_anneal_folding, inputs):
        average_times[out[0]] = out[1]
        with open('outputs/mc_steps.yml','w') as outfile:
            yaml.dump(average_times, outfile)

    with open('outputs/mc_steps.yml','w') as outfile:
        yaml.dump(average_times, outfile)

