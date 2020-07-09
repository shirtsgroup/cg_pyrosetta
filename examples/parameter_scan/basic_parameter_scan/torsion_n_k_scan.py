#! /home/tfobe/anaconda3/envs/python36/bin/python

# This code will be used to perform Rosetta Parameter
# Searches on simple molecules

import numpy as np
from multiprocessing import Pool
import cg_pyrosetta
import pyrosetta
import mdtraj as md
import os


def initializeParameters():
    cg_pyrosetta.change_parameters.changeTorsionParameters(
    {'CG1 CG1 CG1 CG1':[0,0,0],
    'CG2 CG1 CG1 CG2':[0,0,0],
    'CG2 CG1 CG1 CG1':[0,0,0],
    'X CG2 CG1 CG1':[0,0,0]}
    )
    cg_pyrosetta.change_parameters.changeAngleParameters(
    {'CG1 CG1 CG1':[0,0],
     'CG2 CG1 CG1':[0,0],
     'X CG2 CG1':[0,0]},
    )


def buildTorsionParamDict(name, force_constant, periodicity):
    param_dict = {name:[force_constant, periodicity, 0]}
    return(param_dict)


def updateParameters(param_dict):
    cg_pyrosetta.change_parameters.changeTorsionParameters(param_dict)
    cg_pyrosetta.builder.buildCGPyRosetta()
    

def runAnnealingProcess(param_dict, rep, kts, k, p):
    folding_object = cg_pyrosetta.CG_folding.CGFoldingAlgorithm('X[CG11x3]X[CG11x3]X[CG11x3]X[CG11x3]X[CG11x3]')
    # folding_object.scorefxn.set_weight(pyrosetta.rosetta.core.scoring.mm_twist, 1)
    # folding_object.scorefxn.set_weight(pyrosetta.rosetta.core.scoring.mm_bend, 1)
    folding_object.add_folding_move('default', folding_object.pymol)
    folding_object.run_anneal_fold('default', 100, kts)
    folding_object.mc.lowest_score_pose().dump_pdb(os.path.join('outputs', 'force_constant_'+str(round(k, 3)), 'period_'+str(p), 'CG31_rep_'+str(rep)+'.pdb'))



def multiprocess_function(list_of_params):
    runAnnealingProcess(list_of_params[0], list_of_params[1], list_of_params[2], list_of_params[3], list_of_params[4])
    return(list_of_params[1])
    
def combineReplicas(k, p, reps):
    file_names = [os.path.join('outputs', 'force_constant_'+str(round(k, 3)), 'period_'+str(p), 'CG31_rep_'+str(rep)+'.pdb') for rep in range(1, reps+1)]
    traj = md.load(file_names)
    traj.save_dcd(os.path.join('outputs', 'force_constant_'+str(round(k, 3)), 'period_'+str(p), 'CG31.dcd'))
    for rm_file in file_names[1:]:
        os.remove(rm_file)


def main():
    # Reset all parameters to zero

    # Would like to generalize this work flow to:

    # 1) select parameters you would like to scan (modularly) (e.g. [torsion_force_k, BB_LJ_epsilon, SC_LJ_sigma])
    # 2) select ranges for parameter scan (e.g. np.linspace(), np.linspace(), np.linspace(), np.linspace())
    # 3) automatically builds loops/simulations for this
    reps = 30
    initializeParameters()
    torsion_name = 'CG2 CG1 CG1 CG2'
    k_torsions = np.linspace(0, 15, 20)
    periodicities  = np.arange(1,6)

    kts = [0.59616*(1/0.9)**i for i in range(50)]   # ends at 300K starting at a kt ~ 100 --> ~ 50K K
    kts.reverse()

    for k in k_torsions:
        for p in periodicities:
                if not os.path.exists(os.path.join('outputs', 'force_constant_'+str(round(k, 3)), 'period_'+str(p))):
                    os.makedirs(os.path.join('outputs', 'force_constant_'+str(round(k, 3)), 'period_'+str(p)))
                print('Force Constant : '+str(k))
                print('Periodicity : '+str(p))
                torsion_params = buildTorsionParamDict(torsion_name, k, p)
                updateParameters(torsion_params)
                multiprocess_params = [[torsion_params, rep, kts, k, p] for rep in range(1,reps+1)]
                pool = Pool(5)
                for i in pool.imap_unordered(multiprocess_function, multiprocess_params):
                    print('Working on rep:', i)

                

if __name__ == '__main__':
    main()