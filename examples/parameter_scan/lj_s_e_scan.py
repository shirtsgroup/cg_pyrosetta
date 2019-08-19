#! /home/tfobe/anaconda3/envs/python36/bin/python

# This code will be used to perform Rosetta Parameter
# Searches on simple molecules

import numpy as np
from multiprocessing import Pool
import mdtraj as md
from cg_pyrosetta.change_parameters import changeAtomParameters
import os
from importlib import reload

def initializeParameters():
    import cg_pyrosetta
    import pyrosetta

    cg_pyrosetta.change_parameters.changeTorsionParameters(
    {'CG1 CG1 CG1 CG1':[0,0,0],
    'CG2 CG1 CG1 CG2':[0,0,0],
    'CG2 CG1 CG1 CG1':[0,0,0],
    'X CG2 CG1 CG1':[0,0,0]}
    )
    cg_pyrosetta.change_parameters.changeAngleParameters(
    {'CG1 CG1 CG1':[0,0],
     'CG2 CG1 CG1':[0,0],
     'X CG2 CG1':[0,0]}
    )

    cg_pyrosetta.change_parameters.changeAtomParameters(
        {'CG1':['X', 1.0, 0.2],
         'CG2':['X', 1.0, 0.2],
         'CG3':['X', 1.0, 0.2],
         }
    )

def runAnnealingProcess(rep, kts, sigma, epsilon):
    import cg_pyrosetta
    import pyrosetta

    pyrosetta.init()
    
    folding_object = cg_pyrosetta.CG_folding.CGFoldingAlgorithm('X[CG11]X[CG11]X[CG11]X[CG11]X[CG11]X[CG11]X[CG11]X[CG11]X[CG11]X[CG11]X[CG11]X[CG11]X[CG11]X[CG11]X[CG11]')
    # folding_object.scorefxn.set_weight(pyrosetta.rosetta.core.scoring.mm_twist, 1)
    # folding_object.scorefxn.set_weight(pyrosetta.rosetta.core.scoring.mm_bend, 1)
    folding_object.add_folding_move('default', folding_object.pymol)
    folding_object.run_anneal_fold('default', 100, kts)
    folding_object.mc.lowest_score_pose().dump_pdb(os.path.join('outputs', 'sigma_'+str(round(sigma, 3)), 'epsilon_'+str(round(epsilon,3)), 'CG11_rep_'+str(rep)+'.pdb'))

def runMultiProcess(list_of_params):
    runAnnealingProcess(list_of_params[0], list_of_params[1], list_of_params[2], list_of_params[3])

def main():
    # Reset all parameters to zero

    # Would like to generalize this work flow to:

    # 1) select parameters you would like to scan (modularly) (e.g. [torsion_force_k, BB_LJ_epsilon, SC_LJ_sigma])
    # 2) select ranges for parameter scan (e.g. np.linspace(), np.linspace(), np.linspace(), np.linspace())
    # 3) automatically builds loops/simulations for this
    reps = 1
    initializeParameters()
    atom_name = 'CG2'
    sigmas = np.linspace(1, 5, 3)
    epsilons  = np.linspace(0.2, 10, 1)

    kts = [0.59616*(1/0.9)**i for i in range(10)]   # ends at 300K starting at a kt ~ 100 --> ~ 50K K
    kts.reverse()

    for sigma in sigmas:
        for epsilon in epsilons:
            if not os.path.exists(os.path.join('outputs', 'sigma_'+str(round(sigma, 3)), 'epsilon_'+str(round(epsilon,3)))):
                os.makedirs(os.path.join('outputs', 'sigma_'+str(round(sigma, 3)), 'epsilon_'+str(round(epsilon,3))))
            changeAtomParameters({atom_name:['X', sigma, epsilon]})
            multiprocess_params = [[rep, kts, sigma, epsilon] for rep in range(1,reps+1)]
            pool = Pool(1)

            for i in pool.imap_unordered(runMultiProcess, multiprocess_params):
                print("Working on rep:", i)
            
            file_names = [os.path.join('outputs', 'sigma_'+str(round(sigma, 3)), 'epsilon_'+str(round(epsilon, 3)), 'CG11_rep_'+str(rep)+'.pdb') for rep in range(1, reps+1)]
            traj = md.load(file_names)
            traj.save_dcd(os.path.join('outputs', 'sigma_'+str(round(sigma, 3)), 'epsilon_'+str(round(epsilon, 3)), 'CG11.dcd'))
            for rm_file in file_names[1:]:
                os.remove(rm_file)

if __name__ == "__main__":
    main()
