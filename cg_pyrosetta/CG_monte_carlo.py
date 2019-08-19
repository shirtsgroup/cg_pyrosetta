# This file will contain the object used for running simple CG MC simulations
# Will require :
# 1) ScoreFucntion object
# 2) Pose
# 3) Set of Movers (SequenceMover object)
# 4) kT
# 5) n_steps
# 6) output_freq?
from abc import ABC, abstractmethod
import CG_movers
import os
import sys
import warnings

current_path = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.abspath(current_path + '/../PyRosetta4.modified'))

import pyrosetta
import numpy as np


class CGMonteCarlo:
    """
    Docstring here
    """
    def __init__(self,
        pose : object = None, # pyrosetta.rosetta.core.pose.Pose, Pose of starting structure
        score : object = None, # pyrosetta.ScoreFunction, # scorefunction from pyrosetta
        seq_mover : object = None, #pyrosetta.rosetta.protocols.moves.SequenceMover, # sequence of moves between MC evaluations
        n_steps : int = 1000000,
        kT : float = 1,
        output : bool = True,):

        # initialize input values
        self.pose = pose
        self.score = score
        self.seq_mover = seq_mover
        self._kT = kT
        self.n_steps = n_steps
        self.output = output

        # Build MC Object
        self.mc = pyrosetta.MonteCarlo(self.pose, self.score, self._kT)
        self.mc_trial = pyrosetta.TrialMover(self.seq_mover, self.mc)

        if self.output:
            self.pymol = pyrosetta.PyMOLMover()
            print("Initial Energy :", self.get_energy())
    @property
    def kT(self):
        return(self._kT)
    
    @kT.setter
    def kT(self, kT):
        self._kT = kT
        self.mc.set_temperature(self.kT)
        self.mc_trial = pyrosetta.TrialMover(self.seq_mover, self.mc)
    
    def get_energy(self):
        return(self.score(self.pose))

    # def __call__():
    def run(self):
        run = pyrosetta.RepeatMover(self.mc_trial, 1)
        for i in range(self.n_steps):
            run.apply(self.pose)
            if self.output:

                print("Step :", i)
                print("Energy : ", self.get_energy())
                self.pymol.apply(self.pose)
    
    def get_pose(self):
        return(self.pose)

    def get_minimum_energy_pose(self):
        return(self.mc.lowest_score_pose())

    
class CGMonteCarloScheduler:
    """
    Docstring here
    """
    def __init__(self,
            seq_mover_maker : object = None,
            energy_builder : object = None,
            pose_adapter : object = None,
            param_file_object : object = None: 
            ):

        self.seq_builder = seq_mover_maker
        self.energy_builder = energy_builder
        self.pose = pose_adapter.get_pose()
        self.param_file_object = param_file_object 
        # ^^^^^^^^
        # These will hold the set of instructions the scheduler will have to
        # follow


        
    def _get_score_function():
        pass

    def _get_seq_mover():
        pass
    
    def _build_MC_job():
        pass

    def run_schedule():
        pass

    def _add_to_schedule():
        pass





# class PoseAdapter(ABC):
#     @abstractmethod
#     def get_pose(self):
#         pass

class SequenceMoverFactory:

    def __init__(self):
        self.methods = {
                        'small_dihe': CG_movers.CGSmallMover,
                        'small_angle': CG_movers.CGSmallAngleMover,
                        'shear_dihe': CG_movers.CGShearMover,
                        'sc_small_dihe': CG_movers.CGSmallSCMover,
                        'sc_small_angle': CG_movers.CGSmallAngleMover,
                        }

    def build_seq_mover(self, pose, mover_list, freq_list):
        seq_mover = pyrosetta.SequenceMover()
        for mover, freq in zip(mover_list, freq_list):
            if mover in self.methods.keys():
                seq_mover.add_mover(pyrosetta.RepeatMover(self.methods[mover](pose), freq))
            else:
                warnings.warn("Unimplemented Mover : "+mover+"\n Skipping mover", UserWarning)
    
        return(seq_mover)



class EnergyFunctionFactory:

    def __init__(self):
        methods = {}
        for method in dir(pyrosetta.rosetta.core.scoring):
            methods[method] = eval('pyrosetta.rosetta.core.scoring.'+method)
        
    def build_energy_function(self, score_terms, term_weights):
        score = pyrosetta.ScoreFunction()
        
        for term, weight in zip(score_terms, term_weights):
            if term in self.methods.keys():
                score.set_weight(eval('pyrosetta.rosetta.core.scoring.'+term), weight)
            else:
                warnings.warn("Energy Term not implemented :"+term+"\n Skipping term")

        return(score)




    