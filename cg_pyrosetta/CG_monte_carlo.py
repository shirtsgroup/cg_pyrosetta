# This file will contain the object used for running simple CG MC simulations
# Will require :
# 1) ScoreFucntion object
# 2) Pose
# 3) Set of Movers (SequenceMover object)
# 4) kT
# 5) n_steps
# 6) output_freq?

import os
import sys
import warnings

current_path = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.abspath(current_path + '/../PyRosetta4.modified'))

import pyrosetta
import numpy as np


class CGMonteCarlo:
    def __init__(self,
        pose : object = pyrosetta.rosetta.core.pose.Pose, # Pose of starting structure
        score : object = pyrosetta.ScoreFunction, # scorefunction from pyrosetta
        seq_mover : object = pyrosetta.rosetta.protocols.moves.SequenceMover, # sequence of moves between MC evaluations
        n_steps : int,
        kT : float):

        self.pose = pose
        self.score = score
        self.seq_mover = seq_mover
        self.kT = kT
        self.n_steps = n_steps

        pass

    # def __call__():
    def run(self):
        for i in range(self.n_steps):
            
        pass

    def get_kt(self):
        return(self.kT)

    def _set_kt(self, kT):
        self.kT = kT

    def get_energy(self):
        return(self.score(self.pose))


    