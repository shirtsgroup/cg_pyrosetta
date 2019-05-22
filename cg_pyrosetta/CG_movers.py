import os
import sys

current_path = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.abspath(current_path + '/../PyRosetta4.modified'))

import pyrosetta
import numpy as np

class CG11SmallMover(pyrosetta.rosetta.protocols.moves.Mover):
    """
    Analogous CG mover to the orginal "small" mover
    Will make a random dihedral move
    """
    def __init__(self, angle = 10):
        """
        Build Small Mover for CG polymers
        """
        self.angle = angle
        pyrosetta.rosetta.protocols.moves.Mover.__init__(self)
    
    def __str__(self):
        # pretty sure we don't need this
        pass
    
    def apply(self, pose):
        conf = pose.conformation()
        res = np.random.randint(1, pose.size()-3+1) # has to be able to move only bb atoms
        # print(res, res+1, res+2, res+3)
        d_angle = (np.random.rand()-0.5)*self.angle/180*np.pi
        # print(d_angle)
        old = conf.torsion_angle(pyrosetta.AtomID(1,res), pyrosetta.AtomID(1,res+1),  pyrosetta.AtomID(1,res+2),  pyrosetta.AtomID(1,res+3))
        new = old + d_angle
        conf.set_torsion_angle(pyrosetta.AtomID(1,res), pyrosetta.AtomID(1,res+1),  pyrosetta.AtomID(1,res+2),  pyrosetta.AtomID(1,res+3), new)


class CG11ShearMover(pyrosetta.rosetta.protocols.moves.Mover):
    """
    Analogous CG mover to the orginal "small" mover
    Will make a random dihedral move
    """
    def __init__(self, angle = 10):
        """
        Build Small Mover for CG polymers
        """
        self.angle = angle
        pyrosetta.rosetta.protocols.moves.Mover.__init__(self)
    
    def __str__(self):
        # pretty sure we don't need this
        pass
    
    def apply(self, pose):
        conf = pose.conformation()
        res = np.random.randint(1, pose.size()-4+1) # has to be able to move only bb atoms
        # print(res, res+1, res+2, res+3)
        d_angle = (np.random.rand()-0.5)*self.angle/180*np.pi
        # print(d_angle)
        old_1 = conf.torsion_angle(pyrosetta.AtomID(1,res), pyrosetta.AtomID(1,res+1),  pyrosetta.AtomID(1,res+2),  pyrosetta.AtomID(1,res+3))
        old_2 = conf.torsion_angle(pyrosetta.AtomID(1,res+1), pyrosetta.AtomID(1,res+2),  pyrosetta.AtomID(1,res+3),  pyrosetta.AtomID(1,res+4))
        new_1 = old_1 + d_angle
        new_2 = old_2 - d_angle
        conf.set_torsion_angle(pyrosetta.AtomID(1,res), pyrosetta.AtomID(1,res+1),  pyrosetta.AtomID(1,res+2),  pyrosetta.AtomID(1,res+3), new_1)
        conf.set_torsion_angle(pyrosetta.AtomID(1,res+1), pyrosetta.AtomID(1,res+2),  pyrosetta.AtomID(1,res+3),  pyrosetta.AtomID(1,res+4), new_2)


class CG11SmallAngleMover(pyrosetta.rosetta.protocols.moves.Mover):
    """
    Analogous CG mover to the orginal "small" mover
    Will make a random angle move
    """
    def __init__(self, angle = 10):
        """
        Build Small Mover for CG polymers
        """
        self.angle = angle
        pyrosetta.rosetta.protocols.moves.Mover.__init__(self)
    
    def __str__(self):
        # pretty sure we don't need this
        pass
    
    def apply(self, pose):
        conf = pose.conformation()
        res = np.random.randint(1, pose.size()-1+1) # has to be able to move only bb atoms
        # print(res, res+1, res+2)
        d_angle = (np.random.rand()-0.5)*self.angle/180*np.pi
        # print(d_angle)
        old = conf.bond_angle(pyrosetta.AtomID(1,res), pyrosetta.AtomID(1,res+1),  pyrosetta.AtomID(2,res+1))
        new = old + d_angle
        conf.set_bond_angle(pyrosetta.AtomID(1,res), pyrosetta.AtomID(1,res+1),  pyrosetta.AtomID(2,res+1), new)


class CG21SmallMover(pyrosetta.rosetta.protocols.moves.Mover):
    """
    Analogous CG mover to the orginal "small" mover for a 2-1 (bb-sc) model
    Will make a random dihedral move
    """
    def __init__(self, angle = 10):
        """
        Build Small Mover for CG polymers
        """
        self.angle = angle
        pyrosetta.rosetta.protocols.moves.Mover.__init__(self)
    
    def __str__(self):
        # pretty sure we don't need this
        pass
    
    def apply(self, pose):
        conf = pose.conformation()
        res = np.random.randint(1, pose.size()-3+1) # has to be able to move only bb atoms
        # print(res, res+1, res+2, res+3)
        d_angle = (np.random.rand()-0.5)*self.angle/180*np.pi
        # print(d_angle)
        old = conf.torsion_angle(pyrosetta.AtomID(1,res), pyrosetta.AtomID(1,res+1),  pyrosetta.AtomID(1,res+2),  pyrosetta.AtomID(1,res+3))
        new = old + d_angle
        conf.set_torsion_angle(pyrosetta.AtomID(1,res), pyrosetta.AtomID(1,res+1),  pyrosetta.AtomID(1,res+2),  pyrosetta.AtomID(1,res+3), new)


class CGSmallMover(pyrosetta.rosetta.protocols.moves.Mover):
    """
    Generalized CG mover analogous to the "small" mover in PyRosetta.
    """
    def __init__(self, pose, angle = 180):
        """
        Build Small Mover for CG polymers

        Arguments
        ---------

        angle : float
            maximum angle the mover can change an angle
        bb_model: int
            describes how many backbone beads are in a given CG model
        pose: pyrosetta.Pose()
            used to generate list of possible atoms
        """
        self.angle = angle
        pyrosetta.rosetta.protocols.moves.Mover.__init__(self)
        # Generate a list of all possible dihe in a provided pose with the bb_model given
        #  1 bb model [[1 1] [2 1] [3 1] [4 1]][]
        #  2 bb model [[1 1] [1 2] [2 1] [2 2]][ [3 1] [3 2]]
        #  3 bb model [[1 1] [1 2] [1 3] [2 1] [2 2] [2 3] [3 1] [3 2] [3 3]]
        conf = pose.conformation()

        # build list of all torsion angles in backbone
        self.atoms = []
        for i in range(1, pose.size()+1): # i == residue
            res_name = pose.residue(i).name()
            bb = int(res_name[2])  # hard coded... since we shouldn't need more than CG99 model 
            for j in range(1, bb+1): # j == atom in backbone

                if j == 1 and i == 1: # first residue requires different start atom 
                    self.atoms.append(pyrosetta.AtomID(bb+1, i))
                self.atoms.append(pyrosetta.AtomID(j, i))

                if i == pose.size() and bb == 1: # 1-1 specific modification for last atom
                    self.atoms.append(pyrosetta.AtomID(bb+1, i))  
        self.dihes = []
        for i in range(len(self.atoms)-3):  # builds dihedrals from atom list
            self.dihes.append([self.atoms[i], self.atoms[i+1], self.atoms[i+2], self.atoms[i+3]])

        # builds list of all atoms in backbone
        # self.atoms = []
        # for i in range(1, pose.size()+1):  # should be genearlizeable to any combination of residues
        #     res_name = pose.residue(i).name()
        #    bb = int(res_name[2])  # hard coded... since we shouldn't need more than CG99 model            
        #    for j in range(1, bb+1):
        #        self.atoms.append(pyrosetta.AtomID(j, i))

        # initializing conformation 
        for i in range(0, len(self.dihes)):
            old = conf.torsion_angle(self.dihes[i][0], self.dihes[i][1], self.dihes[i][2], self.dihes[i][3])
            pose.conformation().set_torsion_angle(self.dihes[i][0], self.dihes[i][1], self.dihes[i][2], self.dihes[i][3], old)

    def __str__(self):
        # pretty sure we don't need this
        pass

    def apply(self, pose):
        d_angle = (np.random.rand()-0.5)*self.angle/180*np.pi

        conf = pose.conformation()
        dihe_start = np.random.randint(0, len(self.atoms)-3) # has to be able to move only bb atoms
        old = conf.torsion_angle(self.atoms[dihe_start], self.atoms[dihe_start+1], self.atoms[dihe_start+2], self.atoms[dihe_start+3])
        new = old + d_angle
        # print('Changing Torsion',dihe_start, 'from', old, 'to', new)
        conf.set_torsion_angle(self.atoms[dihe_start], self.atoms[dihe_start+1], self.atoms[dihe_start+2], self.atoms[dihe_start+3], new)


class CGSmallSCMover(pyrosetta.rosetta.protocols.moves.Mover):
    """
    Generalized CG mover analogous to the "small" mover in PyRosetta.
    """
    def __init__(self, sc_model, pose, angle = 30):
        """
        Build Small Mover for CG polymers

        Arguments
        ---------

        angle : float, maximum angle the mover can change an angle
        bb_model: int, describes how many backbone beads are in a given CG model
        pose: pyrosetta.Pose(), used to generate list of possible atoms
        """
        self.angle = angle
        pyrosetta.rosetta.protocols.moves.Mover.__init__(self)
        # Generate a list of all possible dihe in a provided pose with the bb_model given
        #
        #  Legend : atoms --> [atomno, resno] residues --> {all possible atoms} 
        #
        #  1 sc model assert no backbone dihedral
        #
        #  2 sc model [{[1 2][1 1][2 1][3 1]} {[1 1][1 2][2 2][3 2]} {[1 2][1 3][3 2][3 3]}]
        #                First entry is a Unique case
        #  3 bb model [[1 1] [1 2] [1 3] [2 1] [2 2] [2 3] [3 1] [3 2] [3 3]]
        #

        conf = pose.conformation()
        self.atoms = []
        
        for i in range(1, pose.size()+1):  # should be genearlizeable to any combination of residues
            res_dihes = []
            res_name = pose.residue(i).name()
            sc = int(res_name[3])  # hard coded... since we shouldn't need more than CG99 model            
            if i == 1:
                # edge case for first residue
                res_dihes.append(pyrosetta.AtomID(1, 2))
            else:
                # else use first atom of previous residue
                res_dihes.append(pyrosetta.AtomID(1, i-1))

            for j in range(1, sc+1):                    
                res_dihes.append(pyrosetta.AtomID(j, i))
                # print(res, res+1, res+2, res+3)
        
        for i in range(0, len(self.atoms)-3):
            old = conf.torsion_angle(self.atoms[i], self.atoms[i+1], self.atoms[i+2], self.atoms[i+3])
            pose.conformation().set_torsion_angle(self.atoms[i], self.atoms[i+1], self.atoms[i+2], self.atoms[i+3], old)

    def __str__(self):
        # pretty sure we don't need this
        pass

    def apply(self, pose):
        d_angle = (np.random.rand()-0.5)*self.angle/180*np.pi

        conf = pose.conformation()
        dihe_start = np.random.randint(0, len(self.atoms)-3) # has to be able to move only bb atoms
        old = conf.torsion_angle(self.atoms[dihe_start], self.atoms[dihe_start+1], self.atoms[dihe_start+2], self.atoms[dihe_start+3])
        new = old + d_angle
        # print('Changing Torsion',dihe_start, 'from', old, 'to', new)
        conf.set_torsion_angle(self.atoms[dihe_start], self.atoms[dihe_start+1], self.atoms[dihe_start+2], self.atoms[dihe_start+3], new)


