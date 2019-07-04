import os
import sys
import warnings

current_path = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.abspath(current_path + '/../PyRosetta4.modified'))

import pyrosetta
import numpy as np

# Random Movers, these movers will randomly change an internal coordinate
# of a CG model

class CGSmallMoverOld(pyrosetta.rosetta.protocols.moves.Mover):
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
        pose: pyrosetta.Pose()
            used to generate list of possible atoms
        """
        self.angle = angle
        pyrosetta.rosetta.protocols.moves.Mover.__init__(self)
        # Generate a list of all possible dihe in a provided pose with the bb_model given
        #  1 bb model [[1 1] [2 1] [3 1] [4 1]][]
        #  2 bb model [[1 1] [1 2] [2 1] [2 2]][ [3 1] [3 2]]
        #  3 bb model [[1 1] [1 2] [1 3] [2 1] [2 2] [2 3] [3 1] [3 2] [3 3]]
        self.conf = pose.conformation()

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

        # initializing conformation 
        for i in range(0, len(self.dihes)):
            old = self.conf.torsion_angle(self.dihes[i][0], self.dihes[i][1], self.dihes[i][2], self.dihes[i][3])
            pose.conformation().set_torsion_angle(self.dihes[i][0], self.dihes[i][1], self.dihes[i][2], self.dihes[i][3], old)

    def __str__(self):
        # pretty sure we don't need this
        pass

    def apply(self, pose):
        d_angle = (np.random.rand()-0.5)*self.angle/180*np.pi
        dihe_start = np.random.randint(0, len(self.dihes)) # has to be able to move only bb atoms
        old = self.conf.torsion_angle(self.dihes[dihe_start][0], self.dihes[dihe_start][1], self.dihes[dihe_start][2], self.dihes[dihe_start][3])
        new = old + d_angle
        # print('Changing Torsion',dihe_start, 'from', old, 'to', new)
        self.conf.set_torsion_angle(self.dihes[dihe_start][0], self.dihes[dihe_start][1], self.dihes[dihe_start][2], self.dihes[dihe_start][3], new)

class CGSmallMover(pyrosetta.rosetta.protocols.moves.Mover):
    """
    Geberalized CG Mover
    """
    def __init__(self, pose, angle = 180):
        """
        Build a Small Mover for CG polymer backbone dihedrals

        Arguments
        ---------

        angle : float
            maximum angle the mover can change an angle
        pose : pyrosetta.Pose()
            used to generate list of possible atoms
        """

        self.angle = angle
        pyrosetta.rosetta.protocols.moves.Mover.__init__(self)
        # Generate a list of all possible backbone dihes in a provided pose

        self.conf = pose.conformation()

        # will populate these lists wil all backbone dihedrals
        self.bb_atoms = []
        self.dihes = []


        # Will build list of backbones by iterating through neighbors

        for i in range(1, pose.size()+1):
            
            # Start backbone search with fist atom of each residue
            start_atom = pyrosetta.AtomID(1, i)
            self.bb_atoms.append(start_atom)

            # Will keep track of residue bb atoms and perviously seen atoms

            residue_bb = [start_atom]
        
            for atom in residue_bb:
                neighbors = self.conf.bonded_neighbor_all_res(atom, virt=False)
                for n in neighbors:
                    info  = [n.rsd(), n.atomno()]
                    isBackbone = self.conf.atom_is_backbone_norefold(info[0], info[1])

                    # exluding atoms from other residues and previously seen atoms
                    if info in [[prev_atom.rsd(), prev_atom.atomno()] for prev_atom in residue_bb] or info[0] != i:  #skip atom if already visted or not in residue
                        continue

                    # Skip atoms not in backbone
                    if isBackbone == False:
                        continue
                    
                    else:
                        residue_bb.append(n)
                        self.bb_atoms.append(n)

        bb_atoms_info = [[bb_atom.rsd(), bb_atom.atomno()] for bb_atom in self.bb_atoms]


            # Check to see if bb atom is first atom and will try and 
            # add the SC atom to include the dihedral between bb atom 1 and bb atom 2
            # Shown below as ABCD
            #       A        0        0
            #       |        |        |
            #       B--C--D--0--0--0--0

            # Get neighbors of first atom (only have to do this once)
            
        end_cases = [self.bb_atoms[0], self.bb_atoms[-1]]
        pos = [0 , len(self.bb_atoms)+1]
        for atom, p in zip(end_cases, pos):
            neighs = self.conf.bonded_neighbor_all_res(atom, virt=False)
            init_add_sc = None
            
            # Iterate through neighbors and position we want to insert SC atoms
            for n in neighs:
                # If in bb atoms skip
                if [n.rsd(), n.atomno()] in bb_atoms_info:
                    continue
                # If add_sc already defined, skip
                elif init_add_sc != None:
                    continue
                # Add first sc atom to pop up
                else:
                    init_add_sc = n
                    self.bb_atoms.insert(p, init_add_sc)

        # Will now create a list of list of atoms for each dihedral

        for i in range(len(self.bb_atoms)-3):
            dihe = [self.bb_atoms[i], self.bb_atoms[i+1], self.bb_atoms[i+2], self.bb_atoms[i+3]]
            self.dihes.append(dihe)


    def __str__(self):
        # pretty sure we don't need this
        pass

    def apply(self, pose):
        d_angle = (np.random.rand()-0.5)*self.angle/180*np.pi
        dihe_start = np.random.randint(0, len(self.dihes)) # has to be able to move only bb atoms
        old = self.conf.torsion_angle(self.dihes[dihe_start][0], self.dihes[dihe_start][1], self.dihes[dihe_start][2], self.dihes[dihe_start][3])
        new = old + d_angle
        # print('Changing Torsion',dihe_start, 'from', old, 'to', new)
        self.conf.set_torsion_angle(self.dihes[dihe_start][0], self.dihes[dihe_start][1], self.dihes[dihe_start][2], self.dihes[dihe_start][3], new)  

class CGShearMover(CGSmallMover):
    """
    Generalized CG mover analogous to the "shear" mover in PyRosetta.
    """
    def __init__(self, pose, angle = 180):
        """
        Builds a Shear Mover for CG set of polymers

        Arguments
        ---------

        angle : float
            maximum angle the mover can change an angle
        bb_model: int
            describes how many backbone beads are in a given CG model
        pose: pyrosetta.Pose()
            used to generate list of possible atoms
        """
        CGSmallMover.__init__(self, pose, angle)

    def __str__(self):
        # pretty sure we don't need this
        pass

    def apply(self, pose):
        d_angle = (np.random.rand()-0.5)*self.angle/180*np.pi

        conf = pose.conformation()
        dihe_start = np.random.randint(0, len(self.dihes)-1) # has to be able to move only bb atoms
        old_1 = conf.torsion_angle(self.dihes[dihe_start][0], self.dihes[dihe_start][1], self.dihes[dihe_start][2], self.dihes[dihe_start][3])
        old_2 = conf.torsion_angle(self.dihes[dihe_start+1][0], self.dihes[dihe_start+1][1], self.dihes[dihe_start+1][2], self.dihes[dihe_start+1][3])
        new_1 = old_1 + d_angle
        new_2 = old_2 - d_angle
        # print('Changing Torsion',dihe_start, 'from', old, 'to', new)
        conf.set_torsion_angle(self.dihes[dihe_start][0], self.dihes[dihe_start][1], self.dihes[dihe_start][2], self.dihes[dihe_start][3], new_1)
        conf.set_torsion_angle(self.dihes[dihe_start+1][0], self.dihes[dihe_start+1][1], self.dihes[dihe_start+1][2], self.dihes[dihe_start+1][3], new_2)

class CGSmallSCMover(pyrosetta.rosetta.protocols.moves.Mover):
    """
    Generalized SC torsion mover 
    """

    def __init__(self, pose, angle = 30):
        """
        Build Small Mover for CG polymers sidechains

        Arguments
        ---------

        angle : float
            maximum angle the mover can change an angle
        pose : pyrosetta.Pose()
             used to generate list of possible atoms
        """
        self.angle = angle
        pyrosetta.rosetta.protocols.moves.Mover.__init__(self)
        # Generate a list of all possible dihe in a provided pose with the bb_model given


        # All of the following logic depends on the conformation object

        conf = pose.conformation()
        
        
        # will populate this list with all side chain dihedrals
        self.dihes = []


        # Take a residue based approach
        for i in range(1, pose.size()+1):  # should be genearlizeable to any combination of residues
            # Start each search for dihedrals at the first atom of each CG residue
            start_atom = pyrosetta.AtomID(1, i)


            # Getting 2 bb atoms used for building SC dihedral
            bb1 = start_atom
            bb2 = None

            # Will determine 
            residue_atoms = [start_atom]
            sidechain_atoms = []
            # get side chain atoms in residue and define bb2  for SC dihedral definition

            for atom in residue_atoms:
                neighbors = conf.bonded_neighbor_all_res(atom, virt=False)
                for neigh in neighbors:
                    info = [neigh.rsd(), neigh.atomno()]
                    isBackbone = conf.atom_is_backbone_norefold(info[0], info[1])

                    # extracting 2nd backbone atom
                    if isBackbone and bb2 == None:
                        bb2 = neigh
                    
                    # exluding atoms from other residues and previously seen atoms
                    if info in [[prev_atom.rsd(), prev_atom.atomno()] for prev_atom in residue_atoms] or info[0] != i:  #skip atom if already visted or not in residue
                        continue

                    # atoms that haven't been seen before are added to the residue_atoms list and side chain atoms are added to the sidechain_atoms list
                    else:
                        residue_atoms.append(neigh)
                        if not isBackbone:
                            sidechain_atoms.append(neigh)
            
            # generating sidechain dihedrals from sidechain atoms
            
            for sc in sidechain_atoms:
                # Build dihedral candidates from SC atoms
                potential_dihe = [sc]
                pos = sc

                # Recusively populate potential_dihe by moving down SC
                # If we reach backbone, appends the two backbone atoms describe above 
                 
                while len(potential_dihe) < 4:
                    neighbors = conf.bonded_neighbor_all_res(pos, virt=False)

                    # Next atom will always be that with the lowest atomno in the residue
                    atom_ids = [atom.atomno() for atom in neighbors]
                    i_next = atom_ids.index(min(atom_ids))

                    # If atom is in the backbone, append bb1 and bb2 to the potential dihe
                    if conf.atom_is_backbone_norefold(i, atom_ids[i_next]):
                        potential_dihe.append(bb1)
                        potential_dihe.append(bb2)
                        break

                    # Otherwise add the neighbor w/ lowest atomno and change the current pos to that atom    
                    else:
                        pos = neighbors[i_next+1]
                        potential_dihe.append(pos)
                

                # Quality control : Only potential_dihes over length 4 are considered and 
                # only use the first 4 entries of the dihedral. 
                # (Ex. This would be generated for CG13 model --> [SC3 SC2 SC1 BB1 BB2] )
                if len(potential_dihe) >= 4:
                    self.dihes.append(potential_dihe[0:4])



                                        

        #exit()

    def __str__(self):
        # pretty sure we don't need this
        pass

    def apply(self, pose):
        if self.dihes:
            d_angle = (np.random.rand()-0.5)*self.angle/180*np.pi
            conf = pose.conformation()
            dihe = np.random.randint(0, len(self.dihes)) # has to be able to move only bb atoms
            old = conf.torsion_angle(self.dihes[dihe][0], self.dihes[dihe][1], self.dihes[dihe][2], self.dihes[dihe][3])
            new = old + d_angle
            # print('Changing Torsion',dihe_start, 'from', old, 'to', new)
            conf.set_torsion_angle(self.dihes[dihe][0], self.dihes[dihe][1], self.dihes[dihe][2], self.dihes[dihe][3], new)
        else:
            warnings.warn('No Sidechain Dihedrals Available')


class CGSmallAngleMover(CGSmallMover):
    """
    Generalized CG mover analogous to the "small" mover in PyRosetta.
    """
    def __init__(self, pose, angle = 10):
        """
        Build Small Mover for CG polymers specifically for moving backbone bond angles

        Arguments
        ---------
        pose: pyrosetta.Pose()
            used to generate list of possible atoms

        angle : float
            maximum angle the mover can change an angle


        """
        CGSmallMover.__init__(self, pose, angle)
        self.bb_atoms = self.bb_atoms[1:-1]
        self.angles = []

        for i in range(len(self.bb_atoms)-2):
            angle = [self.bb_atoms[i], self.bb_atoms[i+1], self.bb_atoms[i+2]]
            self.angles.append(angle)



    def __str__(self):
        # pretty sure we don't need this
        pass

    def apply(self, pose):
        d_angle = (np.random.rand()-0.5)*self.angle/180*np.pi
        angle_start = np.random.randint(0, len(self.angles)) # has to be able to move only bb atoms
        old = self.conf.bond_angle(self.angles[angle_start][0], self.angles[angle_start][1], self.angles[angle_start][2])
        new = old + d_angle
        # print('Changing Torsion',angle_start, 'from', old, 'to', new)
        self.conf.set_bond_angle(self.angles[angle_start][0], self.angles[angle_start][1], self.angles[angle_start][2], new)


# Set Movers, these movers are used to uniformly change internal coordinates
# for a CG model. ex. change all bb dihedral angles to a value.
#                     randomize all bb dihedral angles

class randomizeBackBone(CGSmallMover):
    """
    Generalized CG model dihedral randomizer. Used for random initial starting configurations
    """
    def __init__(self, pose):
        """
        Build randomizeBackBone Mover for CG models

        Arguments
        ---------

        pose: pyrosetta.Pose()
            used to generate list of possible dihedrals used in randomizing pose
        

        Example
        -------

        >>>pose = pyrosetta.pose_from_seqence('X[CG11]X[CG11]X[CG11]X[CG11]')
        >>>randomizer = cg_pyrosetta.CG_movers.randomizeBackBone(pose)
        >>>randomizer.apply(pose)

        """
        # inherits CGSmallMover's init call, but runs with a different apply call
        CGSmallMover.__init__(self, pose)

    def apply(self, pose):
        for i in range(len(self.dihes)):
            angle = np.random.uniform(-180, 180)
            self.conf.set_torsion_angle(self.dihes[i][0], self.dihes[i][1], self.dihes[i][2], self.dihes[i][3], angle)
            
class randomizeBackBoneAngles(CGSmallAngleMover):
    """
    Generalized CG model dihedral randomizer. Used for random initial starting configurations
    """
    def __init__(self, pose):
        """
        Build randomizeBackBoneAngles Mover for CG models

        Arguments
        ---------

        pose: pyrosetta.Pose()
            used to generate list of possible dihedrals used in randomizing pose
        

        Example
        -------

        >>>pose = pyrosetta.pose_from_seqence('X[CG11]X[CG11]X[CG11]X[CG11]')
        >>>randomizer = cg_pyrosetta.CG_movers.randomizeBackBoneAngles(pose)
        >>>randomizer.apply(pose)

        """