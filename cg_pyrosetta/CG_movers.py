import numpy as np
import os
import sys
import warnings
import random

current_path = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.abspath(current_path + '/../PyRosetta4.modified'))

import pyrosetta


# General Movers
# Torsion, Bondangle and Bondlength movers that have all posible DOFs
# in a 

class CGSmallAngleMover(pyrosetta.rosetta.protocols.moves.Mover):
    """
    Small angle mover for moving all angles within a CG model. Will move a random angle
    by an angle uniformly distributed between `self.angle/2` and `-self.angle/2`.
    """
    def __init__(self, pose, angle = 5):
        """
        CG SmallAngleMover constructor

        Parameters
        ----------
        pose : pyrosetta.rosetta.protocols.moves.Mover
            Pose to apply this mover to. Mover needs the pose object to find 
            all the bond-angles present in the model
        angle : float (5)
            2x the max angle moves can be
        """
        pyrosetta.rosetta.protocols.moves.Mover.__init__(self)
        self.pose = pose
        self.angle = angle
        self.conf = pose.conformation()
        self.bond_angles = []
        self.atoms = [pyrosetta.AtomID(1, 1)]
        for atom_1 in self.atoms:
            # print("Working on", atom_1)
            for atom_2 in self.get_neighbors(atom_1):
                if self.is_new_atom(atom_2):
                    self.atoms.append(atom_2)
                for atom_3 in self.get_neighbors(atom_2):
                    # print("Angle Candidate:", atom_1, atom_2, atom_3)
                    if [atom_1.rsd(), atom_1.atomno()] == [atom_3.rsd(), atom_3.atomno()]:
                        continue
                    if pose.has_dof(self.conf.atom_tree().bond_angle_dof_id(atom_1, atom_2, atom_3, 0)):
                        if self.is_new_angle([atom_1, atom_2, atom_3]):
                            self.bond_angles.append([atom_1, atom_2, atom_3])
                            # print("Adding Bond Angle!")
                            # print("A1:", atom_1)
                            # print("A2:", atom_2)
                            # print("A3:", atom_3)
                        else:
                            continue
            

    def is_new_angle(self, angle):
        """
        Function which check if an angle is present in the mover's angle list

        Parameters
        ----------
        angle : list of pyrosetta.rosetta.core.id.AtomID
            List of atoms present in angle
        """
        for old_angle in self.bond_angles:
            if old_angle[0].rsd() == angle[0].rsd() and \
                old_angle[0].atomno() == angle[0].atomno() and \
                old_angle[1].rsd() == angle[1].rsd() and \
                old_angle[1].atomno() == angle[1].atomno() and \
                old_angle[2].rsd() == angle[2].rsd() and \
                old_angle[2].atomno() == angle[2].atomno():
                return(False)
            if old_angle[0].rsd() == angle[2].rsd() and \
                old_angle[0].atomno() == angle[2].atomno() and \
                old_angle[1].rsd() == angle[1].rsd() and \
                old_angle[1].atomno() == angle[1].atomno() and \
                old_angle[2].rsd() == angle[0].rsd() and \
                old_angle[2].atomno() == angle[0].atomno():
                return(False)
        return(True)

    def is_new_atom(self, atom):
        """
        Function which check if an atom is present in the mover's atom list

        Parameters
        ----------
        atom : pyrosetta.rosetta.core.id.AtomID
            AtomID to check
        """
        if [atom.rsd(), atom.atomno()] in [[old_atom.rsd(), old_atom.atomno()] for old_atom in self.atoms]:
            return(False)
        else:
            return(True)
    
    def get_neighbors(self, atom):
        """
        Get neighbors of an atom within the pose object

        Parameters
        ----------
        atom : pyrosetta.rosetta.core.id.AtomID
            AtomID to check
        """
        return(self.conf.bonded_neighbor_all_res(pyrosetta.AtomID(atom.atomno(), atom.rsd())))

    def apply(self, pose):
        """
        Apply mover to pose intially given to mover object
        """
        d_angle = (np.random.rand()-0.5)*self.angle/180*np.pi
        angle_i = np.random.randint(0, len(self.bond_angles))  # has to be able to move only bb atoms
        old = self.conf.bond_angle(*self.bond_angles[angle_i])
        new = old + d_angle
        # print('Changing Torsion',angle_start, 'from', old, 'to', new)
        self.conf.set_bond_angle(*self.bond_angles[angle_i], new)

class CGSmallMover(pyrosetta.rosetta.protocols.moves.Mover):
    """
    Small torsion mover for moving all angles within a CG model. Will move a random angle
    by an angle uniformly distributed between `self.angle/2` and `-self.angle/2`.
    """
    def __init__(self, pose, angle = 180):
        """
        CG SmallMover constructor

        Parameters
        ----------
        pose : pyrosetta.rosetta.protocols.moves.Mover
            Pose to apply this mover to. Mover needs the pose object to find 
            all the bond-angles present in the model
        angle : float (180)
            2x the max angle moves can be
        """
        pyrosetta.rosetta.protocols.moves.Mover.__init__(self)
        self.pose = pose
        self.angle = angle
        self.conf = pose.conformation()
        self.torsions = []
        self.atoms = [pyrosetta.AtomID(1, 1)]
        for atom_1 in self.atoms:
            # print("Working on", atom_1)
            for atom_2 in self.get_neighbors(atom_1):
                if self.is_new_atom(atom_2):
                    self.atoms.append(atom_2)
                for atom_3 in self.get_neighbors(atom_2):
                    if [atom_1.rsd(), atom_1.atomno()] == [atom_3.rsd(), atom_3.atomno()]:
                        continue
                    for atom_4 in self.get_neighbors(atom_3):
                        if [atom_2.rsd(), atom_2.atomno()] == [atom_4.rsd(), atom_4.atomno()]:
                            continue
                        if pose.has_dof(self.conf.atom_tree().torsion_angle_dof_id(atom_1, atom_2, atom_3, atom_4, 0)):
                            if self.is_new_torsion([atom_1, atom_2, atom_3, atom_4]):
                                self.torsions.append([atom_1, atom_2, atom_3, atom_4])
                            else:
                                continue
            

    def is_new_torsion(self, torsion):
        """
        Function which check if a torsion is present in the mover's torsion list

        Parameters
        ----------
        torsion : list of pyrosetta.rosetta.core.id.AtomID
            List of atoms present in torsion
        """
        for old_torsion in self.torsions:
            if old_torsion[0].rsd() == torsion[0].rsd() and \
                old_torsion[0].atomno() == torsion[0].atomno() and \
                old_torsion[1].rsd() == torsion[1].rsd() and \
                old_torsion[1].atomno() == torsion[1].atomno() and \
                old_torsion[2].rsd() == torsion[2].rsd() and \
                old_torsion[2].atomno() == torsion[2].atomno() and \
                old_torsion[3].atomno() == torsion[3].atomno() and \
                old_torsion[3].rsd() == torsion[3].rsd():
                return(False)
            if old_torsion[0].rsd() == torsion[3].rsd() and \
                old_torsion[0].atomno() == torsion[3].atomno() and \
                old_torsion[1].rsd() == torsion[2].rsd() and \
                old_torsion[1].atomno() == torsion[2].atomno() and \
                old_torsion[2].rsd() == torsion[1].rsd() and \
                old_torsion[2].atomno() == torsion[1].atomno() and \
                old_torsion[3].atomno() == torsion[0].atomno() and \
                old_torsion[3].rsd() == torsion[0].rsd():
                return(False)
        return(True)
    def is_new_atom(self, atom):
        """
        Function which check if an atom is present in the mover's atom list

        Parameters
        ----------
        atom : pyrosetta.rosetta.core.id.AtomID
            AtomID to check
        """
        if [atom.rsd(), atom.atomno()] in [[old_atom.rsd(), old_atom.atomno()] for old_atom in self.atoms]:
            return(False)
        else:
            return(True)
    
    def get_neighbors(self, atom):
        """
        Get neighbors of an atom within the pose object

        Parameters
        ----------
        atom : pyrosetta.rosetta.core.id.AtomID
            AtomID to check
        """
        return(self.conf.bonded_neighbor_all_res(pyrosetta.AtomID(atom.atomno(), atom.rsd())))

    def apply(self, pose):
        """
        Apply mover to pose intially given to mover object
        """
        d_angle = (np.random.rand()-0.5)*self.angle/180*np.pi
        angle_i = np.random.randint(0, len(self.torsions))  # has to be able to move only bb atoms
        old = self.conf.torsion_angle(*self.torsions[angle_i])
        new = old + d_angle
        # print('Changing Torsion',angle_start, 'from', old, 'to', new)
        self.conf.set_torsion_angle(*self.torsions[angle_i], new)

class CGBondLengthMover(pyrosetta.rosetta.protocols.moves.Mover):
    """
    Small bond-lenght mover for a CG pose model
    """
    def __init__(self, pose, d_bond = 0.1):
        """
        CG SmallBondLengthMover constructor

        Parameters
        ----------
        pose : pyrosetta.rosetta.protocols.moves.Mover
            Pose to apply this mover to. Mover will the pose object to find 
            all the bond-lenghts present in the model
        d_bond : float (0.1)
            2x the max distance bond-length moves can be
        """
        pyrosetta.rosetta.protocols.moves.Mover.__init__(self)
        self.pose = pose
        self.d_bond = d_bond
        self.conf = pose.conformation()
        self.bond_lengths = []
        self.atoms = [pyrosetta.AtomID(1, 1)]
        for atom_1 in self.atoms:
            for atom_2 in self.get_neighbors(atom_1):
                if self.is_new_atom(atom_2):
                    self.atoms.append(atom_2)
                if pose.has_dof(self.conf.atom_tree().bond_length_dof_id(atom_1, atom_2)):
                    if self.is_new_bond([atom_1, atom_2]):
                        self.bond_lengths.append([atom_1, atom_2])
                    else:
                        continue
            

    def is_new_bond(self, bond):
        for old_bond in self.bond_lengths:
            if old_bond[0].rsd() == bond[0].rsd() and \
                old_bond[0].atomno() == bond[0].atomno() and \
                old_bond[1].rsd() == bond[1].rsd() and \
                old_bond[1].atomno() == bond[1].atomno():
                return(False)
            if old_bond[0].rsd() == bond[1].rsd() and \
                old_bond[0].atomno() == bond[1].atomno() and \
                old_bond[1].rsd() == bond[0].rsd() and \
                old_bond[1].atomno() == bond[0].atomno():
                return(False)
        return(True)
    def is_new_atom(self, atom):
        if [atom.rsd(), atom.atomno()] in [[old_atom.rsd(), old_atom.atomno()] for old_atom in self.atoms]:
            return(False)
        else:
            return(True)
    
    def get_neighbors(self, atom):
        return(self.conf.bonded_neighbor_all_res(pyrosetta.AtomID(atom.atomno(), atom.rsd())))

    def apply(self, pose):
        """
        Apply mover to pose intially given to mover object
        """
        d_bond = (np.random.rand()-0.5)*self.d_bond
        bond_i = np.random.randint(0, len(self.bond_lengths))  # has to be able to move only bb atoms
        old = self.conf.bond_angle(*self.bond_lengths[bond_i])
        new = old + d_bond
        # print('Changing Torsion',angle_start, 'from', old, 'to', new)
        self.conf.set_bond_length(*self.bond_lengths[bond_i], new)

# Set Movers, these movers are used to uniformly change internal coordinates
# for a CG model. ex. change all bb dihedral angles to a value.
# randomize all bb dihedral angles
class randomizeAngles(CGSmallAngleMover):
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
        CGSmallAngleMover.__init__(self, pose)

    def apply(self, pose):
        """
        applies mover to a given pose

        Arguments
        ---------

        pose : pyrosetta.Pose()
            pose which you would like to apply randomizeBackbone mover to
        """
        for i in range(len(self.angles)):
            angle = np.random.uniform(-np.pi, np.pi)
            # print("Position :", i, ":", angle)
            self.conf.set_bond_angle(self.angles[i][0], self.angles[i][1], self.angles[i][2], angle)

class randomizeTorsions(CGSmallMover):
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
        CGSmallMover.__init__(self, pose)

    def apply(self, pose):
        """
        applies mover to a given pose

        Arguments
        ---------

        pose : pyrosetta.Pose()
            pose which you would like to apply randomizeBackbone mover to
        """
        for i in range(len(self.dihes)):
            angle = np.random.uniform(-np.pi, np.pi)
            # print("Position :", i, ":", angle)
            self.conf.set_torsion_angle(self.dihes[i][0], self.dihes[i][1], self.dihes[i][2], self.dihes[i][3], angle)

class setBondLengths(CGBondLengthMover):
    def __init__(self, pose, bond_length_dict):
        """
        Build setBackBone Mover for CG models

        Arguments
        ---------

        pose : pyrosetta.Pose()
            used to generate list of possible dihedrals used in randomizing pose
        bond_length_dicts : dict
            dictionary with bond names and their desired length in r_bb

        Example
        -------

        >>>pose = pyrosetta.pose_from_seqence('X[CG11]X[CG11]X[CG11]X[CG11]')
        >>>set_bond_length = cg_pyrosetta.CG_movers.setBondLengths(pose, {"BB1 SC1":sc, "BB2 SC2":sc, "BB3 SC3":sc})
        >>>set_bond_length.apply(pose)

        """
        CGBondLengthMover.__init__(self, pose)
        self.bond_length_dict = bond_length_dict
        self.bond_names = []
        for i in range(len(self.bond_lengths)):
            atom_1_name = pose.residue(self.bond_lengths[i][0].rsd()).atom_name(self.bond_lengths[i][0].atomno()).rstrip()
            atom_2_name = pose.residue(self.bond_lengths[i][1].rsd()).atom_name(self.bond_lengths[i][1].atomno()).rstrip()
            self.bond_names.append(atom_1_name+" "+atom_2_name)

    def apply(self, pose):
        """
        Apply bondlength changer to desired pose
        """
        conf = pose.conformation()
        for bond_atoms, bond_name in zip(self.bond_lengths, self.bond_names):
            # print(bond_name)
            # print(self.bond_length_dict.keys())
            rev_bond_name = bond_name.split(" ")
            rev_bond_name.reverse()
            rev_bond_name = " ".join(rev_bond_name)
            if bond_name in self.bond_length_dict.keys():
                # print(self.bond_length_dict[bond_name])
                # print("Changing", bond_name, "to a length of:", self.bond_length_dict[bond_name])
                conf.set_bond_length(bond_atoms[0], bond_atoms[1], self.bond_length_dict[bond_name])
            elif rev_bond_name in self.bond_length_dict.keys():
                # print(self.bond_length_dict[rev_bond_name])
                # print("Changing", rev_bond_name, "to a length of:", self.bond_length_dict[rev_bond_name])
                conf.set_bond_length(bond_atoms[1], bond_atoms[0], self.bond_length_dict[rev_bond_name])
                
            else:
                continue

class setBondAngle(CGSmallAngleMover):
    def __init__(self, pose, bond_angle_dict, radians = False):
        """
        Build setBackBone Mover for CG models

        Arguments
        ---------

        pose : pyrosetta.Pose()
            used to generate list of possible dihedrals used in randomizing pose
        bond_length_dicts : dict
            dictionary with bond names and their desired length in r_bb

        Example
        -------

        >>>pose = pyrosetta.pose_from_seqence('X[CG11]X[CG11]X[CG11]X[CG11]')
        >>>set_bond_angle = cg_pyrosetta.CG_movers.set_BondAngle(pose)
        >>>set_bond_angle.apply(pose)

        """
        CGSmallAngleMover.__init__(self, pose)
        self.bond_angle_dict = bond_angle_dict
        if radians is False:
            for key in self.bond_angle_dict.keys():
                self.bond_angle_dict[key] = self.bond_angle_dict[key] * np.pi / 180
        self.angle_names = []
        for i in range(len(self.bond_angles)):
            atom_1_name = pose.residue(self.bond_angles[i][0].rsd()).atom_name(self.bond_angles[i][0].atomno()).rstrip()
            atom_2_name = pose.residue(self.bond_angles[i][1].rsd()).atom_name(self.bond_angles[i][1].atomno()).rstrip()
            atom_3_name = pose.residue(self.bond_angles[i][2].rsd()).atom_name(self.bond_angles[i][2].atomno()).rstrip()
            self.angle_names.append(atom_1_name+" "+atom_2_name+" "+atom_3_name)

    def apply(self, pose):
        """
        Apply bondlength changer to desired pose
        """
        conf = pose.conformation()
        for angle_atoms, angle_name in zip(self.bond_angles, self.angle_names):
            print(angle_name)
            print(self.bond_angle_dict.keys())
            rev_angle_name = angle_name.split(" ")
            rev_angle_name.reverse()
            rev_angle_name = " ".join(rev_angle_name)
            if angle_name in self.bond_angle_dict.keys():
                print(self.bond_angle_dict[angle_name])
                print("Changing", angle_name, "to a length of:", self.bond_angle_dict[angle_name])
                conf.set_bond_angle(angle_atoms[0], angle_atoms[1], angle_atoms[2], self.bond_angle_dict[angle_name])
            elif rev_angle_name in self.bond_angle_dict.keys():
                print(self.bond_angle_dict[rev_angle_name])
                print("Changing", rev_angle_name, "to a length of:", self.bond_angle_dict[rev_angle_name])
                conf.set_bond_angle(angle_atoms[2], angle_atoms[1], angle_atoms[0], self.bond_angle_dict[rev_angle_name])
                
            else:
                continue

# Backbone/Sidechain specific movers

class CGBBSmallMover(CGSmallMover):
    """
    Variant of CGSmall mover which only moves backbone torsions
    """
    def __init__(self, pose, angle = 180):
        super().__init__(pose, angle)
        for torsion in self.torsions:
            is_bb = []
            for atom in torsion:
                is_bb.append(self.conf.atom_is_backbone_norefold(atom))
            if not all(is_bb):
                self.torsions.remove(torsion)

class CGSCSmallMover(CGSmallMover):
    """
    Variant of CGSmall mover which only moves torsions with SC atoms
    """
    def __init__(self, pose, angle = 180):
        super().__init__(pose, angle)
        for torsion in self.torsions:
            is_sc = []
            for atom in torsion:
                is_sc.append(not self.conf.atom_is_backbone_norefold(atom))
            if not any(is_sc):
                self.torsions.remove(torsion)

class CGBBSmallAngleMover(CGSmallAngleMover):
    """
    Variant of CGSmallAngleMover which only moves backbone angles
    """
    def __init__(self, pose, angle = 10):
        super().__init__(pose, angle)
        for angle in self.bond_angles:
            is_bb = []
            for atom in angle:
                is_bb.append(self.conf.atom_is_backbone_norefold(atom))
            if not all(is_bb):
                self.bond_angles.remove(angle)

class CGSCSmallAngleMover(CGSmallAngleMover):
    """
    Variant of CGSmallAngleMover which only moves angles with SC atoms
    """
    def __init__(self, pose, angle = 10):
        super().__init__(pose, angle)
        for angle in self.bond_angles:
            is_sc = []
            for atom in angle:
                is_sc.append(not self.conf.atom_is_backbone_norefold(atom))
            if not any(is_sc):
                self.bond_angles.remove(angle)
