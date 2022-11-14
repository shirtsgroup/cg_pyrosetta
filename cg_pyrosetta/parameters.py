import shutil
import os
import re
import argparse
import numpy as np
import sys
import getpass
from datetime import datetime
import platform
from itertools import combinations, compress, permutations
import mdtraj as md
import cg_pyrosetta.build_cg_pyrosetta

current_path = os.path.dirname(os.path.abspath(__file__))
data_path = os.path.join(current_path, 'data')


def changeAtomParameters(param_dict, atom_types_path = None, mm_atom_types_path = None):
    """
    function to change atom parameters on the fly

    Arguments
    ---------

    dict : dict
        Dictionary containing name of atom and list of specific parameters for
        that atom type. (Parameters are expected in the following order
        ATOM LJ_RADIUS LJ_WDEPTH LK_DGFREE LK_LAMBDA LK_VOLUME)

    Examples
    --------

    params = {'CG1':['X', 1.0, 0.2, 1.0, 3.5, 23.7]}
    cg_pyrosetta.parameters.changeAtomParameters(params)

    """
    if atom_types_path is None:
        atom_types_path = os.path.join(data_path, 'atom_type_sets', 'atom_properties.txt')
    
    if mm_atom_types_path is None:
        mm_atom_types_path = os.path.join(data_path, 'mm_atom_type_sets', 'mm_atom_properties.txt')


    with open(atom_types_path, 'r') as f:
        atom_lines = f.readlines()

    atom_params_list = [line.rstrip('\n').split() for line in atom_lines[1:]]
    prev_atoms = [atom_list[0] for atom_list in atom_params_list]

    names = param_dict.keys()

    for name in names:
        assert len(param_dict[name]) >= 3

        # extract LJ parameters and atom name
        atom_type = param_dict[name][0]
        lj_radius = param_dict[name][1]
        lj_wdepth = param_dict[name][2]

        # can take either a list of length 3 or 6

        # if LK parameters are provided, will use what is provided
        if len(param_dict[name]) == 6:
            lk_dgfree = param_dict[name][3]
            lk_lambda = param_dict[name][4]
            lk_volume = param_dict[name][5]
        # otherwise use default values
        else:
            lk_dgfree = 1.0
            lk_lambda = 3.5
            lk_volume = 23.7

        if name in prev_atoms:
            i_atom = prev_atoms.index(name)
            atom_params_list[i_atom][0] = name
            atom_params_list[i_atom][1] = atom_type
            atom_params_list[i_atom][2] = lj_radius
            atom_params_list[i_atom][3] = lj_wdepth
            atom_params_list[i_atom][4] = lk_dgfree
            atom_params_list[i_atom][5] = lk_lambda
            atom_params_list[i_atom][6] = lk_volume
        else:
            new_entry = [name,
                         atom_type,
                         lj_radius,
                         lj_wdepth,
                         lk_dgfree,
                         lk_lambda,
                         lk_volume]
            atom_params_list.append(new_entry)

    # Rewrite parameters to atom_properties.txt file
    with open(atom_types_path, 'w') as f:
        # write header
        f.write('NAME  ATOM LJ_RADIUS LJ_WDEPTH LK_DGFREE LK_LAMBDA LK_VOLUME\n')

        # rewrite updated parameters
        for param in atom_params_list:
            # print param line with specific formating
            f.write("%s%6.1s%10.4f%10.4f%10.4f%10.4f%10.4f\n" % (param[0], param[1], float(
                param[2]), float(param[3]), float(param[4]), float(param[5]), float(param[6])))

    with open(mm_atom_types_path, 'w') as f:
        # write header
        f.write('NAME    LJ_WDEPTH   LJ_RADIUS   LJ_3B_WDEPTH    LJ_3B_RADIUS\n')

        # rewrite updated parameters
        for param in atom_params_list:
            # print param line with specific formating
            f.write("%s%10.4f%10.4f%10.4f%10.4f\n" %
                    (param[0], -float(param[3]), float(param[2]), -float(param[3]), float(param[2])))
    # cg_pyrosetta.builder.buildCGPyRosetta()


def changeTorsionParameters(param_dict, torsion_file = None, mode = "single"):
    """
    function to change torsion parameters on the fly

    Arguments
    ---------

    dict : dict
        Dictionary containing name of torsion and list of specific parameters for
        that atom type. (Parameters are expected in the following order
        ATOM LJ_RADIUS LJ_WDEPTH LK_DGFREE LK_LAMBDA LK_VOLUME)
    torsion_file : string
        Path to file to change. If left blank will default to data files within CG PyRosetta package
    mode : string
        whether to insert single or multiple instances of torsion parameters (to create complex torsion patterns)

    Examples
    --------

    params = {'CG1 CG1 CG1 CG1':[k_force, periodicity, phase_shift]}
    cg_pyrosetta.parameters.changeAtomParameters(params)

    """

    if torsion_file is None:
        torsion_file = os.path.join(data_path, 'mm_atom_type_sets', 'mm_torsion_params.txt')


    with open(torsion_file, 'r') as f:
        torsion_lines = f.readlines()

    # Build torsion params list [name_of_torsion k_constant periodicity phase_shift]
    torsion_params_list = []
    for line in torsion_lines[1:]:
        split_line = line.rstrip('\n').split()
        name = ' '.join(split_line[:4])
        torsion_params_list.append([name, split_line[4], split_line[5], split_line[6]])

    prev_torsion = [torsion[0] for torsion in torsion_params_list]
    names = param_dict.keys()

    for name in names:

        # extract LJ parameters and atom name
        force_constant = param_dict[name][0]
        periodicity = param_dict[name][1]
        phase_shift = param_dict[name][2]

        if name in prev_torsion:
            i_torsion = prev_torsion.index(name)
            torsion_params_list[i_torsion][0] = name
            torsion_params_list[i_torsion][1] = force_constant
            torsion_params_list[i_torsion][2] = periodicity
            torsion_params_list[i_torsion][3] = phase_shift

        else:
            new_entry = [name,
                         force_constant,
                         periodicity,
                         phase_shift]
            torsion_params_list.append(new_entry)

    # Rewrite parameters to atom_properties.txt file
    with open(torsion_file, 'w') as f:
        # write header
        f.write('# CG torsion parameters\n')

        # rewrite updated parameters
        for param in torsion_params_list:
            # print param line with specific formating
            atom_names = param[0].split()
            f.write("%s %s %s %s %.4f %.1i %.4f\n" % (
                atom_names[0], atom_names[1], atom_names[2], atom_names[3], float(param[1]), int(param[2]), float(param[3])))
    # cg_pyrosetta.builder.buildCGPyRosetta()


def changeAngleParameters(param_dict, angle_file = None):
    """
    function to change angle parameters (mm_bend)on the fly

    Arguments
    ---------

    dict : dict
        Dictionary containing name of atom and list of specific parameters for
        that atom type. (Parameters are expected in the following order
        !atom types     Ktheta    Theta0   Kub     S0)

    Examples
    --------

    params = {'CG1 CG1 CG1':[10, 120]}
     cg_pyrosetta.parameters.changeAtomParameters(params)

    """

    if angle_file is None:
        angle_file = os.path.join(data_path, 'mm_atom_type_sets', 'mm_angle_params.txt')

    with open(angle_file, 'r') as f:
        angle_lines = f.readlines()

    # Build torsion params list [name_of_torsion k_constant periodicity phase_shift]
    angle_params_list = []
    for line in angle_lines[2:]:
        split_line = line.rstrip('\n').split()
        name = ' '.join(split_line[:3])
        angle_params_list.append([name, split_line[3], split_line[4]])

    prev_angles = [angle[0] for angle in angle_params_list]
    names = param_dict.keys()
    for name in names:

        # extract LJ parameters and atom name
        force_constant = param_dict[name][0]
        angle = param_dict[name][1]

        if name in prev_angles:
            i_angle = prev_angles.index(name)
            angle_params_list[i_angle][0] = name
            angle_params_list[i_angle][1] = force_constant
            angle_params_list[i_angle][2] = angle

        else:
            new_entry = [name,
                         force_constant,
                         angle, ]
            angle_params_list.append(new_entry)

    # Rewrite parameters to atom_properties.txt file
    with open(angle_file, 'w') as f:
        # write header
        f.write('ANGLES\n')
        f.write('!atom types     Ktheta    Theta0   Kub     S0\n')

        # rewrite updated parameters
        for param in angle_params_list:
            # print param line with specific formating
            atom_names = param[0].split()
            f.write("%-4.4s%4.3s%4.3s%10.4f%10.4f\n" %
                    (atom_names[0], atom_names[1], atom_names[2], float(param[1]), float(param[2])))
    
    # cg_pyrosetta.builder.buildCGPyRosetta()

class ItpFileObject:
    """
    Object to convert between gromacs itp files and rosetta parameter files
    """
    def __init__(self, filename):
        """
        Instantiate an ItpFileObject

        Parameters
        ----------
        filename : str
            string of the filename to load into the object
        """
        self.filename = filename
        self.itp_file = []
        with open(self.filename, 'r') as f:
            for line in f:
                if line[0] != ";":
                    self.itp_file.append(line)
                else:
                    pass
        self.get_atomtypes()
        self.get_atoms()
        self.get_bonds()
        self.get_angles()
        self.get_torsions()
        self.get_molecule()

    def get_molecule(self):
        """
        Function for saving molecule name from itp file
        """
        moltype_index = 0
        while "[ moleculetype ]" not in self.itp_file[moltype_index]:
            moltype_index += 1
        
        end_section = moltype_index + 1
        while "[" not in self.itp_file[end_section]:
            end_section += 1
            if len(self.itp_file) == end_section:
                break

        name_line = moltype_index + 1
        line_entries = self.itp_file[name_line].split()
        self.mol_name = line_entries[0]
    
    def get_atomtypes(self):
        """
        Function for saving atomtypes from itp file
        """
        atomtype_index = 0
        while "[ atomtypes ]" not in self.itp_file[atomtype_index]:
            atomtype_index += 1
        
        end_section = atomtype_index + 1
        while "[" not in self.itp_file[end_section]:
            if len(self.itp_file) == end_section:
                break
            else:
                end_section += 1
            

        self.atomtypes = {}
        for i in range(atomtype_index + 1, end_section):
            if len(self.itp_file[i]) > 1:
                line = self.itp_file[i]
                line_entries = line.split() # split line based on whitespace 
                self.atomtypes[line_entries[0]] = {
                    "atom number" : int(line_entries[1]),
                    "mass" : float(line_entries[2]),
                    "charge" : float(line_entries[3]),
                    "ptype" : line_entries[4],
                    "sigma" : float(line_entries[5]),
                    "epsilon" : float(line_entries[6]),
                }

    def get_atoms(self):
        """
        Function for saving atom names and other properties from itp file
        """
        atom_index = 0
        while "[ atoms ]" not in self.itp_file[atom_index]:
            atom_index += 1
        
        end_section = atom_index + 1
        while "[" not in self.itp_file[end_section]:
            end_section += 1
            if len(self.itp_file) == end_section:
                break

        self.atoms = []

        for i in range(atom_index + 1, end_section):
            if len(self.itp_file[i]) > 1:
                line = self.itp_file[i]
                line_entries = line.split() # split line based on whitespace 
                self.atoms.append(
                        {
                    "nr" : int(line_entries[0]),
                    "type" : line_entries[1].lower(),
                    "resnr" : int(line_entries[2]),
                    "residue" : line_entries[3],
                    "atom" : line_entries[4],
                    "cgnr" : int(line_entries[5]),
                    "charge" : float(line_entries[6]),
                    "mass" : float(line_entries[7]),
                    }
                )
        
        # Rosetta only reads in at most 3 character atoms...
        if len(self.atoms) > 100:
            characters = 'abcdefghijklmnopqrstuvwxyz'
            for i in range(len(self.atoms)):
                element = self.atoms[i]["atom"][0]
                number = int(self.atoms[i]["atom"][1:]) - 1
                abc_id = element + characters[int(np.floor(number/26))] + characters[number % 26]
                self.atoms[i]["atom"] = abc_id


    def get_bonds(self):
        """
        Function for saving bond information from itp file
        """
        bond_index = 0
        while "[ bonds ]" not in self.itp_file[bond_index]:
            bond_index += 1
        
        end_section = bond_index + 1
        while "[" not in self.itp_file[end_section]:
            end_section += 1
            if len(self.itp_file) == end_section:
                break

        self.bonds = []
        for i in range(bond_index + 1, end_section):
            if len(self.itp_file[i]) > 1:
                line = self.itp_file[i]
                line_entries = line.split() # split line based on whitespace 
                self.bonds.append(
                        {
                    "ai" : int(line_entries[0]),
                    "aj" : int(line_entries[1]),
                    "funct" : int(line_entries[2]),
                    "length" : float(line_entries[3]),
                    "k" : float(line_entries[4]),
                    }
                )
    def get_angles(self):
        """
        Function for saving atom information from itp file
        """
        angle_index = 0
        while "[ angles ]" not in self.itp_file[angle_index]:
            angle_index += 1
        
        end_section = angle_index + 1
        while "[" not in self.itp_file[end_section]:
            end_section += 1
            if len(self.itp_file) == end_section:
                break

        self.angles = []
        for i in range(angle_index + 1, end_section):
            if len(self.itp_file[i]) > 1:
                line = self.itp_file[i]
                line_entries = line.split() # split line based on whitespace 
                self.angles.append(
                        {
                    "ai" : int(line_entries[0]),
                    "aj" : int(line_entries[1]),
                    "ak" : int(line_entries[2]),
                    "funct" : int(line_entries[3]),
                    "angle" : float(line_entries[4]),
                    "k" : float(line_entries[5]),
                    }
                )

    def get_torsions(self):
        """
        Function for saving torsion information from itp file
        """
        torsion_index = 0
        while "[ dihedrals ]" not in self.itp_file[torsion_index]:
            torsion_index += 1
        
        end_section = torsion_index + 1
        while "[" not in self.itp_file[end_section]:
            end_section += 1
            if len(self.itp_file) == end_section:
                break


        self.torsions = []
        for i in range(torsion_index + 1, end_section):
            if len(self.itp_file[i]) > 1:
                line = self.itp_file[i]
                line_entries = line.split() # split line based on whitespace 
                self.torsions.append(
                        {
                    "ai" : int(line_entries[0]),
                    "aj" : int(line_entries[1]),
                    "ak" : int(line_entries[2]),
                    "al" : int(line_entries[3]),
                    "funct" : int(line_entries[4]),
                    "angle" : float(line_entries[5]),
                    "k" : float(line_entries[6]),
                    "periodicity" : int(line_entries[7]),
                    }
                )

    def find_bond(self, ai, aj):
        """
        Function for identifying a bond entry by atom type

        Parameters
        ----------
        ai : int
            atom index of atom i
        aj : int
            atom index of atom j

        Returns
        -------
        bond : dict
            bond information dictionary. If no bond is found, returns None
        """
        for perm in permutations([ai, aj]):
            for bond in self.bonds:
                if bond["ai"] == perm[0]:
                    if bond["aj"] == perm[1]:
                        return bond
                    else:
                        pass
                else:
                    pass
        return None
    
    def is_hydrogen(self, atom):
        """
        Function to check if a provided atom index is a hydrogen

        Parameters
        ----------
        atom : int
            atom index to check
        """
        return''.join([i.lower() for i in atom["type"] if not i.isdigit()]) == "h"
    
    def is_bonded(self, ai, aj):
        """
        Function to check if two provided atom indices are bonded

        Parameters
        ----------
        ai : int
            atom index i
        aj : int
            atom index j

        Returns
        -------
        out : bool
        """
        check_bond = self.find_bond(ai, aj)
        if check_bond is None:
            return False
        else:
            return True

    def get_neighbors(self, ai):
        """
        Function returns list of neighboring atoms

        Parameters
        ----------
        ai : int
            index of atom
        """
        n_filter = [self.is_bonded(ai, n + 1) or self.is_bonded(n + 1, ai) for n in range(len(self.atoms))]
        neighbors_filterd = compress(self.atoms, n_filter)
        return [atom["nr"] for atom in neighbors_filterd]


    def find_angle(self, ai, aj, ak):
        """
        Function for identifying angle information by atom number. Checks for
        permutations of the provided atoms.

        Parameters
        ----------
        ai : int
            atom index of atom i
        aj : int
            atom index of atom j
        ak : int
            atom index of atom k

        Returns
        -------
        angle : dict
            angle information dictionary. If no bond is found, returns None
        """
        for perm in permutations([ai, aj, ak]):
            for angle in self.angles:
                if angle["ai"] == perm[0]:
                    if angle["aj"] == perm[1]:
                        if angle["ak"] == perm[2]:
                            return angle
                        else:
                            pass
                    else:
                        pass
                else:
                    pass
        return None

    def find_torsion(self, ai, aj, ak, al):
        """
        Function for identifying torsion information by atom number. Checks for permutations
        of the provided atoms.

        Parameters
        ----------
        ai : int
            atom index of atom i
        aj : int
            atom index of atom j
        ak : int
            atom index of atom k
        al : int
            atom index of atom l

        Returns
        -------
        angle : dict
            torsion information dictionary. If no bond is found, returns None
        """
        for perm in permutations([ai, aj, ak, al]):
            for torsion in self.torsions:
                if torsion["ai"] == perm[0]:
                    if torsion["aj"] == perm[1]:
                        if torsion["ak"] == perm[2]:
                            if torsion["al"] == perm[3]:
                                return torsion
                            else:
                                pass
                        else:
                            pass
                    else:
                        pass
                else:
                    pass
        return None

    
    def write_param_file(self, param_filename, structure_file, nbr_radius = 20):
        """
        Function to write Rosetta parameter file

        Parameters
        ----------
        param_filename : str
            File name to write parameter file
        structure_file : str
            Structure file used for placing atoms in 3D space
        nbr_radius : float, default 20
        """

        md_structure = md.load(structure_file)

        with open(param_filename, "w") as param_f:
            # write header
            param_f.write("# Rosetta residue topology file\n")
            param_f.write("# Generated using itp_to_param.py\n")
            param_f.write("# Original file: " + self.filename + "\n")
            param_f.write("# Author: " + getpass.getuser() + "\n")
            param_f.write("# Date: " + datetime.now().strftime("%A, %d. %B %Y") + "\n")
            param_f.write("# Time: " + datetime.now().strftime("%I:%M%p")+ "\n")
            param_f.write("# System: " + platform.platform() + "\n\n")

            # Name section
            param_f.write("NAME " + self.mol_name + "\n")
            param_f.write("IO_STRING " + self.mol_name + " X\n")
            param_f.write("TYPE POLYMER\n")
            param_f.write("AA UNK\n\n")

            # Atom section
            for atom in self.atoms:
                entry = [
                    "ATOM", atom["atom"], atom["type"], atom["type"],atom["charge"], 5
                ]
                param_f.write('{0:<5} {1:<5} {2:<5} {3:<5} {4:10.5f}\n'.format(*entry))
            param_f.write("\n")

            # Bond section
            for bond in self.bonds:
                entry = [
                    "BOND", self.atoms[bond["ai"]-1]["atom"], self.atoms[bond["aj"]-1]["atom"]
                ]
                param_f.write('{0:<5} {1:<5} {2:<5}\n'.format(*entry))
            param_f.write("\n")

            # non-bonded atom
            heavy_atoms = []
            for atom in self.atoms:
                if self.is_hydrogen(atom):
                    continue
                else:
                    heavy_atoms.append(atom)
            
            middle_heavy_atoms = int(len(heavy_atoms)/2)
            param_f.write("NBR_ATOM " + heavy_atoms[middle_heavy_atoms]["atom"] + "\n")
            param_f.write("NBR_RADIUS " + str(nbr_radius) + "\n\n")

            # Find first 3 heavy atoms that are bonded

            start_atoms = []

            # Get the first heavy atom
            for atom in self.atoms:
                if self.is_hydrogen(atom): # checks if Hydrogen
                    continue
                else:
                    start_atoms.append(atom["nr"])
                    break
            
            # Find 2 more heavy atoms bonded to

            for atom_2 in self.get_neighbors(start_atoms[0]):
                if self.is_hydrogen(self.atoms[atom_2-1]):
                    continue
                else:
                    start_atoms.append(atom_2)
                    break
            
            for atom_3 in self.get_neighbors(start_atoms[1]):
                if self.is_hydrogen(self.atoms[atom_3-1]):
                    continue
                else:
                    if atom_3 != start_atoms[0]:
                        start_atoms.append(atom_3)
                        break
                    else:
                        continue
            
            # Z-matrix entry 1
            z_entry = ["ICOOR_INTERNAL",
                       self.atoms[start_atoms[0]-1]["atom"],
                       0,
                       0,
                       0,
                       self.atoms[start_atoms[0]-1]["atom"],
                       self.atoms[start_atoms[1]-1]["atom"],
                       self.atoms[start_atoms[2]-1]["atom"]
                    ]

            param_f.write('{0:<13} {1:>8} {2:10.5f} {3:10.5f} {4:10.5f} {5:>6} {6:>6} {7:>6}\n'.format(*z_entry))
            
            # Z-matrix entry 2
            
            z_entry = ["ICOOR_INTERNAL",
                       self.atoms[start_atoms[1]-1]["atom"],
                       0,
                       round(180,5),
                       round(10*self.find_bond(start_atoms[0], start_atoms[1])["length"], 5),
                       self.atoms[start_atoms[0]-1]["atom"],
                       self.atoms[start_atoms[1]-1]["atom"],
                       self.atoms[start_atoms[2]-1]["atom"]
                    ]

            param_f.write('{0:<13} {1:>8} {2:10.5f} {3:10.5f} {4:10.5f} {5:>6} {6:>6} {7:>6}\n'.format(*z_entry))


            # Z-matrix entry 3
            z_entry = ["ICOOR_INTERNAL",
                       self.atoms[start_atoms[2]-1]["atom"],
                       0,
                       round(180 - self.find_angle(start_atoms[0], start_atoms[1], start_atoms[2])["angle"], 5),
                       round(10*self.find_bond(start_atoms[1], start_atoms[2])["length"], 5),
                       self.atoms[start_atoms[1]-1]["atom"],
                       self.atoms[start_atoms[0]-1]["atom"],
                       self.atoms[start_atoms[2]-1]["atom"]
                    ]
            param_f.write('{0:<13} {1:>8} {2:10.5f} {3:10.5f} {4:10.5f} {5:>6} {6:>6} {7:>6}\n'.format(*z_entry))

            # Writing remainder Z-matrix entries
            existing_atoms = [*start_atoms]

            # Grow a torsion list from existing atoms
            torsion_list = []
            atom_tree_iter = 0
            print("Building atom tree...")
            while len(existing_atoms) != len(self.atoms):
                atom_tree_iter += 1
                print("Iteration:", atom_tree_iter)
                print("existing atoms:", len(existing_atoms))
                print("self.atoms:", len(self.atoms))
                for atom_1 in existing_atoms:
                    for atom_2 in self.get_neighbors(atom_1):
                        # print(self.get_neighbors(atom_1))
                        if atom_2 not in existing_atoms:
                            continue
                        for atom_3 in self.get_neighbors(atom_2):
                            if atom_3 not in existing_atoms:
                                continue
                            if atom_1 == atom_3:
                                continue
                            for atom_4 in self.get_neighbors(atom_3):
                                potential_torsion = [atom_1, atom_2, atom_3, atom_4]
                                if atom_4 not in existing_atoms:
                                    existing_atoms.append(atom_4)
                                else:
                                    continue
                                if atom_2 == atom_4 or \
                                    potential_torsion in torsion_list or \
                                    potential_torsion[::-1] in torsion_list:  
                                    continue
                                else:
                                    if self.find_torsion(*potential_torsion) is not None:
                                        torsion = potential_torsion
                                        b_angle = torsion[1:]
                                        bond_length = torsion[2:]
                                        
                                    elif self.find_torsion(*potential_torsion[::-1]) is not None:
                                        torsion = potential_torsion[::-1]
                                        b_angle = torsion[1:]
                                        bond_length = torsion[2:]
                                    else:
                                        sys.exit("Error! Torsion between atoms" + ", ".join(potential_torsion) + "was not found in the .itp file.")
                                    
                                    # calculate torsion/bond angle/bond length from structure                                    
                                    torsion_index = [t_i - 1 for t_i in torsion]
                                    b_angle_index = [ba_i - 1 for ba_i in b_angle]
                                    bond_length_index = [bl_i - 1 for bl_i in bond_length]
                                    structure_torsion = md.compute_dihedrals(md_structure, np.array([*torsion_index]).reshape(1,4), periodic=False)[0][0]
                                    structure_torsion *= 180 / np.pi
                                    structure_angle = md.compute_angles(md_structure, np.array([*b_angle_index]).reshape(1,3), periodic=False)[0][0]
                                    structure_angle *= 180 / np.pi
                                    structure_bl = md.compute_distances(md_structure, np.array([*bond_length_index]).reshape(1,2), periodic=False)[0][0]
                                    
                                    
                                    z_entry = ["ICOOR_INTERNAL",
                                            self.atoms[torsion_index[3]]["atom"],
                                            round(structure_torsion, 5),
                                            round(180 - structure_angle, 5),
                                            round(10*structure_bl, 5),
                                            self.atoms[torsion_index[2]]["atom"],
                                            self.atoms[torsion_index[1]]["atom"],
                                            self.atoms[torsion_index[0]]["atom"]
                                            ]

                                    param_f.write('{0:<13} {1:>8} {2:10.5f} {3:10.5f} {4:10.5f} {5:>6} {6:>6} {7:>6}\n'.format(*z_entry))

    def write_atom_properties(self, param_filename):
        """
        Function to write rosetta atom_propeties file

        Parameters
        ----------
        param_filename : str
            File name to write atom_property file
        """
        with open(param_filename, "w") as atom_f:
            # write header
            atom_f.write("NAME  ATOM LJ_RADIUS LJ_WDEPTH LK_DGFREE LK_LAMBDA LK_VOLUME\n")
            atom_f.write("# Rosetta atom_properties file\n")
            atom_f.write("# Generated using itp_to_param.py\n")
            atom_f.write("# Original file: " + self.filename + "\n")
            atom_f.write("# Author: " + getpass.getuser() + "\n")
            atom_f.write("# Date: " + datetime.now().strftime("%A, %d. %B %Y") + "\n")
            atom_f.write("# Time: " + datetime.now().strftime("%I:%M%p")+ "\n")
            atom_f.write("# System: " + platform.platform() + "\n")
            at_in_atoms = list(set([a["type"] for a in self.atoms])) 
            for atom_type in self.atomtypes.keys():
                if atom_type.lower() in at_in_atoms:
                    entry = [
                        atom_type.lower(),
                        ''.join([i for i in atom_type if not i.isdigit()]),  # remove numbers from atomtype
                        round(self.atomtypes[atom_type]["sigma"] * 10 * 2 ** (1/6) / 2, 5) , # Converting nm to Anstroms and r_min and divide by 2 for Rosetta conventions
                        round(self.atomtypes[atom_type]["epsilon"], 5),
                        round(1.0, 5),                                   # Placeholders since we don't use these
                        round(3.5, 5),                                   # Placeholders since we don't use these
                        round(13.5, 5)
                    ]
                    atom_f.write('{0:<5}{1:>5}{2:>10}{3:>10}{4:>10}{5:>10}{6:>10}\n'.format(*entry))

    def write_mm_atom_properties(self, mm_param_filename):
        """
        Function to write rosetta mm_atom_properties file

        Parameters
        ----------
        mm_param_filename : str
            File name to write mm_atom_property file
        """
        with open(mm_param_filename, "w") as atom_f:
            # write header
            atom_f.write("NAME    LJ_WDEPTH   LJ_RADIUS   LJ_3B_WDEPTH    LJ_3B_RADIUS\n")
            atom_f.write("# Rosetta atom_properties file\n")
            atom_f.write("# Generated using itp_to_param.py\n")
            atom_f.write("# Original file: " + self.filename + "\n")
            atom_f.write("# Author: " + getpass.getuser() + "\n")
            atom_f.write("# Date: " + datetime.now().strftime("%A, %d. %B %Y") + "\n")
            atom_f.write("# Time: " + datetime.now().strftime("%I:%M%p")+ "\n")
            atom_f.write("# System: " + platform.platform() + "\n")
            at_in_atoms = list(set([a["type"] for a in self.atoms])) 
            for atom_type in self.atomtypes.keys():
                if atom_type.lower() in at_in_atoms:
                    entry = [
                        atom_type.lower(),
                        -round(self.atomtypes[atom_type]["epsilon"], 5),                        
                        round(self.atomtypes[atom_type]["sigma"] * 10 * 2 ** (1/6) / 2, 5) , # Converting nm to Anstroms and r_min and divide by 2 for Rosetta conventions
                        -round(self.atomtypes[atom_type]["epsilon"], 5),                        
                        round(self.atomtypes[atom_type]["sigma"] * 10 * 2 ** (1/6) / 2, 5) , # Converting nm to Anstroms and r_min and divide by 2 for Rosetta conventions
                    ]
                    atom_f.write('{0:<5}{1:>5}{2:>10}{3:>10}{4:>10}\n'.format(*entry))



    def write_mm_bond_lengths(self, mm_bond_length_filename):
        """
        Function to write rosetta mm_bond_length parameter file

        Parameters
        ----------
        mm_bond_length_filename : str
            File name to write mm_bond_length file
        """
        with open(mm_bond_length_filename, "w") as bl_f:
            # write header
            bl_f.write("ANGLES\n")
            bl_f.write("! Rosetta atom_properties file\n")
            bl_f.write("! Generated using itp_to_param.py\n")
            bl_f.write("! Original file: " + self.filename + "\n")
            bl_f.write("! Author: " + getpass.getuser() + "\n")
            bl_f.write("! Date: " + datetime.now().strftime("%A, %d. %B %Y") + "\n")
            bl_f.write("! Time: " + datetime.now().strftime("%I:%M%p")+ "\n")
            bl_f.write("! System: " + platform.platform() + "\n")
            bl_f.write("! atom types     Kb    b0\n")

            # write angle parameters 
            exisiting_angles = {}
            for bond in self.bonds:
                type_i, type_j = self.atoms[bond["ai"]-1]["type"], self.atoms[bond["aj"]-1]["type"]
                bond_id = type_i + type_j
                if bond_id not in exisiting_angles.keys():
                    entry = [
                        type_i,
                        type_j,
                        round(bond["k"] * 2 * 0.239 / 100, 5), # kJ to kcal and nm to A
                        round(bond["length"]*10, 5)
                    ]
                    bl_f.write('{0:<5}{1:<5}{2:>10}{3:>10}\n'.format(*entry))
                    exisiting_angles[bond_id] = [bond["k"], bond["length"]]
                else:
                    print("Angle Type:", bond_id, "was already found. Skipping!")
                    assert(exisiting_angles[bond_id][0] == bond["k"])
                    assert(exisiting_angles[bond_id][1] == bond["angle"])
                    continue


    def write_mm_bond_angles(self, mm_bond_angle_filename):
        """
        Function to write rosetta mm_bond_angle parameter file

        Parameters
        ----------
        mm_bond_angle_filename : str
            File name to write mm_bond_angle file
        """
        with open(mm_bond_angle_filename, "w") as angle_f:
            # write header
            angle_f.write("ANGLES\n")
            angle_f.write("! Rosetta atom_properties file\n")
            angle_f.write("! Generated using itp_to_param.py\n")
            angle_f.write("! Original file: " + self.filename + "\n")
            angle_f.write("! Author: " + getpass.getuser() + "\n")
            angle_f.write("! Date: " + datetime.now().strftime("%A, %d. %B %Y") + "\n")
            angle_f.write("! Time: " + datetime.now().strftime("%I:%M%p")+ "\n")
            angle_f.write("! System: " + platform.platform() + "\n")
            angle_f.write("! atom types     Ktheta    Theta0   Kub     S0\n")

            # write angle parameters 
            exisiting_angles = {}
            for angle in self.angles:
                type_i, type_j, type_k = self.atoms[angle["ai"]-1]["type"], self.atoms[angle["aj"]-1]["type"], self.atoms[angle["ak"]-1]["type"]
                angle_id = type_i + type_j + type_k
                if angle_id not in exisiting_angles.keys():
                    entry = [
                        type_i,
                        type_j,
                        type_k,
                        round(angle["k"] * 2 * 0.239, 5), # kJ to kcal
                        round(angle["angle"], 4)
                    ]
                    angle_f.write('{0:<5}{1:<5}{2:<5}{3:>10}{4:>10}\n'.format(*entry))
                    exisiting_angles[angle_id] = [angle["k"], angle["angle"]]
                else:
                    print("Angle Type:", angle_id, "was already found. Skipping!")
                    assert(exisiting_angles[angle_id][0] == angle["k"])
                    assert(exisiting_angles[angle_id][1] == angle["angle"])
                    continue



    def write_mm_torsions(self, mm_torsion_filename):
        """
        Function to write rosetta mm_bond_torsion parameter file

        Parameters
        ----------
        mm_torsion_filename : str
            File name to write mm_bond_torsion file
        """
        with open(mm_torsion_filename, "w") as torsion_f:
            # write header
            torsion_f.write("# Rosetta atom_properties file\n")
            torsion_f.write("# Generated using itp_to_param.py\n")
            torsion_f.write("# Original file: " + self.filename + "\n")
            torsion_f.write("# Author: " + getpass.getuser() + "\n")
            torsion_f.write("# Date: " + datetime.now().strftime("%A, %d. %B %Y") + "\n")
            torsion_f.write("# Time: " + datetime.now().strftime("%I:%M%p")+ "\n")
            torsion_f.write("# System: " + platform.platform() + "\n")
            torsion_f.write("# atom types     Ktheta  periodicity  Theta0\n")

            # write angle parameters 
            exisiting_torsions = {}
            error_count = 0
            for torsion in self.torsions:
                type_i, type_j, type_k, type_l = self.atoms[torsion["ai"]-1]["type"], \
                                                self.atoms[torsion["aj"]-1]["type"], \
                                                self.atoms[torsion["ak"]-1]["type"], \
                                                self.atoms[torsion["al"]-1]["type"]
                torsion_id = type_i + type_j + type_k + type_l + "_" + str(torsion["periodicity"]) 
                if torsion_id not in exisiting_torsions.keys():
                    entry = [
                        type_i,
                        type_j,
                        type_k,
                        type_l,
                        round(torsion["k"] * 2 *0.239, 5), # kJ to kcal
                        torsion["periodicity"],
                        round(torsion["angle"], 5)
                    ]
                    torsion_f.write('{0:<5}{1:<5}{2:<5}{3:<5}{4:>5}{5:>10}{6:>10}\n'.format(*entry))
                    exisiting_torsions[torsion_id] = [torsion["k"], torsion["angle"]]
                else:
                    print("Torsion Type:", torsion_id, "was already found. Skipping!")
                    print("k:", exisiting_torsions[torsion_id][0], torsion["k"] )
                    print("torsion:", exisiting_torsions[torsion_id][1], torsion["angle"] )
                    if exisiting_torsions[torsion_id][0] != torsion["k"] or \
                            exisiting_torsions[torsion_id][1] != torsion["angle"]:
                        print("ERROR!!", torsion_id)
                        error_count += 1
                    continue
            
            print("There were", error_count, "instances of oversubscribed torsion ids")