import argparse
import numpy as np
import sys
import getpass
from datetime import datetime
import platform
from itertools import combinations, compress, permutations
import mdtraj as md

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
                print(abc_id)
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


def parse_args():
    parser = argparse.ArgumentParser(
        description = "This script converts between Gromacs .itp files \
                        and Rosetta .param files."
    )
    parser.add_argument(
        "-f", "--file",
        type = str,
        help = ".itp file to convert to a .param file"
    )

    args = parser.parse_args()
    return args


def main():
    args = parse_args()
    itp_file = ItpFileObject(args.file, )
    id = args.file.split(".")[0]
    # itp_file.write_param_file("parameters/tph_tetramer.param", "em_terphenyl_pmp_octamer.gro")
    itp_file.write_atom_properties("parameters/atom_properties.txt")
    itp_file.write_mm_atom_properties("parameters/mm_atom_properties.txt")
    itp_file.write_mm_bond_lengths("parameters/mm_bond_length_params.txt")
    itp_file.write_mm_bond_angles("parameters/mm_bond_angle_params.txt")
    itp_file.write_mm_torsions("parameters/mm_torsion_params.txt")

    
    # Test some basic CG PyRosetta functionality
    
    import cg_pyrosetta

    cg_pyrosetta.init(add_atom_types = "fa_standard parameters/atom_properties.txt", 
                      add_mm_atom_type_set_parameters = "fa_standard parameters/mm_atom_properties.txt",
                      extra_res_fa = "parameters/tph_tetramer.param",
                      extra_mm_params_dir = "parameters",
                      mute = "no"
    )
    pymol = cg_pyrosetta.pyrosetta.PyMOLMover()
    pose = cg_pyrosetta.pyrosetta.pose_from_sequence("X[MOL]")
    pymol.apply(pose)
    pose.dump_pdb("test.pdb")
    small = cg_pyrosetta.CG_movers.CGSmallMover(pose)

    energy_function = cg_pyrosetta.CG_monte_carlo.EnergyFunctionFactory().build_energy_function(
        {
            "mm_twist" : 1,
            "mm_bend" : 1,
            "mm_stretch" : 1,
            "ring_close" : 10,
            "fa_atr" : 1,
            "fa_rep" : 1,
            "fa_intra_rep" : 1,
            "fa_intra_atr" : 1,
        }
    )
    
    print("Total Energy:", energy_function(pose))
    

    # Build Minimizer
    mini = cg_pyrosetta.pyrosetta.rosetta.protocols.minimization_packing.MinMover()
    mini.min_type('lbfgs_armijo_nonmonotone')
    mini.score_function(energy_function)

    movemap = cg_pyrosetta.pyrosetta.MoveMap()
    movemap.set_bb_true_range(1, pose.size())
    movemap.set_branches(True)
    movemap.set(cg_pyrosetta.pyrosetta.rosetta.core.id.PHI, True)
    movemap.set(cg_pyrosetta.pyrosetta.rosetta.core.id.THETA, True)

    mini.apply(pose)
    pymol.apply(pose)

    small.angle = 10

    for i in range(5000):
        small.apply(pose)
        pymol.apply(pose)
        mini.apply(pose)
        pymol.apply(pose)
        print("Step: ", i)
        print("Total Energy:", energy_function(pose))

    pose.dump_pdb("test.pdb")




if __name__ == "__main__":
    main()