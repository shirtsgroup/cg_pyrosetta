import argparse
import numpy as np
import sys
import getpass
from datetime import datetime
import platform
from itertools import combinations, compress, permutations

class ItpFileObject:
    def __init__(self, filename):
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
        atomtype_index = 0
        while "[ atomtypes ]" not in self.itp_file[atomtype_index]:
            atomtype_index += 1
        
        end_section = atomtype_index + 1
        while "[" not in self.itp_file[end_section]:
            end_section += 1
            if len(self.itp_file) == end_section:
                break

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
                    "type" : line_entries[1],
                    "resnr" : int(line_entries[2]),
                    "residue" : line_entries[3],
                    "atom" : line_entries[4],
                    "cgnr" : int(line_entries[5]),
                    "charge" : float(line_entries[6]),
                    "mass" : float(line_entries[7]),
                    }
                )

    def get_bonds(self):
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
    
    def is_bonded(self, ai, aj):
        check_bond = self.find_bond(ai, aj)
        if check_bond is None:
            return False
        else:
            return True

    def get_neighbors(self, ai):
        n_filter = [self.is_bonded(ai, n + 1) or self.is_bonded(n + 1, ai) for n in range(len(self.atoms))]
        neighbors_filterd = compress(self.atoms, n_filter)
        return [atom["nr"] for atom in neighbors_filterd]


    def find_angle(self, ai, aj, ak):
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

    
    def write_param_file(self, param_filename):
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
                param_f.write("ATOM " + atom["atom"] + " " + atom["type"] + " " + atom["type"] + " " + str(atom["charge"]) + "\n")
            param_f.write("\n")

            # Bond section
            for bond in self.bonds:
                param_f.write("BOND " + self.atoms[bond["ai"]-1]["atom"] + " " + self.atoms[bond["aj"]-1]["atom"] + "\n")
            param_f.write("\n")

            # Z-matrix entry 1
            z_entry = ["ICOOR_INTERNAL",
                       self.atoms[0]["atom"],
                       0,
                       0,
                       0,
                       self.atoms[0]["atom"],
                       self.atoms[1]["atom"],
                       self.atoms[2]["atom"]
                    ]

            param_f.write('{0:<13} {1:^8} {2:>10} {3:>10} {4:>10} {5:>6} {6:>6} {7:>6}\n'.format(*z_entry))
            
            # Z-matrix entry 2
            z_entry = ["ICOOR_INTERNAL",
                       self.atoms[1]["atom"],
                       0,
                       round(180,5),
                       round(10*self.find_bond(1, 2)["length"], 5),
                       self.atoms[0]["atom"],
                       self.atoms[1]["atom"],
                       self.atoms[2]["atom"]
                    ]

            param_f.write('{0:<13} {1:^8} {2:>10} {3:>10} {4:>10} {5:>6} {6:>6} {7:>6}\n'.format(*z_entry))


            # Z-matrix entry 3
            z_entry = ["ICOOR_INTERNAL",
                       self.atoms[2]["atom"],
                       0,
                       round(self.find_angle(1,2,3)["angle"], 5),
                       round(10*self.find_bond(2,3)["length"], 5),
                       self.atoms[0]["atom"],
                       self.atoms[1]["atom"],
                       self.atoms[2]["atom"]
                    ]

            param_f.write('{0:<13} {1:^8} {2:>10} {3:>10} {4:>10} {5:>6} {6:>6} {7:>6}\n'.format(*z_entry))

            # Writing remainder Z-matrix entries
            existing_atoms = [self.atoms[0]["nr"],
                              self.atoms[1]["nr"],
                              self.atoms[2]["nr"]
                        ]

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
                                                
                                    z_entry = ["ICOOR_INTERNAL",
                                            self.atoms[torsion[3]-1]["atom"],
                                            round(self.find_torsion(*torsion)["angle"], 5),
                                            round(self.find_angle(*b_angle)["angle"], 5),
                                            round(10*self.find_bond(*bond_length)["length"], 5),
                                            self.atoms[torsion[2]-1]["atom"],
                                            self.atoms[torsion[1]-1]["atom"],
                                            self.atoms[torsion[0]-1]["atom"]
                                            ]

                                    param_f.write('{0:<13} {1:^8} {2:>10} {3:>10} {4:>10} {5:>6} {6:>6} {7:>6}\n'.format(*z_entry))

    def write_atom_properties(self, param_filename):
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
                if atom_type in at_in_atoms:
                    entry = [
                        atom_type,
                        ''.join([i for i in atom_type if not i.isdigit()]),  # remove numbers from atomtype
                        round(self.atomtypes[atom_type]["sigma"] * 10 * 2 ** (1/6) / 2, 5) , # Converting nm to Anstroms and r_min and divide by 2 for Rosetta conventions
                        round(self.atomtypes[atom_type]["epsilon"], 5),
                        round(1.0, 5),                                   # Placeholders since we don't use these
                        round(3.5, 5),                                   # Placeholders since we don't use these
                        round(13.5, 5)
                    ]
                    atom_f.write('{0:<5}{1:>5}{2:>10}{3:>10}{4:>10}{5:>10}{6:>10}\n'.format(*entry))

    def write_mm_atom_properties(self, mm_param_filename):
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
                if atom_type in at_in_atoms:
                    entry = [
                        atom_type,
                        -round(self.atomtypes[atom_type]["epsilon"], 5),                        
                        round(self.atomtypes[atom_type]["sigma"] * 10 * 2 ** (1/6) / 2, 5) , # Converting nm to Anstroms and r_min and divide by 2 for Rosetta conventions
                        -round(self.atomtypes[atom_type]["epsilon"], 5),                        
                        round(self.atomtypes[atom_type]["sigma"] * 10 * 2 ** (1/6) / 2, 5) , # Converting nm to Anstroms and r_min and divide by 2 for Rosetta conventions
                    ]
                    atom_f.write('{0:<5}{1:>5}{2:>10}{3:>10}{4:>10}\n'.format(*entry))




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
    itp_file = ItpFileObject(args.file)
    id = args.file.split(".")[0]
    itp_file.write_param_file(id + ".params")
    itp_file.write_atom_properties("atom_properties.txt")
    itp_file.write_mm_atom_properties("mm_atom_properties.txt")

    import cg_pyrosetta
    cg_pyrosetta.init(add_atom_types = "fa_standard atom_properties.txt",
                      add_mm_atom_type_set_parameters = "fa_standard mm_atom_properties.txt",
                      extra_res_fa = "OCT.params"
                )
    cg_pyrosetta.pyrosetta.pose_from_sequence("X[OCT]")





if __name__ == "__main__":
    main()