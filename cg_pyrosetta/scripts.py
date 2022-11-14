import argparse
import os
from .parameters import ItpFileObject


def itp_to_rosetta_params():
    def parse_args():
        parser = argparse.ArgumentParser(
            description = "This script converts between Gromacs .itp files \
                            and Rosetta .param files."
        )
        parser.add_argument(
            "-f", "--file",
            type = str,
            help = ".itp file to convert to a .param file",
            required = True
        )
        parser.add_argument(
            "-s", "--structure",
            type = str,
            help = "Structure file (.pdb, .gro) to build structure from",
            required = True
        )
        parser.add_argument(
            "-o", "--output",
            type = str,
            help = "output name for residue file",
            required = True
        )
        parser.add_argument(\
            "-p", "--prefix",
            type = str,
            help = "prefix for filenames. Default parameters",
            default = "parameters",
            nargs = "?"
        )

        args = parser.parse_args()
        return args
    args = parse_args()
    itp_file = ItpFileObject(args.file)
    id = args.file.split(".")[0]
    print(args.prefix)
    if not os.path.exists(args.prefix):
        os.makedirs(args.prefix)
    itp_file.write_param_file(args.prefix + "/" + args.output + ".param", args.structure)
    itp_file.write_atom_properties(args.prefix + "/" + "atom_properties.txt")
    itp_file.write_mm_atom_properties(args.prefix + "/" + "mm_atom_properties.txt")
    itp_file.write_mm_bond_lengths(args.prefix + "/" + "mm_bond_length_params.txt")
    itp_file.write_mm_bond_angles(args.prefix + "/" + "mm_bond_angle_params.txt")
    itp_file.write_mm_torsions(args.prefix + "/" + "mm_torsion_params.txt")