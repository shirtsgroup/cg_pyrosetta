import cg_pyrosetta
import pyrosetta
import argparse
import yaml
import sys
import numpy as np

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--model",
                        required=False,
                        help="Type of model that you would like to generate \
                             random distribution of structures and energy for",
                        type=str,
                        default="CG11x3",
                        )
    parser.add_argument("--mer",
                        required=False,
                        help = "Number of monomers for this model",
                        type = int,
                        default = 5,
                        )
    parser.add_argument("--kt",
                        required=True,
                        default=100,
                        help="kT value used to generate a random ensemble",
                        type=int,
                        )
    parser.add_argument("--params",
                        required=False,
                        type=str,
                        help="YAML file where specific torsions are defined",
                        default="params.yml")

    parser.add_argument("--e_output",
                        required=False,
                        type=str,
                        default="outputs/random_distribution.npy",
                        help="output filename of the energies",
                        )
    parser.add_argument("--s_output",
                        required=False,
                        type=str,
                        default="outputs/random_distribution.pdb",
                        help="output filename of the energies",
                        )

    parser.add_argument("--steps",
                        required=False,
                        default=1000000,
                        type=int,
                        help="Number of steps to run the high T simulation",
                        )

    parser.add_argument("--stride",
                        required=False,
                        default=1000,
                        type=int,
                        help="Number of steps to run the high T simulation",
                        )

    # Parse Arguments
    args = parser.parse_args()
    
    # Get parameters from .yml file
    param_file = open(args.params, 'r')
    params = yaml.load(param_file)
    print(params)

    # for param_type in params:
    #     if param_type == "atoms":
    #         cg_pyrosetta.change_parameters.changeAtomParameters(params["atoms"])
    #     if param_type == "dihedrals":
    #         cg_pyrosetta.change_parameters.changeTorsionParameters(params["dihedrals"])
    #     if param_type == "angles":
    #         cg_pyrosetta.change_parameters.changeAngleParameters(params["angles"])
    #     else:
    #         print("Input YAML file had a key for", param_type+".", "This key is not a valid parameter type for this model.", file=sys.stderr)
    #         print("Ignoring", param_type, "and continuing!", file=sys.stderr)
    
    cg_pyrosetta.init()
    
    # Build sequence for folder object
    monomer = "X["+args.model+"]"
    sequence = monomer*args.mer

    # Build Folder Object
    folder = cg_pyrosetta.CG_folding.CGFoldingAlgorithm(sequence)
    
    # Change PDBWriter to specific outputs
    folder.PDB_writer.stride(args.stride)
    folder.PDB_writer.file_name(args.s_output)
    
    # Build MC algorithm to get unfolded ensemble
    folder.build_fold_alg('no_min')
    folder.add_folding_move('no_min', pyrosetta.RepeatMover(folder.small, 10))
    folder.add_folding_move('no_min', folder.PDB_writer)

    # Set kT
    folder.kT = args.kt

    # Run algorithm 
    folder.run_folding_alg('no_min', args.steps)


    energy_list = []
    for line in open(args.s_output):
        if "pose" in line:
            # print(line)
            termwise_energy = line.rstrip().split(' ')
            energy_list.append(float(termwise_energy[1]))

        else:
            continue

    np.save(args.e_output, np.array(energy_list))


    


if __name__ == "__main__":
    main()
