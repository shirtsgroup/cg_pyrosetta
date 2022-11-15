import argparse
import mdtraj as md
import copy
import os

def parse_args():

    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--medoid_ids",
        required = True,
        help = "id of medoids to mirror",
        type = int,
        nargs = '+'
    )

    parser.add_argument(
        "--cluster_output_dir",
        required = False,
        default = "cluster_output"
    )
    
    return parser.parse_args()

def main():
    args = parse_args()
    for medoid_id in args.medoid_ids:
        medoid = md.load(os.path.join(args.cluster_output_dir, "medoid_" + str(medoid_id) + ".pdb"))
        mirror = copy.deepcopy(medoid)
        mirror.xyz[0, :, 0] = -mirror.xyz[0, :, 0]
        mirror.save(os.path.join(args.cluster_output_dir, "mirror_" + "medoid_" + str(medoid_id) + ".pdb"))

if __name__ == "__main__":
    main()

