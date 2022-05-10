import argparse
import mdtraj as md
import copy
import os

def parse_args():

    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--file_paths",
        required = True,
        help = "id of medoids to mirror",
        type = str,
        nargs = '+'
    )
    
    return parser.parse_args()

def main():
    args = parse_args()
    for f in args.file_paths:
        medoid = md.load(f)
        mirror = copy.deepcopy(medoid)
        mirror.xyz[0, :, 0] = -mirror.xyz[0, :, 0]
        file_array = f.split("/")
        if len(file_array) > 1:
            mirror.save("/".join(file_array[:-1]) + "/mirror_" + file_array[-1])
        else:
            mirror.save("mirror_" + file_array[0])

if __name__ == "__main__":
    main()

