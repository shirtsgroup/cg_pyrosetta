import signac
import numpy as np
import pandas as pd
import analyze_foldamers
import os
import argparse

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--signac_dir",
        help="Directory containing signac project",
        required=False,
        type=str,
        default=os.path.abspath(""),
    )
    
    args = parser.parse_args()
    return args

def main():
    args = parse_args()
    original_path = os.path.abspath("")

    # Getting Signac project schema
    print("Fetching signac schema information...")
    os.chdir(args.signac_dir)
    original_project = signac.get_project(root=args.signac_dir)
    schema = original_project.detect_schema()
    sc_sizes, reps = schema.items()
    sc_sizes = list(sc_sizes[1][float])
    sc_sizes.sort()
    select_sc_sizes = [str(round(a, 3)) for a in sc_sizes]

    # Creating new signac proejct
    os.chdir(original_path)
    project = signac.get_project()

    for sc in sc_sizes:
        print("New job:", sc)
        job = project.open_job({'sc_size':sc, 'production_path':args.signac_dir})
        job.init()

if __name__ == "__main__":
    main()