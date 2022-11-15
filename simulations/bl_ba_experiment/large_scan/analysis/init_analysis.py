import signac
import numpy as np
import pandas as pd
import analyze_foldamers
import os
import yaml
import shutil as sh

def main():
    stream = open("cluster_parameters.yml", "r")
    params = yaml.safe_load(stream)

    # Getting Signac project schema
    print("Fetching signac schema information...")
    original_project = signac.get_project(root=params["signac_dir"])
    schema = original_project.detect_schema()
    bond_angles, reps = schema.items()
    bond_angles = list(bond_angles[1][float])
    bond_angles.sort()
    select_bond_angles = [str(round(a, 3)) for a in bond_angles]

    # Creating new signac proejct
    project = signac.get_project()

    # Clustering parameters to scan over
    eps_params = np.linspace(0.05, 0.5, 30)

    for ba in bond_angles:
        for eps in eps_params:
            exists = False
            for job in project.find_jobs({"eps":eps, "bond_angle":ba}):
                exists = True
            if exists is False:
                print("New job\n", "bond_angle:", ba, "\neps", eps)
                job = project.open_job({"bond_angle":ba, "eps":eps})
                job.init()
            else:
                print("Skipping\n", "bond_angle:", ba, "\neps", eps)
                pass
if __name__ == "__main__":
    main()