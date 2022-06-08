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
    eps_params = np.linspace(0.05, 0.15, 5)
    min_samples = np.linspace(50, 350, 7, dtype = int)
    for ba in bond_angles:
        for eps in eps_params:
            for min_sample in min_samples:
                exists = False
                for job in project.find_jobs({"eps":eps, "min_samples":min_sample, "bond_angle":ba}):
                    exists = True
                if exists is False:
                    print("New job", "\nbond_angle:", ba, "\neps", eps, "\nmin_samples", min_sample)
                    job = project.open_job({"eps":eps, "min_samples":min_sample, "bond_angle":ba})
                    job.init()
                    sh.copyfile("cluster_parameters.yml", job.fn("cluster_parameters.yml"))
                else:
                    print("Skipping\n", "bond_angle:", ba, "\neps", eps)
                    pass
if __name__ == "__main__":
    main()
