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
    sc_values, reps = schema.items()
    sc_values = list(sc_values[1][float])
    sc_values.sort()
    select_sc_values = [str(round(a, 3)) for a in sc_values]

    # Creating new signac proejct
    project = signac.get_project()

    # Clustering parameters to scan over
    eps_params = np.linspace(0.05, 0.4, 15)
    min_samples = np.linspace(50, 350, 7, dtype = int)
    for sc in sc_values:
        for eps in eps_params:
            for min_sample in min_samples:
                exists = False
                for job in project.find_jobs({"eps":eps, "min_samples":min_sample, "sc_size":sc}):
                    exists = True
                if exists is False:
                    print("New job", "\nsc_size:", sc, "\neps", eps, "\nmin_samples", min_sample)
                    job = project.open_job({"eps":eps, "min_samples":min_sample, "sc_size":sc})
                    job.init()
                    sh.copyfile("cluster_parameters.yml", job.fn("cluster_parameters.yml"))
                else:
                    print("Skipping\n", "sc_size:", sc, "\neps", eps)
                    pass
if __name__ == "__main__":
    main()
