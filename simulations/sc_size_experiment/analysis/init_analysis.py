import signac
import numpy as np
import pandas as pd
import analyze_foldamers
import os
import yaml

def main():
    stream = open("cluster_parameters.yml", "r")
    params = yaml.safe_load(stream)

    # Getting Signac project schema
    print("Fetching signac schema information...")
    original_project = signac.get_project(root=params["signac_dir"])
    schema = original_project.detect_schema()
    sc_sizes, reps = schema.items()
    sc_sizes = list(sc_sizes[1][float])
    sc_sizes.sort()
    select_sc_sizes = [str(round(a, 3)) for a in sc_sizes]

    # Creating new signac proejct
    project = signac.get_project()

    for sc in sc_sizes:
        print("New job:", sc)
        job = project.open_job({'sc_size':sc})
        job.init()

if __name__ == "__main__":
    main()