import signac
import numpy as np
import os
import shutil as sh

project = signac.get_project()

bond_lengths = np.linspace(0.5, 4, 10)
bond_angles = np.linspace(45, 135, 10)
reps = 100

for bond_length in bond_lengths:
    for bond_angle in bond_angles:
        for rep in range(reps):
            job = project.open_job({'bond_length':bond_length, 'bond_angle': bond_angle, 'rep':rep})
            sh.copytree("parameters", os.path.join(job.workspace(), "parameters"))
            job.init()
