import signac
import numpy as np
import os
import shutil as sh

project = signac.get_project()

bond_angles = np.linspace(100, 160, 31)
reps = 500

for bond_angle in bond_angles:
    for rep in range(reps):
        job = project.open_job({'bond_angle': bond_angle, 'rep':rep})
        sh.copytree("parameters", os.path.join(job.workspace(), "parameters"))
        job.init()
