import signac
import numpy as np
import os
import shutil as sh

project = signac.get_project()
bond_angles = [60, 120, 180]
bond_lengths = [1.0, 1.1, 1.2, 1.3, 1.4, 1.5]

for rep in range(100):
    for bl in bond_lengths:
        for angle in bond_angles:
            job = project.open_job({'rep':rep, "bond_length":bl, "bond_angle":angle})
            sh.copytree("parameters", os.path.join(job.workspace(), "parameters"))
            job.init()
