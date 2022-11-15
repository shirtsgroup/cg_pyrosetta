import signac
import numpy as np
import os
import shutil as sh

project = signac.get_project()

sc_sizes = np.linspace(0.5, 3, 11)
bond_angles = np.linspace(100, 160, 16)

for sc_size in sc_sizes:
    for ba in bond_angles:
        for rep in range(100):
            job = project.open_job({'sc_size':sc_size, 'bond_angle':ba, 'rep':rep})
            sh.copytree("parameters", os.path.join(job.workspace(), "parameters"))
            job.init()
