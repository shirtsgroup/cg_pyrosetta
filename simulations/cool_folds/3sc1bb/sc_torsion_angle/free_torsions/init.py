import signac
import numpy as np
import os
import shutil as sh

project = signac.get_project()

sc_bond_angles = np.linspace(90, 180, 19)

for sc_ba in sc_bond_angles:
    for rep in range(100):
        job = project.open_job({'sc_bond_angle':sc_ba, 'rep':rep})
        sh.copytree("parameters", os.path.join(job.workspace(), "parameters"))
        job.init()
