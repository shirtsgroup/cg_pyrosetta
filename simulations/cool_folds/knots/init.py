import signac
import numpy as np
import os
import shutil as sh

project = signac.get_project()
bond_angles = [145, 147.5, 150, 152.5, 155, 157.5, 160]

for rep in range(200):
    for angle in bond_angles:
        job = project.open_job({'rep':rep, "bond_angle":angle})
        sh.copytree("parameters", os.path.join(job.workspace(), "parameters"))
        job.init()
