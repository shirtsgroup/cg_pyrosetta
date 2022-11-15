import signac
import numpy as np
import os
import shutil as sh

project = signac.get_project()
hinge_angles = np.linspace(60, 160, 20)

for rep in range(100):
    for angle in hinge_angles:
        job = project.open_job({'rep':rep, "hinge_angle":angle})
        sh.copytree("parameters", os.path.join(job.workspace(), "parameters"))
        job.init()
