import signac
import numpy as np
import os
import shutil as sh

project = signac.get_project()

bb_lengths = np.linspace(0.5, 5, 5)
bbb_angles = np.linspace(45, 135, 5)

for bb_l in bb_lengths:
    for bbb_angle in bbb_angles:
        job = project.open_job({'bb_length' : bb_l, 'bbb_angle' : bbb_angle})
        sh.copytree("parameters", os.path.join(job.workspace(), "parameters"))
        job.init()
