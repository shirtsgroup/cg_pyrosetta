import signac
import numpy as np
import os
import shutil as sh

project = signac.get_project()

nmer_lengths = np.array([3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18])
reps = 20

for nmer_length in nmer_lengths:
    for rep in range(reps):
        job = project.open_job({'nmer' : nmer_length, 'rep' : rep})
        sh.copytree("parameters", os.path.join(job.workspace(), "parameters"))
        job.init()
