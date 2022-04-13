import signac
import numpy as np
import os
import shutil as sh

project = signac.get_project()

sc_sizes = np.linspace(0.5, 4, 96)

for sc_size in sc_sizes:
    for rep in range(100):
        job = project.open_job({'sc_size':sc_size, 'rep':rep})
        sh.copytree("parameters", os.path.join(job.workspace(), "parameters"))
        job.init()
