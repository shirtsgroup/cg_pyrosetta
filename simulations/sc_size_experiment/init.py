import signac
import numpy as np
import os
import shutil as sh

project = signac.get_project()

sc_sizes = np.linspace(0.5, 5, 50)

for sc_size in sc_sizes:
        job = project.open_job({'sc_size':sc_size})
        sh.copytree("parameters", os.path.join(job.workspace(), "parameters"))
        job.init()
