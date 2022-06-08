import signac
import numpy as np
import os
import shutil as sh

project = signac.get_project()

for rep in range(100):
    job = project.open_job({'rep':rep})
    sh.copytree("parameters", os.path.join(job.workspace(), "parameters"))
    job.init()