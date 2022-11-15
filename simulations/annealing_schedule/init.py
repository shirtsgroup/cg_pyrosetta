import signac
import numpy as np
import os
import shutil as sh

project = signac.get_project()

kt_initials = np.logspace(1, 3, base=10, num=10)
geometric_rates = np.linspace(0.5, 0.999, 10)
replicas = 50

for kt_i in kt_initials:
    for r in geometric_rates:
        for i in range(replicas):
            job = project.open_job({'kt_initial': kt_i, 'geometric_rate' : r, 'rep':i})
            sh.copytree("parameters", os.path.join(job.workspace(), "parameters"))
            job.init()
