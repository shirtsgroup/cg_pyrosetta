import signac
import numpy as np

project = signac.get_project()

bb_lengths = np.linspace(0.5, 5, 10)
bbb_angles = np.linspace(45, 135, 10)

for bb_l in 