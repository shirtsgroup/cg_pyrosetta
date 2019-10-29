"""
Coarse Grained PyRosetta
An implementation of Pyrosetta with CG capabilities
"""

# Add imports here
from .cg_pyrosetta import *
import cg_pyrosetta.build_cg_pyrosetta
import os
import yaml

current_path = os.path.dirname(os.path.abspath(__file__))
pyrosetta_path = os.path.join(current_path, '../PyRosetta4.modified')
data_path    = os.path.join(current_path, 'data')
#configs_file = open(os.path.join(current_path.split('/cg_pyrosetta/')[1], '.configs.yml'), 'r')
#configs = yaml.load(configs_file)
clean_pyrosetta_path = "/home/gmeek/software/pyrosetta"

builder = cg_pyrosetta.build_cg_pyrosetta.PyRosettaBuilder(clean_pyrosetta_path, pyrosetta_path, data_path)
builder.buildCGPyRosetta()

import cg_pyrosetta.CG_movers
import cg_pyrosetta.CG_folding
import cg_pyrosetta.change_parameters
import pyrosetta

# Handle versioneer
from ._version import get_versions
versions = get_versions()
__version__ = versions['version']
__git_revision__ = versions['full-revisionid']
del get_versions, versions
