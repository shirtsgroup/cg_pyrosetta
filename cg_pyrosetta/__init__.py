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
configs_file = open(os.path.join(current_path, '../.configs.yml'), 'r')
configs = yaml.safe_load(configs_file)
clean_pyrosetta_path = configs['clean_pyrosetta_path']

builder = cg_pyrosetta.build_cg_pyrosetta.PyRosettaBuilder(clean_pyrosetta_path, pyrosetta_path, data_path)
builder.buildCGPyRosetta()

import cg_pyrosetta.CG_movers
import cg_pyrosetta.CG_folding
import cg_pyrosetta.change_parameters
import cg_pyrosetta.CG_monte_carlo
import pyrosetta

# pyrosetta.init("--add_mm_atom_type_set_parameters fa_standard mm_atom_type_sets/mm_atom_properties.txt " +
#                    "--extra_mm_params_dir mm_atom_type_sets")
# Handle versioneer
from ._version import get_versions
versions = get_versions()
__version__ = versions['version']
__git_revision__ = versions['full-revisionid']
del get_versions, versions