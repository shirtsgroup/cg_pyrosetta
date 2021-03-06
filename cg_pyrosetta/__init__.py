"""
Coarse Grained PyRosetta
An implementation of Pyrosetta with CG capabilities
"""

# Add imports here
from .cg_pyrosetta import *
from cg_pyrosetta.build_cg_pyrosetta import PyRosettaBuilder, add_import_path
import os
import yaml

current_path = os.path.dirname(os.path.abspath(__file__))
pyrosetta_path = os.path.join(current_path, '../PyRosetta4.modified')
data_path    = os.path.join(current_path, 'data')
configs_file = open(os.path.join(current_path, '../.configs.yml'), 'r')
configs = yaml.safe_load(configs_file)
clean_pyrosetta_path = configs['clean_pyrosetta_path']

if not os.path.isdir(pyrosetta_path):
    builder = PyRosettaBuilder(clean_pyrosetta_path, pyrosetta_path, data_path)
    builder.buildCGPyRosetta()
else:
    add_import_path(pyrosetta_path)

import cg_pyrosetta.CG_movers
import cg_pyrosetta.CG_folding
import cg_pyrosetta.change_parameters
import cg_pyrosetta.CG_monte_carlo
from cg_pyrosetta.utils import init
import pyrosetta

# pyrosetta.init("--add_mm_atom_type_set_parameters fa_standard mm_atom_type_sets/mm_atom_properties.txt " +
# "--extra_mm_params_dir mm_atom_type_sets")

extra_residue_files = [os.path.abspath(os.path.join(current_path, "data", "residue_type_sets", a)) \
    for a in os.listdir(os.path.join(current_path, "data", "residue_type_sets"))]

for filename in extra_residue_files:
    if ".params" not in filename:
        extra_residue_files.remove(filename)


cmd_line_options_defaults = {"--add_atom_types" : "fa_standard parameters/atom_properties.txt",
"--add_mm_atom_type_set_parameters" : "fa_standard parameters/mm_atom_type_sets/mm_atom_properties.txt ",
"--extra_mm_params_dir" : "parameters/mm_atom_type_sets",
"--mute" : "all",
"--extra_res_fa" : extra_residue_files,
}

# Handle versioneer
from ._version import get_versions
versions = get_versions()
__version__ = versions['version']
__git_revision__ = versions['full-revisionid']
del get_versions, versions