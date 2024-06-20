"""
Coarse Grained PyRosetta
An implementation of Pyrosetta with CG capabilities
"""

# Add imports here
from .cg_pyrosetta import *
from cg_pyrosetta.build_cg_pyrosetta import PyRosettaBuilder, add_import_path
import os
import yaml
import pyrosetta as old_pyrosetta

current_path = os.path.dirname(os.path.abspath(__file__))

import cg_pyrosetta.CG_movers
import cg_pyrosetta.CG_folding
import cg_pyrosetta.parameters
import cg_pyrosetta.CG_monte_carlo
import cg_pyrosetta.scripts
from cg_pyrosetta.utils import init
import pyrosetta_modified as pyrosetta

extra_residue_files = [os.path.abspath(os.path.join(current_path, "data", "residue_type_sets", a)) \
    for a in os.listdir(os.path.join(current_path, "data", "residue_type_sets"))]

for filename in extra_residue_files:
    if ".params" not in filename:
        extra_residue_files.remove(filename)


cmd_line_options_defaults = {"--add_atom_types" : "fa_standard " + os.path.join(current_path, "data", "atom_type_sets", "atom_properties.txt"),
"--add_mm_atom_type_set_parameters" : "fa_standard " + os.path.join(current_path, "data", "mm_atom_type_sets", "mm_atom_properties.txt"),
"--extra_mm_params_dir" : "fa_standard " + os.path.join(current_path, "data", "mm_atom_type_sets"),
"--mute" : "all",
"--extra_res_fa" : extra_residue_files,
}

# Handle versioneer
from ._version import get_versions
versions = get_versions()
__version__ = versions['version']
__git_revision__ = versions['full-revisionid']
del get_versions, versions