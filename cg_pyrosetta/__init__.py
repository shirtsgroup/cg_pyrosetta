"""
Coarse Grained PyRosetta
An implementation of Pyrosetta with CG capabilities
"""

# Add imports here
from .cg_pyrosetta import *
from .build_cg_pyrosetta import *

# pyrosetta_path = os.path.abspath('PyRosetta4.modified')
# input_path = os.path.abspath('cg_pyrosetta/data')
    
# builder = PyRosettaBuilder(pyrosetta_path, input_path)
# builder.buildCGPyRosetta()

from .CG_movers import *
from .CG_folding import *


# Handle versioneer
from ._version import get_versions
versions = get_versions()
__version__ = versions['version']
__git_revision__ = versions['full-revisionid']
del get_versions, versions
