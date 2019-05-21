"""
Coarse Grained PyRosetta
An implementation of Pyrosetta with CG capabilities
"""

# Add imports here
from .cg_pyrosetta import *

# Handle versioneer
from ._version import get_versions
versions = get_versions()
__version__ = versions['version']
__git_revision__ = versions['full-revisionid']
del get_versions, versions
