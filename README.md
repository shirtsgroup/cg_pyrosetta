Coarse Grained PyRosetta
==============================

<!-- [//]: # (Badges)
[![Travis Build Status](https://travis-ci.org/REPLACE_WITH_OWNER_ACCOUNT/Coarse Grained PyRosetta.png)](https://travis-ci.org/REPLACE_WITH_OWNER_ACCOUNT/Coarse Grained PyRosetta)
[![AppVeyor Build status](https://ci.appveyor.com/api/projects/status/REPLACE_WITH_APPVEYOR_LINK/branch/master?svg=true)](https://ci.appveyor.com/project/REPLACE_WITH_OWNER_ACCOUNT/Coarse Grained PyRosetta/branch/master)
[![codecov](https://codecov.io/gh/REPLACE_WITH_OWNER_ACCOUNT/Coarse Grained PyRosetta/branch/master/graph/badge.svg)](https://codecov.io/gh/REPLACE_WITH_OWNER_ACCOUNT/Coarse Grained PyRosetta/branch/master)
[![DOI](https://zenodo.org/badge/187856626.svg)](https://zenodo.org/badge/latestdoi/187856626)

-->
### Description

An implementation of Pyrosetta with CG capabilities. CG PyRosetta can fold arbitrary CG models using a variety of MC minimization methods. In this repository we use the Python wrapper PyRosetta to build a framework to run CG folding simulations. Currently we have:
- working 1-1, 2-1, 3-1, 1-2 and 1-3 models
- working generalized movers for bond angles and torsions
- working MC objects
- working MC annealing folding objects
- working torsion energies (mm_twist)

### Installation

To install a local version of this package, first change the `clean_pyrosetta_path` variable in setup.py to a working version of PyRosetta4.

`clean_pyrosetta_path = /path/to/PyRosetta4`

Then run the following from the top directory:

`pip install -e .`

### Copyright

Copyright (c) 2019, University of Colorado Boulder


#### Acknowledgements
 
Project based on the 
[Computational Molecular Science Python Cookiecutter](https://github.com/molssi/cookiecutter-cms) version 1.0.
