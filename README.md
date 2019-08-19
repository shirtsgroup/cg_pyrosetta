Coarse Grained PyRosetta
==============================
[//]: # (Badges)
### Description

An implementation of Pyrosetta with CG capabilities. This repository is currently underdevelopment. Current features include:

- working 1-1, 2-1, 3-1, 1-2 and 1-3 models
- working generalized backbone for these models (Shear, Small and Min)
- working MC objects
- working MC annealing folding objects.
- working torsion energies (mm_twist)

### Installation

To install a local version of this package, first change the `clean_pyrosetta_path` variable in setup.py to a working version of PyRosetta4.

`clean_pyrosetta_path = /path/to/PyRosetta4`

Then run the following from the top directory:

`pip install -e .`

### Copyright

Copyright (c) 2019, Theodore Fobe


#### Acknowledgements
 
Project based on the 
[Computational Molecular Science Python Cookiecutter](https://github.com/molssi/cookiecutter-cms) version 1.0.
