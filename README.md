Coarse-Grained PyRosetta
==============================

[![DOI](https://zenodo.org/badge/187856626.svg)](https://zenodo.org/badge/latestdoi/187856626)

### Description

A package that implements arbitrary CG models in Pyrosetta. Poses with CG atoms can be folded using the MC objects provided in this package, or using other Rosetta sampling objects. 

### Installation

To install a local version of this package, create a new `conda` environment using the `environment.yml` file included in the root directory of this package:

`conda env create --name cg_pyrosetta --file environment.yml`

Once installed, activate the `conda` environment with:
`conda activate cg_pyrosetta`

Now we can install CG Pyrosetta into this new conda environment with:

`pip install -e .`

The `-e` flag is optional, however, it compiles the package in developer mode, which lets users make and use changes to the code without needing to reinstall the program.

### Copyright

Copyright (c) 2019, University of Colorado Boulder


#### Acknowledgements
 
Project based on the 
[Computational Molecular Science Python Cookiecutter](https://github.com/molssi/cookiecutter-cms) version 1.0.
