Coarse-Grained PyRosetta
==============================

[![DOI](https://zenodo.org/badge/187856626.svg)](https://zenodo.org/badge/latestdoi/187856626)

### Description

A package that implements arbitrary CG models in Pyrosetta. Poses with CG atoms can be folded using the MC objects provided in this package, or using other Rosetta sampling objects. 

### Installation

CG PyRosetta depends on [PyRosetta](https://www.pyrosetta.org), which is distributed by the Rosetta Commons under its
own license. PyRosetta is free for academic and non-commercial use, and commercial use requires a separate license.
Please review the [PyRosetta license](https://www.pyrosetta.org/home/licensing-pyrosetta) before installing.

To install a local version of this package, clone the repository and create a new `conda` environment using the
`environment.yml` file included in the root directory of this package:

`conda env create --name cg_pyrosetta --file environment.yml`

Once installed, activate the `conda` environment with:

`conda activate cg_pyrosetta`

Now we can install CG PyRosetta into this new conda environment with:

`pip install -e .`

The package needs to be installed in developer mode (the `-e` flag), because it builds a modified copy of PyRosetta
next to the source code. Developer mode also lets users make and use changes to the code without needing to reinstall
the program.

The first time `cg_pyrosetta` is imported it copies the installed PyRosetta package into a `pyrosetta_modified`
directory in the root directory of this package and adds the CG patches to the copy. This can take a few minutes and
uses about 3.5 GB of disk space, and only happens once. Deleting `pyrosetta_modified` forces a rebuild on the next import.

The pinned environment (Python 3.12 and PyRosetta 2026.34) is the one the test suite runs against.

### Getting started

```python
import cg_pyrosetta

cg_pyrosetta.init()

pose = cg_pyrosetta.pyrosetta.pose_from_sequence("X[CG11x3:CGLower]X[CG11x3]X[CG11x3]X[CG11x3:CGUpper]")

energy_function = cg_pyrosetta.CG_monte_carlo.EnergyFunctionFactory().build_energy_function(
    {"fa_atr": 1, "fa_rep": 1, "mm_twist": 1, "mm_bend": 1}
)
print(energy_function(pose))
```

`cg_pyrosetta.init()` initializes PyRosetta with the CG atom types, residue types and molecular mechanics parameters
in `cg_pyrosetta/data`, and `cg_pyrosetta.pyrosetta` is the modified copy of PyRosetta that should be used in place
of `import pyrosetta`. Complete Monte Carlo annealing workflows can be found in the `simulations` directory. Those
scripts also require [signac](https://signac.io) and signac-flow.

### Testing

The tests are run from the root directory of this package with:

`pytest -v cg_pyrosetta/tests/`

### Copyright

Copyright (c) 2019, University of Colorado Boulder


#### Acknowledgements
 
Project based on the 
[Computational Molecular Science Python Cookiecutter](https://github.com/molssi/cookiecutter-cms) version 1.0.
