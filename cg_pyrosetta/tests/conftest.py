"""
Shared fixtures for the cg_pyrosetta test suite.

PyRosetta can only be initialized once per process, so initialization is done
here in a single session-scoped fixture, with the default options and the
parameter files packaged in cg_pyrosetta/data.
"""

import os
import shutil

import pytest

import cg_pyrosetta

CG11X3_TETRAMER = 'X[CG11x3:CGLower]X[CG11x3]X[CG11x3]X[CG11x3:CGUpper]'
CG11X3_OCTAMER = 'X[CG11x3:CGLower]' + 'X[CG11x3]' * 6 + 'X[CG11x3:CGUpper]'


@pytest.fixture(scope="session", autouse=True)
def initialized_pyrosetta():
    cg_pyrosetta.init()
    return cg_pyrosetta.pyrosetta


@pytest.fixture
def pyrosetta(initialized_pyrosetta):
    return initialized_pyrosetta


@pytest.fixture
def pose(pyrosetta):
    return pyrosetta.pose_from_sequence(CG11X3_TETRAMER)


@pytest.fixture
def long_pose(pyrosetta):
    return pyrosetta.pose_from_sequence(CG11X3_OCTAMER)


@pytest.fixture
def parameter_files(tmp_path):
    """
    Copy of the packaged parameter files, so tests that rewrite parameters
    never modify the files that ship with the package.
    """
    data_path = os.path.join(os.path.dirname(os.path.abspath(cg_pyrosetta.__file__)), 'data')
    shutil.copytree(os.path.join(data_path, 'atom_type_sets'), tmp_path / 'atom_type_sets')
    shutil.copytree(os.path.join(data_path, 'mm_atom_type_sets'), tmp_path / 'mm_atom_type_sets')
    return {
        'atom_types': str(tmp_path / 'atom_type_sets' / 'atom_properties.txt'),
        'mm_atom_types': str(tmp_path / 'mm_atom_type_sets' / 'mm_atom_properties.txt'),
        'angles': str(tmp_path / 'mm_atom_type_sets' / 'mm_bond_angle_params.txt'),
        'torsions': str(tmp_path / 'mm_atom_type_sets' / 'mm_torsion_params.txt'),
    }
