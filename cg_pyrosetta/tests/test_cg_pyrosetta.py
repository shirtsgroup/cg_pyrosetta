"""
Unit and regression test for the cg_pyrosetta package.
"""

import sys

import cg_pyrosetta


def test_cg_pyrosetta_imported():
    """Sample test, will always pass so long as import statement worked"""
    assert "cg_pyrosetta" in sys.modules


def test_version_is_set():
    assert isinstance(cg_pyrosetta.__version__, str)
    assert cg_pyrosetta.__version__


def test_extra_residue_files_are_params_files():
    assert cg_pyrosetta.extra_residue_files
    assert all(filename.endswith(".params") for filename in cg_pyrosetta.extra_residue_files)
