"""
Unit and regression test for the cg_pyrosetta package.
"""

import subprocess
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


def test_import_from_another_directory(tmp_path):
    """
    The modified copy of PyRosetta must be importable from any working directory
    """
    result = subprocess.run(
        [sys.executable, "-c", "import cg_pyrosetta; print(cg_pyrosetta.pyrosetta.__name__)"],
        cwd=str(tmp_path), capture_output=True, text=True,
    )
    assert result.returncode == 0, result.stderr
    assert result.stdout.strip().endswith("pyrosetta_modified")
