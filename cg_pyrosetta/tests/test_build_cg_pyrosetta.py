"""
Tests confirming the modified copy of PyRosetta was built correctly.
"""

import os

import cg_pyrosetta

FA_STANDARD_ATOMS = ('database', 'chemical', 'atom_type_sets', 'fa_standard')
FA_STANDARD_RESIDUES = ('database', 'chemical', 'residue_type_sets', 'fa_standard')


def test_modified_pyrosetta_exists():
    """
    Confirms the modified copy of PyRosetta is present and has a database
    """
    pyrosetta_path = os.path.abspath(cg_pyrosetta.pyrosetta_path)
    assert os.path.isfile(os.path.join(pyrosetta_path, '__init__.py'))
    assert os.path.isdir(os.path.join(pyrosetta_path, 'database'))


def test_modified_pyrosetta_is_used(pyrosetta):
    """
    Confirms cg_pyrosetta.pyrosetta is the modified copy, not the installed one
    """
    module_dir = os.path.dirname(os.path.abspath(pyrosetta.__file__))
    assert os.path.samefile(module_dir, cg_pyrosetta.pyrosetta_path)


def test_extras_are_turned_off():
    """
    Confirms every line of extras.txt in the atom type set is commented out
    """
    with open(os.path.join(cg_pyrosetta.pyrosetta_path, *FA_STANDARD_ATOMS, 'extras.txt')) as exf:
        lines = exf.readlines()

    assert all(line[0] == '#' for line in lines)


def test_patches_are_added_once():
    """
    Confirms the custom patches are listed in patches.txt exactly once
    """
    with open(os.path.join(cg_pyrosetta.pyrosetta_path, *FA_STANDARD_RESIDUES, 'patches.txt')) as patch_file:
        lines = patch_file.readlines()

    assert lines.count('### custom patches residues\n') == 1
    for patch in ('CG_lower_patch.txt', 'CG_upper_patch.txt'):
        assert sum(line.rstrip('\n').endswith(patch) for line in lines) == 1


def test_pose_from_sequence(pyrosetta):
    """
    Confirms a 1-1 CG model is built with the expected atoms
    """
    pose = pyrosetta.pose_from_sequence('X[CG11]' * 8)
    assert pose.size() == 8
    assert pose.residue(1).atom_name(1) == 'BB1 '
    assert pose.residue(1).atom_name(2) == 'SC1 '


def test_terminal_patches_apply(pose):
    """
    Confirms the custom terminus patches can be applied to a CG residue
    """
    assert pose.residue(1).name() == 'CG11x3:CGLower'
    assert pose.residue(pose.size()).name() == 'CG11x3:CGUpper'
