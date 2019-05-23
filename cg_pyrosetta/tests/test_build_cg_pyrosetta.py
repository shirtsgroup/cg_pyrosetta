"""
Unit and regression test for the cg_pyrosetta package.
"""

# File containing tests to confirm the build_cg_pyrosetta script is working
import pytest
import sys
import os
import cg_pyrosetta


def testPyRosettaModified():
    """
    Confirms presences of working PyRosetta buld in PyRosetta4.modified
    """
    pyrosetta_path = os.path.abspath(cg_pyrosetta.pyrosetta_path)
    assert 'pyrosetta' in os.listdir(pyrosetta_path)


def testAtoms():
    """
    Checks to see if header was written into the atom_properties.txt file
    """
    with open(os.path.join(cg_pyrosetta.pyrosetta_path, 'pyrosetta','database','chemical','atom_type_sets','fa_standard','atom_properties.txt')) as atom_file:
        lines = atom_file.readlines()
    assert '## Custom Added Atom Types ###\n' in lines

def testResidueDir():
    assert os.path.isdir(os.path.join(cg_pyrosetta.pyrosetta_path,'pyrosetta','database','chemical','residue_type_sets','fa_standard','residue_types','custom'))

def testResiduefiles():
    """
    Confirms the 'custom directory' is populated
    """

    assert os.listdir(os.path.join(cg_pyrosetta.pyrosetta_path,'pyrosetta','database','chemical','residue_type_sets','fa_standard','residue_types','custom'))

def testExtras():
    """
    Confirms the comments in extras.txt file in atom_type_set dir
    """
    with open (os.path.join(cg_pyrosetta.pyrosetta_path,'pyrosetta','database','chemical','atom_type_sets','fa_standard','extras.txt')) as exf:
        lines = exf.readlines()
    
    assert all([line[0] == '#' for line in lines])

def testResidueTypetxt():
    """
    confirms atom residue types are written to residue_types.txt by checking the header
    """
    with open(os.path.join(cg_pyrosetta.pyrosetta_path, 'pyrosetta','database','chemical','residue_type_sets','fa_standard','residue_types.txt'), 'r') as rtf:
        lines = rtf.readlines()
    assert '### custom residues\n' in lines


def testPoseSequence():
    """
    confirms 1-1 CG model has correct atom types in them
    """
    seq = cg_pyrosetta.CG_folding.pyrosetta.pose_from_sequence('X[CG11]X[CG11]X[CG11]X[CG11]X[CG11]X[CG11]X[CG11]X[CG11]')
    assert seq.residue(1).atom_name(1) == 'BB1 '


