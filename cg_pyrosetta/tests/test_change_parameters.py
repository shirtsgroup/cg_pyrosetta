import cg_pyrosetta
import pyrosetta
import os
import pytest
import numpy as np

current_path = os.path.dirname(os.path.abspath(__file__))
data_path = os.path.join(current_path, '../data')
pyrosetta_top_path = os.path.join(data_path, '../../PyRosetta4.modified')

def test_change_atom_parameters():
    cg_pyrosetta.change_parameters.changeAtomParameters(
                        {'CG1':['X', 10.0, 10.2, 1.0, 1.5, 23.7]})

    # check local data dir (mm + fa_standard)
    with open(os.path.join(data_path, "atom_type_sets", "atom_properties.txt"), 'r') as fhandle:
        fa_lines = fhandle.readlines()
    with open(os.path.join(data_path, "mm_atom_type_sets", "mm_atom_properties.txt"), 'r') as fhandle:
        mm_lines = fhandle.readlines()
    

    # check pyrosetta data 
    
    atom_prop_line = "%s%6.1s%10.4f%10.4f%10.4f%10.4f%10.4f\n" % ('CG1', 'X', 10.0, 10.2, 1.0, 1.5, 23.7)
    mm_atom_prop_line = "%s%10.4f%10.4f%10.4f%10.4f\n" % ('CG1', -10.2, 10, -10.2, 10)
    
    # confirms changes in data dir
    assert atom_prop_line in fa_lines
    assert mm_atom_prop_line in mm_lines

    # confirms changes in PyRosetta4.modified

    with open(os.path.join(pyrosetta_top_path, "pyrosetta", "database", "chemical", "atom_type_sets", "fa_standard", "atom_properties.txt"), 'r') as fhandle:
        pyrosetta_fa_lines = fhandle.readlines()
    with open(os.path.join(pyrosetta_top_path, "pyrosetta", "database", "chemical", "mm_atom_type_sets", "fa_standard", "mm_atom_properties.txt"), 'r') as fhandle:
        pyrosetta_mm_lines = fhandle.readlines()

    assert atom_prop_line in pyrosetta_fa_lines
    assert mm_atom_prop_line in pyrosetta_mm_lines

    pose = cg_pyrosetta.pyrosetta.pose_from_sequence('X[CG31:CGUpper]X[CG31:CGLower]')
    sf = cg_pyrosetta.pyrosetta.ScoreFunction()
    sf.set_weight(cg_pyrosetta.pyrosetta.rosetta.core.scoring.mm_lj_inter_atr, 1)
    sf.set_weight(cg_pyrosetta.pyrosetta.rosetta.core.scoring.mm_lj_inter_rep, 1)
    sf.set_weight(cg_pyrosetta.pyrosetta.rosetta.core.scoring.mm_lj_intra_atr, 1)
    sf.set_weight(cg_pyrosetta.pyrosetta.rosetta.core.scoring.mm_lj_intra_rep, 1)

    old_score = sf(pose)
    print("Old Score:", old_score)

    cg_pyrosetta.change_parameters.changeAtomParameters(
                    {'CG1':['X', 1.0, 0.2, 1.0, 1.5, 23.7]})

    cg_pyrosetta.pyrosetta.init()
    new_sf = cg_pyrosetta.pyrosetta.ScoreFunction()
    new_sf.set_weight(cg_pyrosetta.pyrosetta.rosetta.core.scoring.mm_lj_inter_atr, 1)
    new_sf.set_weight(cg_pyrosetta.pyrosetta.rosetta.core.scoring.mm_lj_inter_rep, 1)
    new_sf.set_weight(cg_pyrosetta.pyrosetta.rosetta.core.scoring.mm_lj_intra_atr, 1)
    new_sf.set_weight(cg_pyrosetta.pyrosetta.rosetta.core.scoring.mm_lj_intra_rep, 1)

    new_score = new_sf(pose)
    print("New Score:", new_score)

    assert new_score !=  old_score