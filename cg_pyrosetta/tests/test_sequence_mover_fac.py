import pytest
from cg_pyrosetta.CG_monte_carlo import SequenceMoverFactory
import os
import sys
import warnings

current_path = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.abspath(current_path + '/../PyRosetta4.modified'))

import pyrosetta
import numpy as np

@pytest.fixture
def pose():
    return(pyrosetta.pose_from_sequence('X[CG11]X[CG11]X[CG11]X[CG11]X[CG11]'))

def test_sequence_factory(pose):
    mover_list = ['small_dihe', 'shear_dihe']
    freq_list = [10, 10]
    mover_seq_builder = SequenceMoverFactory()

    mover_seq = mover_seq_builder.build_seq_mover(pose, mover_list, freq_list)
    assert(mover_seq.size() == 2)

def test_sequence_factory_warn(pose):
    mover_list = ['unimplemented_mover', 'small_dihe', 'shear_dihe']
    freq_list = [10, 10, 10]
    mover_seq_builder = SequenceMoverFactory()
    with pytest.warns(UserWarning):
        mover_seq_builder.build_seq_mover(pose, mover_list, freq_list)

