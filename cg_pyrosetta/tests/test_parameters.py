"""
Tests for the functions in cg_pyrosetta.parameters that rewrite parameter files.
"""

import numpy as np
import pytest

import cg_pyrosetta
from cg_pyrosetta import parameters


def read_lines(filename):
    with open(filename, 'r') as fhandle:
        return fhandle.readlines()


def test_change_atom_parameters(parameter_files):
    parameters.changeAtomParameters(
        {'CG1': ['X', 1.5, 0.3, 1.0, 1.5, 23.7]},
        atom_types_path=parameter_files['atom_types'],
        mm_atom_types_path=parameter_files['mm_atom_types'],
    )

    atom_prop_line = "%s%6.1s%10.4f%10.4f%10.4f%10.4f%10.4f\n" % ('CG1', 'X', 1.5, 0.3, 1.0, 1.5, 23.7)
    mm_atom_prop_line = "%s%10.4f%10.4f%10.4f%10.4f\n" % ('CG1', -0.3, 1.5, -0.3, 1.5)

    assert atom_prop_line in read_lines(parameter_files['atom_types'])
    assert mm_atom_prop_line in read_lines(parameter_files['mm_atom_types'])


def test_change_atom_parameters_adds_new_atom_type(parameter_files):
    n_lines = len(read_lines(parameter_files['atom_types']))
    parameters.changeAtomParameters(
        {'CGX': ['X', 1.5, 0.3]},
        atom_types_path=parameter_files['atom_types'],
        mm_atom_types_path=parameter_files['mm_atom_types'],
    )
    lines = read_lines(parameter_files['atom_types'])
    assert len(lines) == n_lines + 1
    assert lines[-1].split() == ['CGX', 'X', '1.5000', '0.3000', '1.0000', '3.5000', '23.7000']


@pytest.mark.parametrize("convention, value, expected_radius", [
    ('rmin_half', 1.5, 1.5),
    ('rmin', 3.0, 1.5),
    ('sigma', 3.0 / 2 ** (1 / 6), 1.5),
])
def test_lj_conventions_write_the_same_radius(parameter_files, convention, value, expected_radius):
    assert parameters.lj_size_to_rosetta_radius(value, convention) == pytest.approx(expected_radius)

    parameters.changeAtomParameters(
        {'CG1': ['X', value, 0.2]},
        atom_types_path=parameter_files['atom_types'],
        mm_atom_types_path=parameter_files['mm_atom_types'],
        lj_convention=convention,
    )
    cg1_line = [line.split() for line in read_lines(parameter_files['atom_types']) if line.startswith('CG1 ')][0]
    assert float(cg1_line[2]) == pytest.approx(expected_radius, abs=1e-4)


def test_invalid_lj_convention(parameter_files):
    with pytest.raises(ValueError):
        parameters.lj_size_to_rosetta_radius(1.0, 'diameter')

    before = read_lines(parameter_files['atom_types'])
    with pytest.raises(ValueError):
        parameters.changeAtomParameters(
            {'CG1': ['X', 1.0, 0.2]},
            atom_types_path=parameter_files['atom_types'],
            mm_atom_types_path=parameter_files['mm_atom_types'],
            lj_convention='diameter',
        )
    assert read_lines(parameter_files['atom_types']) == before


def test_lj_radius_is_half_of_rmin(pyrosetta):
    """
    Confirms the convention behind lj_size_to_rosetta_radius. Two CG1 beads
    with the packaged LJ_RADIUS and LJ_WDEPTH have their energy minimum at a
    separation of 2 * LJ_RADIUS with a depth of LJ_WDEPTH, and the energy
    crosses zero at sigma = 2 * LJ_RADIUS / 2**(1/6).
    """
    data_file = cg_pyrosetta.cmd_line_options_defaults["--add_atom_types"].split()[-1]
    cg1_line = [line.split() for line in read_lines(data_file) if line.startswith('CG1 ')][0]
    lj_radius, lj_wdepth = float(cg1_line[2]), float(cg1_line[3])

    pose = pyrosetta.pose_from_sequence('X[CG11]')
    pose.append_pose_by_jump(pyrosetta.pose_from_sequence('X[CG11]'), 1)

    # Keep the sidechain and virtual atoms out of range of the two backbone
    # beads, but close enough that the residues stay neighbors.
    vector = pyrosetta.rosetta.numeric.xyzVector_double_t
    pose.set_xyz(pyrosetta.AtomID(1, 1), vector(0, 0, 0))
    pose.set_xyz(pyrosetta.AtomID(2, 1), vector(0, 0, -6.5))
    pose.set_xyz(pyrosetta.AtomID(3, 1), vector(0, 0, -7.5))
    pose.set_xyz(pyrosetta.AtomID(2, 2), vector(0, 0, 6.5))
    pose.set_xyz(pyrosetta.AtomID(3, 2), vector(0, 0, 7.5))

    score_function = pyrosetta.ScoreFunction()
    score_function.set_weight(pyrosetta.rosetta.core.scoring.fa_atr, 1)
    score_function.set_weight(pyrosetta.rosetta.core.scoring.fa_rep, 1)

    distances = np.arange(0.6 * 2 * lj_radius, 2.0 * 2 * lj_radius, 0.01)
    energies = []
    for distance in distances:
        pose.set_xyz(pyrosetta.AtomID(1, 2), vector(distance, 0, 0))
        energies.append(score_function(pose))
    energies = np.array(energies)

    assert distances[energies.argmin()] == pytest.approx(2 * lj_radius, abs=0.02)
    assert energies.min() == pytest.approx(-lj_wdepth, abs=1e-3)
    first_negative = distances[np.argmax(energies < 0)]
    assert first_negative == pytest.approx(2 * lj_radius / 2 ** (1 / 6), abs=0.02)


def test_change_angle_parameters(parameter_files):
    parameters.changeAngleParameters(
        {'CG2 CG1 CG1': [15, 110], 'CG9 CG9 CG9': [5, 90]},
        angle_file=parameter_files['angles'],
    )
    entries = {tuple(line.split()[:3]): [float(x) for x in line.split()[3:5]]
               for line in read_lines(parameter_files['angles'])[2:] if line.strip()}
    assert entries[('CG2', 'CG1', 'CG1')] == [15.0, 110.0]
    assert entries[('CG9', 'CG9', 'CG9')] == [5.0, 90.0]


def test_change_angle_parameters_default_file_exists():
    """
    The default angle file name must match the file packaged with cg_pyrosetta
    """
    import inspect
    import os

    source = inspect.getsource(parameters.changeAngleParameters)
    assert 'mm_bond_angle_params.txt' in source
    assert os.path.isfile(os.path.join(parameters.data_path, 'mm_atom_type_sets', 'mm_bond_angle_params.txt'))


def test_change_torsion_parameters(parameter_files):
    parameters.changeTorsionParameters(
        {'CG2 CG1 CG1 CG2': [3, 2, 45]},
        torsion_file=parameter_files['torsions'],
    )
    entries = {tuple(line.split()[:4]): [float(x) for x in line.split()[4:7]]
               for line in read_lines(parameter_files['torsions'])[1:] if line.strip()}
    assert entries[('CG2', 'CG1', 'CG1', 'CG2')] == [3.0, 2.0, 45.0]
