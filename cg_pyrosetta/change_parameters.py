import shutil
import os
import re
import cg_pyrosetta.build_cg_pyrosetta

current_path = os.path.dirname(os.path.abspath(__file__))
data_path = os.path.join(current_path, 'data')

def changeAtomParameters(param_dict):
    """
    function to change atom parameters on the fly

    Arguments
    ---------

    dict : dict
        Dictionary containing name of atom and list of specific parameters for
        that atom type. (Parameters are expected in the following order
        ATOM LJ_RADIUS LJ_WDEPTH LK_DGFREE LK_LAMBDA LK_VOLUME)

    Examples
    --------

    params = {'CG1':['X', 1.0, 0.2, 1.0, 3.5, 23.7]}
    cg_pyrosetta.change_parameters.changeAtomParameters(params)

    """

    with open(os.path.join(data_path, 'atom_type_sets', 'atom_properties.txt'), 'r') as f:
        atom_lines = f.readlines()

    atom_params_list = [line.rstrip('\n').split() for line in atom_lines[1:]]
    prev_atoms = [atom_list[0] for atom_list in atom_params_list]

    names = param_dict.keys()

    for name in names:
        assert len(param_dict[name]) >= 3

        # extract LJ parameters and atom name
        atom_type = param_dict[name][0]
        lj_radius = param_dict[name][1]
        lj_wdepth = param_dict[name][2]

        # can take either a list of length 3 or 6

        # if LK parameters are provided, will use what is provided
        if len(param_dict[name]) == 6:
            lk_dgfree = param_dict[name][3]
            lk_lambda = param_dict[name][4]
            lk_volume = param_dict[name][5]
        # otherwise use default values
        else:
            lk_dgfree = 1.0
            lk_lambda = 3.5
            lk_volume = 23.7

        if name in prev_atoms:
            i_atom = prev_atoms.index(name)
            atom_params_list[i_atom][0] = name
            atom_params_list[i_atom][1] = atom_type
            atom_params_list[i_atom][2] = lj_radius
            atom_params_list[i_atom][3] = lj_wdepth
            atom_params_list[i_atom][4] = lk_dgfree
            atom_params_list[i_atom][5] = lk_lambda
            atom_params_list[i_atom][6] = lk_volume
        else:
            new_entry = [name,
                         atom_type,
                         lj_radius,
                         lj_wdepth,
                         lk_dgfree,
                         lk_lambda,
                         lk_volume]
            atom_params_list.append(new_entry)

    # Rewrite parameters to atom_properties.txt file
    with open(os.path.join(data_path, 'atom_type_sets', 'atom_properties.txt'), 'w') as f:
        # write header
        f.write('NAME  ATOM LJ_RADIUS LJ_WDEPTH LK_DGFREE LK_LAMBDA LK_VOLUME\n')

        # rewrite updated parameters
        for param in atom_params_list:
            # print param line with specific formating
            f.write("%s%6.1s%10.4f%10.4f%10.4f%10.4f%10.4f\n" % (param[0], param[1], float(
                param[2]), float(param[3]), float(param[4]), float(param[5]), float(param[6])))

    with open(os.path.join(data_path, 'mm_atom_type_sets', 'mm_atom_properties.txt'), 'w') as f:
        # write header
        f.write('NAME    LJ_WDEPTH   LJ_RADIUS   LJ_3B_WDEPTH    LJ_3B_RADIUS\n')

        # rewrite updated parameters
        for param in atom_params_list:
            # print param line with specific formating
            f.write("%s%10.4f%10.4f%10.4f%10.4f\n" %
                    (param[0], -float(param[3]), float(param[2]), -float(param[3]), float(param[2])))
    cg_pyrosetta.builder.buildCGPyRosetta()

def changeTorsionParameters(param_dict):
    """
    function to change torsion parameters on the fly

    Arguments
    ---------

    dict : dict
        Dictionary containing name of torsion and list of specific parameters for
        that atom type. (Parameters are expected in the following order
        ATOM LJ_RADIUS LJ_WDEPTH LK_DGFREE LK_LAMBDA LK_VOLUME)

    Examples
    --------

    params = {'CG1 CG1 CG1 CG1':[k_force, periodicity, phase_shift]}
    cg_pyrosetta.change_parameters.changeAtomParameters(params)

    """

    with open(os.path.join(data_path, 'mm_atom_type_sets', 'mm_torsion_params.txt'), 'r') as f:
        torsion_lines = f.readlines()

    # Build torsion params list [name_of_torsion k_constant periodicity phase_shift]
    torsion_params_list = []
    for line in torsion_lines[1:]:
        split_line = line.rstrip('\n').split()
        name = ' '.join(split_line[:4])
        torsion_params_list.append([name, split_line[4], split_line[5], split_line[6]])

    prev_torsion = [torsion[0] for torsion in torsion_params_list]
    names = param_dict.keys()

    for name in names:

        # extract LJ parameters and atom name
        force_constant = param_dict[name][0]
        periodicity = param_dict[name][1]
        phase_shift = param_dict[name][2]

        if name in prev_torsion:
            i_torsion = prev_torsion.index(name)
            torsion_params_list[i_torsion][0] = name
            torsion_params_list[i_torsion][1] = force_constant
            torsion_params_list[i_torsion][2] = periodicity
            torsion_params_list[i_torsion][3] = phase_shift

        else:
            new_entry = [name,
                         force_constant,
                         periodicity,
                         phase_shift]
            torsion_params_list.append(new_entry)

    # Rewrite parameters to atom_properties.txt file
    with open(os.path.join(data_path, 'mm_atom_type_sets', 'mm_torsion_params.txt'), 'w') as f:
        # write header
        f.write('# CG torsion parameters\n')

        # rewrite updated parameters
        for param in torsion_params_list:
            # print param line with specific formating
            atom_names = param[0].split()
            f.write("%s %s %s %s %.4f %.1i %.4f\n" % (
                atom_names[0], atom_names[1], atom_names[2], atom_names[3], float(param[1]), int(param[2]), float(param[3])))
    cg_pyrosetta.builder.buildCGPyRosetta()


def changeAngleParameters(param_dict):
    """
    function to change angle parameters (mm_bend)on the fly

    Arguments
    ---------

    dict : dict
        Dictionary containing name of atom and list of specific parameters for
        that atom type. (Parameters are expected in the following order
        !atom types     Ktheta    Theta0   Kub     S0)

    Examples
    --------

    params = {'CG1':['X', 1.0, 0.2, 1.0, 3.5, 23.7]}
    cg_pyrosetta.change_parameters.changeAtomParameters(params)

    """

    with open(os.path.join(data_path, 'mm_atom_type_sets', 'mm_angle_params.txt'), 'r') as f:
        angle_lines = f.readlines()

    # Build torsion params list [name_of_torsion k_constant periodicity phase_shift]
    angle_params_list = []
    for line in angle_lines[2:]:
        split_line = line.rstrip('\n').split()
        name = ' '.join(split_line[:3])
        angle_params_list.append([name, split_line[3], split_line[4]])

    prev_angles = [angle[0] for angle in angle_params_list]
    names = param_dict.keys()
    print(prev_angles)
    for name in names:

        # extract LJ parameters and atom name
        force_constant = param_dict[name][0]
        angle = param_dict[name][1]

        if name in prev_angles:
            i_angle = prev_angles.index(name)
            angle_params_list[i_angle][0] = name
            angle_params_list[i_angle][1] = force_constant
            angle_params_list[i_angle][2] = angle

        else:
            new_entry = [name,
                         force_constant,
                         angle, ]
            angle_params_list.append(new_entry)

    # Rewrite parameters to atom_properties.txt file
    with open(os.path.join(data_path, 'mm_atom_type_sets', 'mm_angle_params.txt'), 'w') as f:
        # write header
        f.write('ANGLES\n')
        f.write('!atom types     Ktheta    Theta0   Kub     S0\n')

        # rewrite updated parameters
        for param in angle_params_list:
            # print param line with specific formating
            atom_names = param[0].split()
<<<<<<< HEAD
            f.write("%-4.4s%4.3s%4.3s%10.4f%10.4f\n" % (atom_names[0], atom_names[1], atom_names[2], float(param[1]), float(param[2])))
=======
            f.write("%-4.4s%4.3s%4.3s%10.4f%10.4f\n" %
                    (atom_names[0], atom_names[1], atom_names[2], float(param[1]), float(param[2])))
    cg_pyrosetta.builder.buildCGPyRosetta()
>>>>>>> 71112bd69634557d3d37866794a0f00e6ffdf80a
