import shutil
import os
import re

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
                         lk_volume ]
            atom_params_list.append(new_entry)
            
    # Rewrite parameters to atom_properties.txt file
    with open(os.path.join(data_path, 'atom_type_sets', 'atom_properties.txt'), 'w') as f:
        # write header
        f.write('NAME  ATOM LJ_RADIUS LJ_WDEPTH LK_DGFREE LK_LAMBDA LK_VOLUME\n')
        
        # rewrite updated parameters
        for param in atom_params_list:
            # print param line with specific formating
            f.write("%s%6.1s%10.4f%10.4f%10.4f%10.4f%10.4f\n" % (param[0], param[1], float(param[2]), float(param[3]), float(param[4]), float(param[5]), float(param[6])))
        
    
        
        
    