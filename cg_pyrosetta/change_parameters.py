import shutil
import os
import re

current_path = os.path.dirname(os.path.abspath(__file__))
data_path = os.path.join(current_path, 'data')

def read_mm_atom_properties_txt(data_path):
        """

        """
        mm_atom_type_sets_directory = str(str(data_path)+"/mm_atom_type_sets")
        mm_atom_properties_file = str(str(mm_atom_type_sets_directory)+"/mm_atom_properties.txt")
        print("Reading "+str(mm_atom_properties_file))
        mm_atom_types = {}
        file = open(mm_atom_properties_file,'r')
        lines = file.readlines()
        for line in lines[1:]:
          mm_atom_data = line.split()
          mm_atom_types.update({'name':mm_atom_data[0],'lj_wdepth':mm_atom_data[1],'lj_radius':mm_atom_data[2],'lj_3b_wdepth':mm_atom_data[3],'lj_3b_radius':mm_atom_data[4]})
        file.close()
        return(mm_atom_types)

def replace_mm_atom_properties_txt(data_path,mm_param_dict):
        """

        """
        mm_atom_type_sets_directory = str(str(data_path)+"/mm_atom_type_sets")
        mm_atom_properties_file = str(str(mm_atom_type_sets_directory)+"/mm_atom_properties.txt")
        file = open(mm_atom_properties_file,'r')
        lines = file.readlines()
        for line in lines[1:]:
          mm_atom_data = line.split()
          mm_atom_types.update({'name':mm_atom_data[0],'lj_wdepth':mm_atom_data[1],'lj_radius':mm_atom_data[2],'lj_3b_wdepth':mm_atom_data[3],'lj_3b_radius':mm_atom_data[4]})
        file.close()
        return(mm_atom_types)


def add_mm_atom_types(mm_param_dict):
        """
        Given a cgmodel and an 'mm' atom type parameter dictionary, this function writes the atoms to 'mm_atom_properties.txt' in the cg_pyrosetta database. 

        """
        mm_atom_type_sets_directory = str(str(data_path)+"/mm_atom_type_sets")
        mm_atom_properties_file = str(str(mm_atom_type_sets_directory)+"/mm_atom_properties.txt")
        mm_atom_types = {}
        file = open(mm_atom_properties_file,'a')
        for atom in mm_atom_types:
          file.write(str(mm_atom_types['name'])+"   "+str(mm_atom_types['lj_wdepth'])+"   "+str(mm_atom_types['lj_radius'])+"   "+str(mm_atom_types['lj_3b_wdepth'])+"   "+str(mm_atom_types['lj_3b_radius']))
        return

def get_param_dict_from_cgopenmm_cgmodel(cgmodel):
    """
    """
    particle_list = cgmodel.get_particle_list()
    mm_param_dict = {}
    for particle_index in range(len(particle_list)):
          particle = particle_list[particle_index]
          lj_wdepth = cgmodel.get_epsilon(particle_index)
          lj_radius = cgmodel.get_sigma(particle_index)*(2.0**(1/6))
          lj_3b_wdepth = lj_wdepth
          lj_3b_radius = lj_radius
          if 'X' in particle:
           particle_name = particle.replace('X','BB')
          if 'A' in particle:
           particle_name = particle.replace('A','SC')
          if 'X' in particle_name or 'A' in particle_name:
           print("ERROR: the particle names are not being re-assigned correctly for cg_pyrosetta")
           exit()

          mm_param_dict.update({'name':particle_name,'lj_wdepth':lj_wdepth,'lj_radius':lj_radius,'lj_3b_wdepth':lj_3b_wdepth,'lj_3b_radius':lj_3b_radius})

    return(mm_param_dict)


def assign_mm_atom_properties_with_cgopenmm_cgmodel(cgmodel):
        """
        """
        existing_mm_atom_data = read_mm_atom_properties_txt(data_path)
        print(existing_mm_atom_data)
        new_mm_atom_data = get_param_dict_from_cgopenmm_cgmodel(cgmodel)
        print(new_mm_atom_data)
        mm_atoms_to_replace = {}
        mm_atoms_to_add = {}
        for new_atom in new_mm_atom_data:
          new = True
          for existing_atom in existing_mm_atom_data:
            if existing_atom['name'] == new_atom['name']:
              new = False
              mm_atoms_to_replace.update({'name':new_atom['name'],'lj_wdepth':new_atom['lj_wdepth'],'lj_radius':new_atom['lj_radius'],'lj_3b_wdepth':new_atom['lj_3b_wdepth'],'lj_3b_radius':new_atom['lj_3b_radius']})
                                
            if new:
              mm_atoms_to_add.update({'name':new_atom['name'],'lj_wdepth':new_atom['lj_wdepth'],'lj_radius':new_atom['lj_radius'],'lj_3b_wdepth':new_atom['lj_3b_wdepth'],'lj_3b_radius':new_atom['lj_3b_radius']})

        if mm_atoms_to_add != {}:
          add_mm_atom_types(mm_atoms_to_add)
        if mm_atoms_to_replace != {}:
          replace_mm_atom_types(mm_atoms_to_replace)

        return

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
    
    with open(os.path.join(data_path, 'mm_atom_type_sets', 'mm_atom_properties.txt'), 'w') as f:
        # write header
        f.write('NAME    LJ_WDEPTH   LJ_RADIUS   LJ_3B_WDEPTH    LJ_3B_RADIUS\n')
        
        # rewrite updated parameters
        for param in atom_params_list:
            # print param line with specific formating
            f.write("%s%10.4f%10.4f%10.4f%10.4f\n" % (param[0], -float(param[3]), float(param[2]), -float(param[3]), float(param[2])))
            
    
        

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

    prev_torsion = [ torsion[0] for torsion in torsion_params_list]
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
            f.write("%s %s %s %s %.4f %.1i %.4f\n" % (atom_names[0], atom_names[1], atom_names[2], atom_names[3], float(param[1]), int(param[2]), float(param[3])))

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

    prev_angles = [ angle[0] for angle in angle_params_list]
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
                         angle,]
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
            f.write("%-4.4s%4.3s%4.3s%10.4f%10.4f\n" % (atom_names[0], atom_names[1], atom_names[2], float(param[1]), float(param[2])))
