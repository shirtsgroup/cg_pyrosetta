import os, shutil, datetime
import cg_pyrosetta
import pyrosetta
import mdtraj as md
from simtk import unit
import foldamers
from foldamers.utilities.iotools import write_pdbfile_without_topology
from foldamers.utilities.util import *
from foldamers.cg_model.cgmodel import *
from foldamers.cg_model.cgmodel import CGModel
# Identify the Rosetta database directory
pyrosetta_database_path = pyrosetta._rosetta_database_from_env()

def get_bonded_particle_list(residue):
        """
        Returns the list of names for particles that are bonded together in this residue.

        Parameters
        ----------

        residue: A dictionary containing information about the particles in a residue/monomer type

        Returns
        -------

        bonded_list: A list of the particles that are bonded together in this residue type.
        List([[atom_1_name(str),atom_2_name(str)]])

        """
        bonded_list = []
        atom_1_name = str("BB1")
        for backbone_bead in residue['backbone_length']:
          if backbone_bead != 0:
            atom_2_name = str("BB"+str(backbone_bead+1))
            bonded_list.append([atom_1_name,atom_2_name])
            atom_1_name = atom_2_name
            if backbone_bead in residue['sidechain_positions']:
              for sidechain_bead in residue['sidechain_length']:
                atom_2_name = str("SC"+str(sidechain_bead+1))
                bonded_list.append([atom_1_name,atom_2_name])
                atom_1_name = atom_2_name
              if residue['sidechain_length'] == 1:
                atom_2_name = str("VIRT")
                bonded_list.append([atom_1_name,atom_2_name])
          atom_1_name = str("BB"+str(backbone_bead+1))
        return(bonded_list)

def get_monomer_internal_coordinates(cgmodel,monomer_type):
        """
        Returns a list of internal coordinates for a monomer/residue, give n a cgmodel class object and a monomer_type dictionary. (Used to construct a residue .params file for PyRosetta.)

        Parameters
        ----------

        cgmodel: CGModel() class object

        monomer_type: A dictionary containing information about the particles in a residue/monomer type

        Returns
        -------

        internal_coordinates: A list of internal coordinates for 'monomer_type'
        List([strings])

        """
        # Building internal coordinates for residue .params files in PyRosetta using the example here:
        # https://www.rosettacommons.org/docs/latest/rosetta_basics/file_types/Residue-Params-file
        monomer_positions = cgmodel.positions[0:monomer_type['num_beads']]
        internal_coordinates = []
        bonded_particle_list = get_bonded_particle_list(monomer_type)

        torsion_index = 0
        # Build monomer internal coordinates by iterating over the list of torsions in the cgmodel
        for torsion in cgmodel.torsion_list:
         if all(torsion) < monomer_type['num_beads']:

          # Construct a list of atoms from which to build an internal coordinate line
          if torsion_index == 0:
           atom_a1,atom_a2,atom_a3,atom_a4 = cgmodel.particle_list[torsion[2]],cgmodel.particle_list[torsion[1]],cgmodel.particle_list[torsion[0]],cgmodel.particle_list[torsion[0]]
          else:
           atom_a1,atom_a2,atom_a3,atom_a4 = cgmodel.particle_list[torsion[3]],cgmodel.particle_list[torsion[2]],cgmodel.particle_list[torsion[1]],cgmodel.particle_list[torsion[0]]

          # Determine the bond length for this internal coordinate line
          if atom_a3 == atom_a4:
           bond_length = 0.0
          else:
           bond_length = get_bond_length_from_names(cgmodel,atom_a3,atom_a4)

          # Determine theta for this internal coordinate line
          if atom_a3 == atom_a4:
           theta = 0.0
          if atom_a2 == atom_a4:
           theta = 180.0
          if atom_a3 != atom_a4 and atom_a2 != atom_a4:
           theta = 180.0 - get_bond_angle_from_names(cgmodel,atom_a4,atom_a3,atom_a2)

          # Determine phi for this internal coordinate line
          if len(torsions) > len(set(torsions)):
           phi = 0.0
          else:
           phi = get_torsion_from_names(cgmodel,atom_a4,atom_a3,atom_a2,atom_a1)
          # Construct a line from the coordinate values chosen.
          line_list = ['ICOOR_INTERNAL    ',atom_a4,phi,theta,bond_length,atom_a3,atom_a2,atom_a1]
          line = '%18s%4s   %10s   %10s   %9s   %4s    %4s    %4s' % (line_list[i] for i in len(line_list))
          internal_coordinates.append(line)

         torsion_index = torsion_index + 1

        return(internal_coordinates)
   
def remove_existing_atom_types(atom_properties_file,list_of_atoms_to_remove):
        """
        Given a 'list_of_atoms_to_remove' and the path to a PyRosetta database file containing atom properties ('atom_properties.txt'), this function removes old atom types.        

        Parameters
        ----------

        atom_properties_file: Path to the PyRosetta database file containing atom type properties ('atom_properties.txt')

        list_of_atoms_to_remove: List of atom types to remove from the database.

        """
        file_obj = open(atom_properties_file,'r')
        lines = file_obj.readlines()
        file_obj.close()
        file_obj = open(atom_properties_file,'w')
        for line in lines:
          atom_type = line.split(' ')[0]
          if atom_type not in list_of_atoms_to_remove:
            file_obj.write(line)
        file_obj.close()
        return
 
def get_existing_atom_types(atom_properties_file):
        """
        Given the path to a PyRosetta database file, 'atom_properties.txt', this function reads the file and returns a list of atom type names.

        Parameters
        ----------

        atom_properties_file: Path to the PyRosetta database file containing atom type properties ('atom_properties.txt')

        Returns
        -------

        existing_atom_types_list: List of existing atom types.

        """
        existing_atom_types_list = []
        file_obj = open(file_name,'r')
        lines = file_obj.readlines()
        for line in lines:
          if line[0] != '#':
            atom_type = line.split(' ')[0]
            existing_atom_types_list.append(atom_type)
        return(existing_atom_types_list)

def write_mm_atom_properties_txt(cgmodel,list_of_atoms_to_add):
        """
        Given a cgmodel and a 'list_of_atoms_to_add', this function adds the atoms to 'mm_atom_properties.txt'.

        Parameters
        ----------

        cgmodel: CGModel() class object

        list_of_atoms_to_add: List of atom types to write to 'mm_atom_properties.txt'
        List([ strings ])

        """
        mm_atom_type_sets_directory = str(str(pyrosetta_database_path)+"mm_atom_type_sets/coarse_grain")
        if not os.path.exists(mm_atom_type_sets_directory): os.mkdir(mm_atom_type_sets_directory)
        atom_properties_file = str(str(atom_type_sets_directory)+"/mm_atom_properties.txt")
        if os.path.exists(atom_properties_file):
         existing_mm_atom_types = get_existing_mm_atom_types(atom_properties_file)
         file_obj = open(atom_properties_file,'a')
         for atom_type in list_of_atoms_to_add:
          if residue_type in existing_residue_types:
           print("WARNING: found an existing atom type with the same name in:\n")
           print(str(str(atom_properties_file)+"\n"))
           print("Removing the existing atom type definition from 'atom_properties.txt'")
           remove_existing_residue_types(file_name,[residue_type])
        else:
         file_obj = open(atom_properties_file,'w')
         file_obj.write("NAME  ATOM LJ_RADIUS LJ_WDEPTH LK_DGFREE LK_LAMBDA LK_VOLUME\n")
         file_obj.write("## Coarse grained residue types to follow\n")
        for atom_type in atom_type_list:
         if len(atom_type) > 4:
          print("ERROR: an atom type with a name longer than 4 characters has been defined for this model.\n")
          print("PyRosetta syntax requires that all atom types have names with for characters or less.")
          exit()
        particle_type = get_particle_type(cgmodel,particle_index=-1,particle_name=atom_type)
        sigma = get_sigma(cgmodel,particle_index=-1,particle_type=particle_type)
        epsilon = get_epsilon(cgmodel,particle_index=-1,particle_type=particle_type)
        lk_dgfree = 0.0000
        lk_lambda = 3.5000
        lk_volume = 0.0000
        symbol = 'X'
        comments = ""
        line_list = [atom_type,symbol,sigma,epsilon,lk_dgfree,lk_lambda,lk_volume,comments]
        line = '%4s     %1s %9s %9s %9s    %6s %9s %s' % (line_list[i] for i in len(line_list))
        file_obj.write(line)
        file_obj.close()
        return


def write_atom_properties_txt(cgmodel,list_of_atoms_to_add):
        """
        Given a cgmodel and a 'list_of_atoms_to_add', this functions writes the atoms to 'atom_properties.txt' in the PyRosetta database. 

        Parameters
        ----------

        cgmodel: CGModel() class object

        list_of_atoms_to_add: List of atom types to write to 'atom_properties.txt'
        List([ strings ])

        """
        atom_type_sets_directory = str(str(pyrosetta_database_path)+"atom_type_sets/coarse_grain")
        if not os.path.exists(atom_type_sets_directory): os.mkdir(atom_type_sets_directory)
        atom_properties_file = str(str(atom_type_sets_directory)+"/atom_properties.txt")
        if os.path.exists(atom_properties_file):
         existing_atom_types = get_existing_atom_types(atom_properties_file)
         file_obj = open(atom_properties_file,'a')
         for residue_type in residue_types_list:
          if residue_type in existing_residue_types:
           print("WARNING: found an existing atom type with the same name in:\n")
           print(str(str(atom_properties_file)+"\n"))
           print("Removing the existing atom type definition from 'atom_properties.txt'")
           remove_existing_residue_types(file_name,[residue_type])
        else:
         file_obj = open(atom_properties_file,'w')
         file_obj.write("NAME  ATOM LJ_RADIUS LJ_WDEPTH LK_DGFREE LK_LAMBDA LK_VOLUME\n")
         file_obj.write("## Coarse grained residue types to follow\n")
         for atom_type in atom_type_list:
          if len(atom_type) > 4:
           print("ERROR: an atom type with a name longer than 4 characters has been defined for this model.\n")
           print("PyRosetta syntax requires that all atom types have names with for characters or less.")
           exit()
         particle_type = get_particle_type(cgmodel,particle_index=-1,particle_name=atom_type)
         sigma = get_sigma(cgmodel,particle_index=-1,particle_type=particle_type)
         epsilon = get_epsilon(cgmodel,particle_index=-1,particle_type=particle_type)
         lk_dgfree = 0.0000
         lk_lambda = 3.5000
         lk_volume = 0.0000
         symbol = 'X'
         comments = ""
         line_list = [atom_type,symbol,sigma,epsilon,lk_dgfree,lk_lambda,lk_volume,comments]
         line = '%4s     %1s %9s %9s %9s    %6s %9s %s' % (line_list[i] for i in len(line_list))
        file_obj.write(line)
        file_obj.close()
        return

def build_patches(cgmodel):
        """
        Builds the patch

        Parameters
        ----------

        cgmodel: CGModel() class object

        particle_index: Index of the particle for which we would like to determine the type
        Type: int()

        Returns
        -------

        particle_type: 'backbone' or 'sidechain'
        Type: str()

        """
        residue_type_sets_directory = str(str(pyrosetta_database_path)+"residue_type_sets/coarse_grain")
        if not os.path.exists(residue_type_sets_directory): 
         os.mkdir(residue_type_sets_directory)
        if not os.path.exists(str(str(residue_type_sets_directory)+"/patches")):
         os.mkdir(str(str(residue_type_sets_directory)+"/patches"))
        patches_list = ['LOWER_TERMINUS_VARIANT','UPPER_TERMINUS_VARIANT']
        return

def remove_existing_residue_types(residue_types_txt_file,list_of_residue_types_to_remove):
        """
        Given a 'list_of_residue_types_to_remove', and a 'residue_types_txt_file', this function will open the residue types file, and remove the target residue types.

        Parameters
        ----------

        residue_types_txt_file: The path to a PyRosetta residue types file ('residue_types.txt')

        list_of_residue_types_to_remove: A list of residue type names to remove from 'residue_types.txt'

        """
        file_obj = open(residue_types_txt_file,'r')
        lines = file_obj.readlines()
        file_obj.close()
        file_obj = open(residue_types_txt_file,'w')
        for line in lines:
          if line[0] != '#':
            if '.params' in line:
              residue_type = line.split('/')[-1].split('.params')[0]
              if residue_type not in list_of_residue_types_to_remove:
                file_obj.write(line)
        file_obj.close() 
        return

def get_existing_residue_types(residue_types_txt_file):
        """
        Given a 'residue_types_txt_file', this function gets a list of residue types that are in the file.
        (Commented out lines/residue types are ignored.)

        Parameters
        ----------

        residue_types_txt_file: The path to a PyRosetta residue types fi
le ('residue_types.txt')

        Returns
        -------

        existing_residue_types_list: A list of residue types that are currently in 'residue_types.txt'.

        """
        existing_residue_types_list = []
        file_obj = open(residue_types_txt_file,'r')
        lines = file_obj.readlines()
        for line in lines:
          if line[0] != '#':
            if '.params' in line:
              residue_type = line.split('/')[-1].split('.params')[0]
              existing_residue_types_list.append(residue_type)
        return(existing_residue_types_list)

def write_residue_types_txt(residue_types_list,residue_type_sets_directory):
        """
        Given a 'residue_type_sets_directory', and a 'residue_types_list', this function writes the residue types to 'residue_types.txt' (creating the file if it doesn't already exist, and appending the residue types to the file if it does exist.

        Parameters
        ----------

        residue_types_sets_directory: The path to a directory containing a list of residue types

        residue_types_list: A list of residue types to write/add to 'residue_types.txt'.

        """

        file_name = str(str(residue_type_sets_directory)+"residue_types.txt")
        if os.path.exists(file_name):
         existing_residue_types = get_existing_residue_types(file_name)
         file_obj = open(file_name,'a')
         file_obj.write("\n")
         file_obj.write("## Coarse grained residue types to follow\n")
         for residue_type in residue_types_list:
          if residue_type in existing_residue_types:
           print("WARNING: found an existing residue type set in the PyRosetta database folder:\n")
           print(str(str(residue_type_sets_directory.split('/residue_types')[0])+"\n"))
           print("with the same name as a residue you wish to add:"+str(str(residue_type)+"\n"))
           print("Removing the existing residue type set's path from 'residue_type_sets.txt'")
           remove_existing_residue_types(file_name,[residue_type])
        else:
           file_obj = open(file_name,'w')
           file_obj.write("## Define atom and mm type sets\n")
           file_obj.write("TYPE_SET_MODE coarse_grain\n")
           file_obj.write("ATOM_TYPE_SET coarse_grain\n")
           file_obj.write("ELEMENT_SET coarse_grain\n")
           file_obj.write("MM_ATOM_TYPE_SET coarse_grain\n")
           file_obj.write("ORBITAL_TYPE_SET coarse_grain\n")
           file_obj.write("##\n")
           file_obj.write("\n")
        file_obj.write("## Coarse grained residue types to follow\n")
        for residue_type in residue_types_list:  
         file_obj.write("residue_types/"+str(residue_type)+".params")
        file_obj.close()
        return

def build_params_files(cgmodel):
        """
        Given a cgmodel class object, this function writes '.params' files for all unique residue types in the class.

        Parameters
        ----------

        cgmodel: CGModel() class object

        """
        print(pyrosetta_database_path)
        residue_type_sets_directory = str(str(pyrosetta_database_path)+"/chemical/residue_type_sets/coarse_grain")
        if not os.path.exists(residue_type_sets_directory): os.mkdir(residue_type_sets_directory)
        if not os.path.exists(str(str(residue_type_sets_directory)+"/residue_types")): 
         os.mkdir(str(str(residue_type_sets_directory)+"/residue_types"))
        residue_code_options = ['A','B','C','D','E','F']
        residue_code_list = []
        residue_code_index = 0
        residue_code = residue_code_options[residue_code_index]
        for residue in cgmodel.monomer_types:
         residue_type_name = str("CG"+str(residue['backbone_length'])+str(residue['sidechain_length']))
         while residue_code in residue_code_list:
          residue_code_index = residue_code_index + 1
          residue_code = residue_code_options[residue_code_index]
         file_name = str(str(residue_code)+".param")
         file_obj = open(file_name,'w')
         file_obj.write("# Rosetta residue topology file\n")
         file_obj.write("# Authors: Lenny T. Fobe, Garrett A. Meek\n")
         file_obj.write("# ( Research group of Professor Michael R. Shirts )\n")
         file_obj.write("# Dept. of Chemical and Biological Engineering\n")
         file_obj.write("# University of Colorado Boulder\n")
         file_obj.write("# This file was written on: "+str(datetime.datetime.today())+"\n")
         file_obj.write("\n")
         file_obj.write("NAME "+str(residue_type_name)+"\n")
         file_obj.write("IO_STRING "+str(residue_type_name)+" "+str(residue_code)+"\n")
         file_obj.write("TYPE POLYMER\n")
         file_obj.write("\n")
         file_obj.write(str("VARIANT\n"))
         file_obj.write("\n")
         for backbone_bead in range(residue['backbone_length']):
           atom_name = str("BB"+str(backbone_bead+1))
           file_obj.write("ATOM "+str(atom_name)+" VIRT X 0.0\n")
           if backbone_bead in [residue['sidechain_positions']]:
             for sidechain_bead in range(residue['sidechain_length']):
               atom_name = str("SC"+str(sidechain_bead+1))
               file_obj.write("ATOM "+str(atom_name)+" VIRT X 0.0\n")
             if residue['sidechain_length'] == 1:
               file_obj.write("ATOM VIRT VIRT X 0.0\n")
         file_obj.write("\n")
         upper_connect = str("BB"+str(residue['backbone_length']))
         file_obj.write("LOWER CONNECT BB1\n")
         file_obj.write("UPPER_CONNECT "+str(upper_connect)+"\n")
         atom_1_name = str("BB1")
         for backbone_bead in range(residue['backbone_length']):
           if backbone_bead != 0:
             atom_2_name = str("BB"+str(backbone_bead+1))
             file_obj.write("BOND "+str(atom_1_name)+" "+str(atom_2_name)+"\n")
             atom_1_name = atom_2_name
           if backbone_bead in [residue['sidechain_positions']]:
             for sidechain_bead in range(residue['sidechain_length']):
               atom_2_name = str("SC"+str(sidechain_bead+1))
               file_obj.write("BOND "+str(atom_1_name)+" "+str(atom_2_name)+"\n")
               atom_1_name = atom_2_name
             if residue['sidechain_length'] == 1:
               atom_2_name = str("VIRT")
               file_obj.write("BOND "+str(atom_1_name)+" "+str(atom_2_name)+"\n")
           atom_1_name = str("BB"+str(backbone_bead+1))    
         file_obj.write("\n")
         file_obj.write("FIRST_SIDECHAIN_ATOM SC1\n")
         file_obj.write("PROPERTIES\n")
         file_obj.write("\n")
         file_obj.close()
         shutil.move(file_name,str(str(residue_type_sets_directory)+"/residue_types/"+str(file_name)))
        write_residue_types_txt(residue_code_list,residue_type_sets_directory)
        return

def build_scorefxn(cgmodel):
        """
        Given a cgmodel class object, this function uses its definitions to build a PyRosetta scoring function.
        (Used to confirm energetic agreement between PyRosetta and OpenMM for identical model parameters and an identical structural ensemble.)

        Parameters
        ----------

        cgmodel: CGModel() class object

        Returns
        -------

        scorefxn: A PyRosetta ScoreFunction() class object, containing all scoring function components that are defined/requested within the cgmodel() class object.

        """
        scorefxn = pyrosetta.ScoreFunction()
        scorefxn.set_weight(pyrosetta.rosetta.core.scoring.fa_atr, 1)
        scorefxn.set_weight(pyrosetta.rosetta.core.scoring.fa_rep, 1)
        scorefxn.set_weight(pyrosetta.rosetta.core.scoring.fa_intra_atr, 1)
        scorefxn.set_weight(pyrosetta.rosetta.core.scoring.fa_intra_rep, 1)
#        scorefxn.set_weight(pyrosetta.rosetta.core.scoring.mm_twist, 0.1)
        return(scorefxn)

def compare_pose_scores(scorefxn,pose_1,pose_2,compare_pdb_sequence=True):
        """
        Given two PyRosetta poses ('pose_1' and 'pose_2') and a scoring function ('scorefxn'), this function compares their scores to confirm agreement.
        (Can also be used to confirm that the scores for poses generated from the sequence and from a PDB file are identical, in order to validate our procedure for building PyRosetta '.params' files from the cgmodel class object.)

        Parameters
        ----------

        scorefxn: A PyRosetta scoring function

        pose_1: A PyRosetta pose (The pose built from a PDB file if 'compare_pdb'=True)

        pose_2: A PyRosetta pose (The pose build from a sequence if 'compare_pdb'=True)

        compare_pdb_sequence: A logical variable determining whether this score comparison is meant to determine agreement between poses built from a PDB file and a sequence.

        """
        if compare_pdb_sequence:
         pdb_score = scorefxn.show(pose_1)
         sequence_score = scorefxn.show(pose_2)
         if sequence_score != pdb_score:
          print("WARNING: Getting different scores when the pose is built from a sequence and from a PDB file.\n")
          print("The score when building a pose from the sequence is: "+str(sequence_score)+"\n")
          print("The score when building a pose from the PDB file is: "+str(pdb_score)+"\n")
         else:
          print("The scores are identical when a pose is built from a PDB file and from the polymer sequence.\n")
        else:
         score_1 = scorefxn.show(pose_1)
         score_2 = scorefxn.show(pose_2)
         if score_1 != score_2:
          print("The scores for these poses are different:\n")
          print("The score for pose 1 is:"+str(score_1)+"\n")
          print("The score for pose 2 is:"+str(score_2)+"\n")
         else:
          print("The scores for these poses are the same.\n") 
        return

def compare_openmm_energy_pyrosetta_score(cgmodel):
        from cg_openmm.cg_mm_tools.cg_openmm import build_mm_simulation

        """
        Given a cgmodel class object, this function determines if PyRosetta and OpenMM give the same score/energy with identical model settings.

        Parameters
        ----------

        cgmodel: Coarse grained model class object.

        """

        build_params_files(cgmodel)
        pyrosetta.init()

        pyrosetta_sequence = ''.join([str('['+str(monomer['monomer_name'])+']') for monomer in cgmodel.sequence])
        # Build a PyRosetta pose
        pose = pyrosetta.pose_from_sequence(pyrosetta_sequence)
        # Define a PyRosetta scoring function
        scorefxn = build_scorefxn(cgmodel)
        # get the PyRosetta score
        score = scorefxn.show(pose)
        # Build an OpenMM simulation object so that we can calculate the energy using a simulation Context()
        simulation = build_mm_simulation(cgmodel.topology,cgmodel.system,cgmodel.positions,temperature=temperature,simulation_time_step=simulation_time_step,total_simulation_time=total_simulation_time,output_pdb=output_pdb,output_data=output_data,print_frequency=print_frequency)
        # Obtain a state for our simulation context
        energy = simulation.context.getState(getEnergy=True).getPotentialEnergy()
        print("The PyRosetta score is: "+str(score)+"\n")
        print("The OpenMM potential energy is: "+str(energy)+"\n")
        return
        

# Set values for the parameters in our coarse grained model:
polymer_length=8
backbone_length=[1,2]
sidechain_length=1
sidechain_positions=[0]
mass = unit.Quantity(10.0,unit.amu)
sigma = unit.Quantity(2.4,unit.angstrom)
bond_length = unit.Quantity(1.0,unit.angstrom)
epsilon = unit.Quantity(0.5,unit.kilocalorie_per_mole)
# charge = unit.Quantity(0.0,unit.elementary_charge)

# Define PDB files to test our PDB writing ability
openmm_pdb_file = 'test_1_1_openmm.pdb'
rosetta_pdb_file = 'test_1_1_rosetta.pdb'

# Build a coarse grained model
cgmodel = basic_cgmodel(polymer_length=polymer_length,backbone_length=backbone_length,sidechain_length=sidechain_length,sidechain_positions=sidechain_positions,mass=mass,bond_length=bond_length,sigma=sigma,epsilon=epsilon)
#write_pdbfile_without_topology(cgmodel,openmm_pdb_file)
pyrosetta_sequence = ''.join([str('['+str(monomer['monomer_name'])+']') for monomer in cgmodel.sequence])
# Compare OpenMM and PyRosetta energies
# (This function is also where we initialize new residue/monomer
#  types in the PyRosetta database.)
compare_openmm_energy_pyrosetta_score(cgmodel)
pose_from_sequence = pyrosetta.pose_from_sequence(pyrosetta_sequence,'coarse_grain')
# Test our ability to write a PDB file using our pose and new residue type sets.
pyrosetta.rosetta.core.io.pdb.dump_pdb(pose,rosetta_pdb_file)
# Test our ability to read a pose from the PDB file we wrote
pose_from_pdb = pyrosetta.pose_from_pdb(rosetta_pdb_file)
# Define scorefunction terms
pyrosetta_scorefxn = build_scorefxn(cgmodel)
# Compare poses built from a PDB file and from the polymer sequence
compare_pose_scores(pyrosetta_scorefxn,pose_from_pdb,pose_from_sequence,compare_pdb_sequence=True)

exit()

