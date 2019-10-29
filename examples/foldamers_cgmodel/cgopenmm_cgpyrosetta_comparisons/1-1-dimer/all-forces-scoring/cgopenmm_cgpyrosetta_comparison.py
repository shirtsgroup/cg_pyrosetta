import cg_pyrosetta
import os
from simtk import unit
from cg_openmm.simulation.tools import *
from cg_openmm.build.cg_build import *
from cg_openmm.utilities.iotools import *
from cg_openmm.utilities.util import *
from foldamers.cg_model.cgmodel import CGModel
#pyrosetta.init()
#pyrosetta_database_path = pyrosetta._rosetta_database_from_env()

def build_params_files(cgmodel):
 for monomer in cgmodel.sequence:
  file_name = str(str(monomer['monomer_name'])+".params")
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
 return

def build_cgmodel(rosetta_scoring):
 # Set values for the parameters in our coarse grained model:
 polymer_length=2
 backbone_lengths=[1]
 sidechain_lengths=[1]
 sidechain_positions=[0]
 include_bond_forces=False # NOTE: By default bonds are constrained to their equilibrium lengths, even when this variable is set to False, unless 'constrain_bonds'=False
 constrain_bonds=True
 include_bond_angle_forces=True
 include_nonbonded_forces=True
 include_torsion_forces=True

 # Particle properties
 mass = 100.0 * unit.amu
 masses = {'backbone_bead_masses': mass, 'sidechain_bead_masses': mass}
 bond_length = 1.0 * unit.angstrom
 bond_lengths = {'bb_bb_bond_length': bond_length,'bb_sc_bond_length': bond_length,'sc_sc_bond_length': bond_length}
 bond_force_constant = 0.0 * unit.kilocalorie_per_mole / unit.nanometer / unit.nanometer
 bond_force_constants = {'bb_bb_bond_k': bond_force_constant, 'bb_sc_bond_k': bond_force_constant, 'sc_sc_bond_k': bond_force_constant}
 r_min = bond_length
 sigma = r_min / (2.0**(1/6))
 sigmas = {'bb_sigma': sigma,'sc_sigma': sigma}
 epsilon = 10.0 * unit.kilocalorie_per_mole
 epsilons = {'bb_eps': epsilon,'sc_eps': epsilon}
 exclusions = True

 # Bond angle properties
 bond_angle_force_constant = 2 * unit.kilocalorie_per_mole / unit.radian / unit.radian
 bond_angle_force_constants = {'bb_bb_sc_angle_k': bond_angle_force_constant}
 equil_bond_angle = 90.0*(3.141/180.0)
 equil_bond_angles = {'bb_bb_sc_angle_0': equil_bond_angle}

 # Torsion properties
 torsion_force_constant = 3 * unit.kilocalorie_per_mole / unit.radian / unit.radian
 torsion_force_constants = {'sc_bb_bb_sc_torsion_k': torsion_force_constant}
 torsion_periodicity = 1
 torsion_periodicities = {'sc_bb_bb_sc_period': torsion_periodicity}
 equil_torsion_angle = 0.0*(3.141/180.0)
 equil_torsion_angles = {'sc_bb_bb_sc_torsion_0': equil_torsion_angle}

 # Get positions from a local PDB file written by PyRosetta, and modified (by hand) to have a geometry where the nonbonded interactions are easy to evaluate
 cgmodel = CGModel(polymer_length=polymer_length,backbone_lengths=backbone_lengths,sidechain_lengths=sidechain_lengths,sidechain_positions=sidechain_positions,masses=masses,sigmas=sigmas,epsilons=epsilons,bond_lengths=bond_lengths,include_nonbonded_forces=include_nonbonded_forces,include_bond_forces=include_bond_forces,include_bond_angle_forces=include_bond_angle_forces,include_torsion_forces=include_torsion_forces,rosetta_scoring=rosetta_scoring,exclusions=exclusions,constrain_bonds=constrain_bonds,equil_torsion_angles=equil_torsion_angles,torsion_force_constants=torsion_force_constants,torsion_periodicities=torsion_periodicities,equil_bond_angles=equil_bond_angles,bond_angle_force_constants=bond_angle_force_constants)
 pyrosetta_sequence = ''.join([str('X['+str(monomer['monomer_name'])+']') for monomer in cgmodel.sequence])
 res_file_list = list(set([str(str(monomer['monomer_name'])+".params") for monomer in cgmodel.sequence]))
 res_file_list = " ".join(res_file_list)
 pyrosetta.init(extra_options = str("-extra_res "+str(res_file_list)))
 pyrosetta_database_path = pyrosetta._rosetta_database_from_env()
 mm_atom_data = read_mm_atom_properties_txt(pyrosetta_database_path)
 particle_list = get_particle_list_from_cgmodel(cgmodel)
 assign_mm_atom_properties(cgmodel,pyrosetta_database_path)
 exit()
 pose = pyrosetta.pose_from_sequence(pyrosetta_sequence)
 for residue_index in range(pose.total_residue()):
  print(pose.residue(residue_index))
 exit()
 pose.dump_pdb("test_pyrosetta.pdb")
 cgmodel.positions = PDBFile("test_pyrosetta.pdb").getPositions()
 cgmodel.topology = build_topology(cgmodel)
 return(cgmodel,pose)

def build_scorefxn(cgmodel,mm=False):
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
        if mm:
          # Create new MM atom types
          file = open("mm_atom_properties.txt","w")
          file.write("NAME    LJ_WDEPTH       LJ_RADIUS       LJ_3B_WDEPTH    LJ_3B_RADIUS\n")
          file.write("BB1      -0.20000       1.00        -0.20000       1.000\n")
          file.write("SC1      -0.20000       1.00        -0.20000       1.000\n")
          file.close()

          mm_atom_type_set = pyrosetta.rosetta.core.chemical.MMAtomTypeSet("CG")
          mm_atom_type_set.read_file("mm_atom_properties.txt")
          pyrosetta.rosetta.core.scoring.mm.MMLJLibrary(mm_atom_type_set)
          #exit()
          #pyrosetta.rosetta.core.scoring.mm.MMLJLibrary() 

          # Create a new MM score function
          file = open("mm.wts","w")
          file.write("UNFOLDED_ENERGIES_TYPE UNFOLDED_MM_STD\n")
          file.write("mm_lj_intra_rep 1.0\n")
          file.write("mm_lj_intra_atr 1.0\n")
          file.write("mm_lj_inter_rep 1.0\n")
          file.write("mm_lj_inter_atr 1.0\n")
          file.write("mm_twist 0.0")
          file.close()

          scorefxn = pyrosetta.create_score_function("mm.wts")
          #print(scorefxn)
          #exit()
          #print(pyrosetta.rosetta.core.scoring.mm.MMLJLibrary.lookup(1))
          #exit()
        
        else:
         scorefxn.set_weight(pyrosetta.rosetta.core.scoring.fa_elec,1)
         scorefxn.set_weight(pyrosetta.rosetta.core.scoring.fa_rep, 1)
         scorefxn.set_weight(pyrosetta.rosetta.core.scoring.fa_atr, 1)
         scorefxn.set_weight(pyrosetta.rosetta.core.scoring.fa_intra_atr, 1)
         scorefxn.set_weight(pyrosetta.rosetta.core.scoring.fa_intra_rep, 1)
         cg_pyrosetta.change_parameters.changeTorsionParameters(
          {'CG1 CG1 CG1 CG1':[0,0,0],
          'CG2 CG1 CG1 CG2':[0,0,0],
          'CG2 CG1 CG1 CG1':[0,0,0],
          'X CG2 CG1 CG1':[0,0,0]}
        )

         cg_pyrosetta.change_parameters.changeAngleParameters(
          {'CG1 CG1 CG1':[0,0],
          'CG2 CG1 CG1':[0,0],
          'CG1 CG1 CG2':[0,0],
          'X CG2 CG1':[0,0]}
        )
        return(scorefxn)

def compare_openmm_energy_pyrosetta_score(cgmodel,pose=None,mm=False):
        """
        Given a cgmodel class object, this function determines if PyRosetta and OpenMM give the same score/energy with identical model settings.

        Parameters
        ----------

        cgmodel: Coarse grained model class object.

        """

        # Build a PyRosetta pose
        #pose = pyrosetta.pose_from_pdb("init.pdb")
        # Define a PyRosetta scoring function
        scorefxn = build_scorefxn(cgmodel,mm=mm)
        # Get the PyRosetta score
        score = scorefxn(pose)
        # Get the cg_openmm energy
        energy = get_mm_energy(cgmodel.topology,cgmodel.system,cgmodel.positions).in_units_of(unit.kilocalorie_per_mole)
        # Obtain a state for our simulation context
        print("The PyRosetta score is: "+str(score))
        #print("The bond list is: "+str(cgmodel.bond_list))
        print("The OpenMM potential energy is: "+str(energy))
        file = open("energies.dat","w")
        file.write("The OpenMM nonbonded interaction list is: "+str(cgmodel.nonbonded_interaction_list)+"\n")
        file.write("The distances between these particles are: "+str([distance(cgmodel.positions[interaction[0]],cgmodel.positions[interaction[1]]) for interaction in cgmodel.nonbonded_interaction_list])+"\n")
        file.write("The LJ sigma value is: "+str(cgmodel.get_sigma(0))+"\n")
        file.write("The LJ epsilon value (in kJ/mol) is: "+str(cgmodel.get_epsilon(0).in_units_of(unit.kilojoules_per_mole))+"\n")
        file.write("The LJ energy calculated by hand is: "+str([lj_v(cgmodel.positions[interaction[0]],cgmodel.positions[interaction[1]],cgmodel.get_sigma(0),cgmodel.get_epsilon(0)) for interaction in cgmodel.nonbonded_interaction_list])+"\n")
        file.write("(in kJ/mol): "+str([lj_v(cgmodel.positions[interaction[0]],cgmodel.positions[interaction[1]],cgmodel.get_sigma(0),cgmodel.get_epsilon(0)).in_units_of(unit.kilojoules_per_mole) for interaction in cgmodel.nonbonded_interaction_list])+"\n")
        file.write("The PyRosetta score is: "+str(score)+"\n")
        file.write("The OpenMM potential energy is: "+str(energy)+"\n")
        file.close()
        return

# Build a cg_openmm cgmodel
mm=True
rosetta_scoring=False
cgmodel,pose = build_cgmodel(rosetta_scoring)

# Compare Rosetta score and OpenMM energy
compare_openmm_energy_pyrosetta_score(cgmodel,pose=pose,mm=mm)

# Perform cgopenmm simulation
temperature = 1.0 * unit.kelvin
total_simulation_time = 100.0 * unit.picosecond
simulation_time_step = 5.0 * unit.femtosecond
print_frequency=5
run_simulation(cgmodel,'output',total_simulation_time,simulation_time_step,temperature,print_frequency)
#plot_simulation_results("simulation.dat",str(os.getcwd()+'/output'),simulation_time_step,total_simulation_time)

exit()
