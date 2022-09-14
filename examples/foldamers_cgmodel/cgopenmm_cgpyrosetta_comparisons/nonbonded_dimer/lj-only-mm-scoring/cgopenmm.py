import cg_pyrosetta
import pyrosetta
from simtk import unit
from cg_openmm.simulation.tools import *
from cg_openmm.build.cg_build import *
from cg_openmm.utilities.iotools import *
from cg_openmm.utilities.util import *
from foldamers.cg_model.cgmodel import CGModel
# Identify the Rosetta database directory
pyrosetta_database_path = pyrosetta._rosetta_database_from_env()

def build_cgmodel(rosetta_scoring):
 # Set values for the parameters in our coarse grained model:
 polymer_length=2
 backbone_lengths=[1]
 sidechain_lengths=[0]
 sidechain_positions=[0]
 include_bond_forces=False # NOTE: By default bonds are constrained to their equilibrium lengths, even when this variable is set to False, unless 'constrain_bonds'=False
 constrain_bonds=False
 include_bond_angle_forces=False
 include_nonbonded_forces=True
 include_torsion_forces=False

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
 epsilon = 0.2 * unit.kilocalorie_per_mole
 epsilons = {'bb_eps': epsilon,'sc_eps': epsilon}
 exclusions = False

 # Get positions from a local PDB file written by PyRosetta, and modified (by hand) to have a geometry where the nonbonded interactions are easy to evaluate
 positions = get_positions_from_pdbfile("openmm_init.pdb")

 # Build a coarse grained model
 cgmodel = CGModel(polymer_length=polymer_length,backbone_lengths=backbone_lengths,sidechain_lengths=sidechain_lengths,sidechain_positions=sidechain_positions,masses=masses,sigmas=sigmas,epsilons=epsilons,bond_lengths=bond_lengths,include_nonbonded_forces=include_nonbonded_forces,include_bond_forces=include_bond_forces,include_bond_angle_forces=include_bond_angle_forces,include_torsion_forces=include_torsion_forces,rosetta_scoring=rosetta_scoring,exclusions=exclusions,constrain_bonds=constrain_bonds,positions=positions)

 return(cgmodel)

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
        scorefxn = pyrosetta.ScoreFunction()
        
        if mm:
         scorefxn.set_weight(pyrosetta.rosetta.core.scoring.mm_lj_intra_rep,1)
         scorefxn.set_weight(pyrosetta.rosetta.core.scoring.mm_lj_intra_atr,1)
         scorefxn.set_weight(pyrosetta.rosetta.core.scoring.mm_lj_inter_atr,1)
         scorefxn.set_weight(pyrosetta.rosetta.core.scoring.mm_lj_inter_rep,1)
        
        else:
         scorefxn.set_weight(pyrosetta.rosetta.core.scoring.fa_elec,1)
         scorefxn.set_weight(pyrosetta.rosetta.core.scoring.fa_rep, 1)
         scorefxn.set_weight(pyrosetta.rosetta.core.scoring.fa_atr, 1)
         scorefxn.set_weight(pyrosetta.rosetta.core.scoring.fa_intra_atr, 1)
         scorefxn.set_weight(pyrosetta.rosetta.core.scoring.fa_intra_rep, 1)
        cg_pyrosetta.parameters.changeTorsionParameters(
          {'CG1 CG1 CG1 CG1':[0,0,0],
          'CG2 CG1 CG1 CG2':[0,0,0],
          'CG2 CG1 CG1 CG1':[0,0,0],
          'X CG2 CG1 CG1':[0,0,0]}
        )

        cg_pyrosetta.parameters.changeAngleParameters(
          {'CG1 CG1 CG1':[0,0],
          'CG2 CG1 CG1':[0,0],
          'CG1 CG1 CG2':[0,0],
          'X CG2 CG1':[0,0]}
        )
#        scorefxn.set_weight(pyrosetta.rosetta.core.scoring.mm_twist, 0.1)
        return(scorefxn)

def compare_openmm_energy_pyrosetta_score(cgmodel,mm=False):
        """
        Given a cgmodel class object, this function determines if PyRosetta and OpenMM give the same score/energy with identical model settings.

        Parameters
        ----------

        cgmodel: Coarse grained model class object.

        """

        # Build a PyRosetta pose
        #pose = pyrosetta.pose_from_pdb("init.pdb")
        # Define a PyRosetta scoring function
        #scorefxn = build_scorefxn(cgmodel,mm=mm)
        # Get the PyRosetta score
        #score = scorefxn(pose)
        # Get the cg_openmm energy
        energy = get_mm_energy(cgmodel.topology,cgmodel.system,cgmodel.positions).in_units_of(unit.kilocalorie_per_mole)
        # Obtain a state for our simulation context
        #print("The PyRosetta score is: "+str(score))
        print("The bond list is: "+str(cgmodel.bond_list))
        print("The OpenMM potential energy is: "+str(energy))
        file = open("energies.dat","w")
        file.write("The nonbonded interaction list is: "+str(cgmodel.nonbonded_interaction_list)+"\n")
        file.write("The distances between these particles are: "+str([distance(cgmodel.positions[interaction[0]],cgmodel.positions[interaction[1]]) for interaction in cgmodel.nonbonded_interaction_list])+"\n")
        file.write("The LJ energy calculated by hand is: "+str([lj_v(cgmodel.positions[interaction[0]],cgmodel.positions[interaction[1]],cgmodel.get_sigma(0),cgmodel.get_epsilon(0)) for interaction in cgmodel.nonbonded_interaction_list])+"\n")
        #file.write("The PyRosetta score is: "+str(score)+"\n")
        file.write("The OpenMM potential energy is: "+str(energy)+"\n")
        file.close()
        return

# Build a cg_openmm cgmodel
mm=True
rosetta_scoring=False
cgmodel = build_cgmodel(rosetta_scoring)

# Compare Rosetta score and OpenMM energy
compare_openmm_energy_pyrosetta_score(cgmodel,mm=mm)

#exit()
