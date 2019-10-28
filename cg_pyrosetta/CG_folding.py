import numpy as np
import cg_pyrosetta.CG_movers as CG_movers
import os
import sys
import math

current_path = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.abspath(current_path + '/../PyRosetta4.modified'))
import pyrosetta
pyrosetta.init()

np.random.seed()


class CGFoldingAlgorithm():

    def __init__(self,
                 sequence,
                 BBB_angle=120,
                 BBBB_dihe=180,
                 file_name='outputs/traj.pdb',
                 energy_graph_output=False):
        """
        folding object used for easily implementing different
        movers into a single folding algorithm.

        Arguments
        ---------

        sequence : str
            Sequence of CG residues
        BBB_angle : float
            Desired angle of all B-B-B angles. Generalizes to all backbone models (not working)
        BBBB_angle : float
            Desired dihedral of all B-B-B-B torsion angles. Generalizes to all backbone models (not working)


        """
        # Build CG model and set desired initial angles
        self.pose = pyrosetta.pose_from_sequence(sequence, auto_termini=False)
        self.energy_graph_output = energy_graph_output
        # self.pose = self.set_BBB_angles(self.pose, BBB_angle)
        # self.pose = self.set_BBBB_dihe(self.pose, BBBB_dihe)
        # PyMOL mover, if wanting to visualize
        self.pymol = pyrosetta.PyMOLMover()
        self.pymol.apply(self.pose)

        randomizer = CG_movers.randomizeBackBone(self.pose)
        randomizer.apply(self.pose)

        self.pymol.apply(self.pose)
        # Building PDBTrajWriter object, used for writing multiple structures
        # to a single file
        self.PDB_writer = pyrosetta.rosetta.protocols.canonical_sampling.PDBTrajectoryRecorder()
        # self.PDB_writer.apply(self.pose) # write initial structure
        self.PDB_writer.file_name('outputs/traj.pdb')
        self.PDB_writer.stride(100)

        # Define scorefunction terms
        self.scorefxn = pyrosetta.ScoreFunction()
        self.scorefxn.set_weight(pyrosetta.rosetta.core.scoring.fa_rep, 1)
        self.scorefxn.set_weight(pyrosetta.rosetta.core.scoring.fa_atr, 1)
        self.scorefxn.set_weight(pyrosetta.rosetta.core.scoring.fa_intra_atr, 1)
        self.scorefxn.set_weight(pyrosetta.rosetta.core.scoring.fa_intra_rep, 1)
        self.scorefxn.set_weight(pyrosetta.rosetta.core.scoring.mm_twist, 1)
        self.scorefxn.set_weight(pyrosetta.rosetta.core.scoring.mm_bend, 1)
        # self.scorefxn.set_weight(pyrosetta.rosetta.core.scoring.mm_lj_inter_rep, 1) segfaults beware!
        # self.scorefxn.set_weight(pyrosetta.rosetta.core.scoring.mm_lj_inter_atr, 1)
        # self.scorefxn.set_weight(pyrosetta.rosetta.core.scoring.mm_lj_intra_rep, 1)
        # self.scorefxn.set_weight(pyrosetta.rosetta.core.scoring.mm_lj_intra_atr, 1)

        # Build standard CG 1-1 movers
        self.small = CG_movers.CGSmallMover(self.pose)
        self.shear = CG_movers.CGShearMover(self.pose)
        self.small_angle = CG_movers.CGSmallAngleMover(self.pose)

        # Build minimization movers
        self.mini = pyrosetta.rosetta.protocols.minimization_packing.MinMover()
        self.mini.min_type('lbfgs_armijo_nonmonotone')
        self.movemap = pyrosetta.MoveMap()
        self.mini.score_function(self.scorefxn)
        # for atom in self.small_angle.bb_atoms:
        #     self.movemap.set(pyrosetta.rosetta.core.id.DOF_ID(atom , pyrosetta.rosetta.core.id.THETA), True)
        self.movemap.set_bb_true_range(1, self.pose.size())
        self.mini.movemap(self.movemap)

        # Build MC object + Trial Mover (empty for now)
        self.mc = pyrosetta.MonteCarlo(self.pose, self.scorefxn, 1)
        self.trial_mc = pyrosetta.TrialMover()

        # Building variable to store various folding algorithms
        self.folding_protocols = {}

        # Adding a default mover
        self.build_fold_alg('default')
        self.add_folding_move('default', pyrosetta.RepeatMover(self.small, 10))
        self.add_folding_move('default', pyrosetta.RepeatMover(self.shear, 10))
        self.add_folding_move('default', pyrosetta.RepeatMover(self.mini, 10))

        # If writing a trajectory file or pymol visualization is desired uncomment these lines

        # self.add_folding_move('default', self.PDB_writer)
        # self.add_folding_move('default', self.pymol)

    def set_BBB_angles(self, pose, angle):
        """
        set initial B-B-B angle for each residue to a certain value
        pose : pyrosetta.Pose() object
        angle : angle in degrees

        works for any object with 1-bead backbone
        """
        for i in range(2, pose.size()):
            pose.conformation().set_bond_angle(pyrosetta.AtomID(1, int(i-1)),
                                               pyrosetta.AtomID(1, int(i)),
                                               pyrosetta.AtomID(1, int(i+1)),
                                               angle*np.pi/180)

        return(pose)

    def set_BBBB_dihe(self, pose, angle):
        """
        set initial B-B-B angle for each residue to a certain value
        pose : pyrosetta.Pose() object
        angle : angle in degrees

        works for any object with 1-bead backbone
        """
        for i in range(2, pose.size()-1):
            pose.conformation().set_torsion_angle(pyrosetta.AtomID(1, int(i-1)),
                                                  pyrosetta.AtomID(1, int(i)),
                                                  pyrosetta.AtomID(1, int(i+1)),
                                                  pyrosetta.AtomID(1, int(i+2)),
                                                  angle*np.pi/180)

        return(pose)

    def random_BBBB_dihe(self, pose):
        """
        set random initial B-B-B-B angle for each residue
        pose : pyrosetta.Pose() object
        angle : angle in degrees

        works for any object with 1-bead backbone
        """
        for i in range(2, pose.size()-1):
            angle = np.random.rand()*2*np.pi
            pose.conformation().set_torsion_angle(pyrosetta.AtomID(1, int(i-1)),
                                                  pyrosetta.AtomID(1, int(i)),
                                                  pyrosetta.AtomID(1, int(i+1)),
                                                  pyrosetta.AtomID(1, int(i+2)),
                                                  angle)

        return(pose)

    def build_fold_alg(self, name):
        """
        Add a new entry to the folding protocols dictionary
        Can be used to create custom sequences of movers

        name : str of name of folding algorithm
        """
        self.folding_protocols[name] = pyrosetta.SequenceMover()

    def add_folding_move(self, name, mover):
        """
        Add new mover to a entry in the folding protocols dictionary
        name : str of folding algorithm
        mover : mover object to be added to this folding sequence (e.g. small, shear, mc_trial)
        """
        self.folding_protocols[name].add_mover(mover)

    def build_trial_mc_alg(self, mover):
        """
        Write any sequence mover (stored in 'folding_protocols') to a trial object for
        MC simulations. Depends on MC mover, so change any mc parameters (kT, pose, etc.)
        before running this.
        mover: mover object to be added to TrialMover object
        """
        self.trial_mc = pyrosetta.TrialMover(mover, self.mc)

    def run_folding_alg(self, name, iter):
        """
        Runs any given mover object iter times as an MC trial steps useful for
        running multiple steps with the same parameters.
        name : str of folding algorithm (key of self.folding_algorithm)
        iter : integer of times to repeat TrialMC object
        """
        print('Building MC Trial Mover...')
        self.build_trial_mc_alg(self.folding_protocols[name])
        run = pyrosetta.RepeatMover(self.trial_mc, iter)
        print('Folding...')
        run.apply(self.pose)
        self.mc.set_last_accepted_pose(self.mc.lowest_score_pose())

    def run_anneal_fold(self, name, iter, kt_range):
        """
        Runs any specified mover object iter times as an MC trial step, while
        also anealling kT to achieve extremely small changes towards the end
        of the simulation.
        name : str of folding algorithm (key of self.folding_algorithm)
        iter : interger of times to repeat TrialMC object
        kt_range : array of kT values to run folding algorithm

        Note will run len(kt_range)*iter trial MC steps
        """
        if self.energy_graph_output:
            e_graph = []
        for kt in kt_range:

            # Updat kt in MC object
            self.mc.set_temperature(kt)

            # Rebuild trail mc object with new kT value
            self.build_trial_mc_alg(self.folding_protocols[name])

            # Build RepeatMover with MC trial to iterate iter times
            run = pyrosetta.RepeatMover(self.trial_mc, iter)

            old_energy = self.scorefxn(self.mc.lowest_score_pose())
            if self.energy_graph_output:
                e_graph.append(old_energy)
            new_energy = None
            counter = 0
            while old_energy != new_energy or counter < 10:
                #  new_energy == None or not math.isclose(old_energy, new_energy) and
                counter += 1
                print('Folding at T =', kt, '...', 'Rep:', counter)
                old_energy = self.scorefxn(self.mc.lowest_score_pose())
                run.apply(self.pose)
                self.mc.show_counters()
                new_energy = self.scorefxn(self.mc.lowest_score_pose())
                print('Old Energy:', old_energy, 'New Energy:', new_energy)
                if self.energy_graph_output:
                    e_graph.append(new_energy)
                if old_energy != new_energy:
                    counter = 0

        if self.energy_graph_output:
            e_graph = np.array(e_graph)
            count = 0
            while os.path.exists(name+"_"+str(iter)+"_"+str(count)+".npy"):
                count += 1
            np.save(name+"_"+str(iter)+"_"+str(count)+".npy", e_graph)


def main():
    for s in range(10):
        kt_i = 100
        kt_anneal = [kt_i*(0.9)**i for i in range(50)]
        obj = CGFoldingAlgorithm('X[CG11]X[CG11]X[CG11]X[CG11]X[CG11]X[CG11]X[CG11]X[CG11]X[CG11]X[CG11]X[CG11]X[CG11]X[CG11]X[CG11]X[CG11]',
                                 file_name='outputs/traj_anneal.pdb')
        obj.add_folding_move('default', obj.pymol)
        obj.PDB_writer.stride(2000)
        obj.small.angle = 180
        obj.shear.angle = 180
        obj.add_folding_move('default', obj.PDB_writer)
        obj.run_anneal_fold('default', 150, kt_anneal)
        obj.mc.lowest_score_pose().dump_pdb('outputs/11_lowest_energy_'+str(s)+'.pdb')

        lj_scorefxn = pyrosetta.ScoreFunction()

        lj_scorefxn.set_weight(pyrosetta.rosetta.core.scoring.mm_twist, 0.1)
        lj_scorefxn.set_weight(pyrosetta.rosetta.core.scoring.mm_lj_inter_rep, 1)
        lj_scorefxn.set_weight(pyrosetta.rosetta.core.scoring.mm_lj_inter_atr, 1)
        lj_scorefxn.set_weight(pyrosetta.rosetta.core.scoring.mm_lj_intra_rep, 1)
        lj_scorefxn.set_weight(pyrosetta.rosetta.core.scoring.mm_lj_intra_atr, 1)

        with open('outputs/11_lowest_energy_'+str(s)+'.pdb', 'a') as f:
            f.write('# LJ Energy : '+str(lj_scorefxn(obj.mc.lowest_score_pose()))
                    )


if __name__ == '__main__':
    main()
