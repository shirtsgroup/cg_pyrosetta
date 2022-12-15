# This file will contain the object used for running simple CG MC simulations
# Will require :
# 1) ScoreFucntion object
# 2) Pose
# 3) Set of Movers (SequenceMover object)
# 4) kT
# 5) n_steps
# 6) output_freq?
import numpy as np
import copy
from abc import ABC, abstractmethod
import os
import sys
import warnings
import cg_pyrosetta.CG_movers


current_path = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.abspath(current_path + '/../PyRosetta4.modified'))

import pyrosetta

class Subject(ABC):
    def __init__(self):
        self._observers = []

    def registerObserver(self, observer):
        self._observers.append(observer)

    def removeObserver(self, observer):
        self._observers.remove(observer)

    def notifyObservers(self):
        for observer in self._observers:
            observer.update()

class Observer(ABC):
    @abstractmethod
    def update(self):
        pass



class EnergyObserver(Observer):
    """
    Observer design pattern to record energy from MC simulation
    """
    def __init__(self, subject, write_file = True, file_name = "energies.txt"):
        """
        EnergyObserver constructor

        Paramters
        ---------
        subject : Subject
            Object the energy observer will observe for energies
        write_file : bool
            Flag wether to write energies to file
        file_name : string
            output name of energy file
        """
        self.energies = []
        self.subject = subject
        self.write_file = write_file
        self.file_name = file_name
        if os.path.isfile(file_name):
            os.remove(file_name)

    def update(self):
        """
        Function pulls energy from subject and writes to file, depending on `self.write_file` flag.
        """
        energy = self.subject.get_energy()
        self.energies.append(energy)
        if self.write_file is True:
            if os.path.isfile(self.file_name):
                with open(self.file_name, 'a') as f:
                    f.write(str(self.subject._total_trials) + ",")
                    f.write(str(energy) + "\n")
            else:
                with open(self.file_name, 'w') as f:
                    f.write(str(self.subject._total_trials) + ",")
                    f.write(str(energy) + "\n")



class StructureObserver(Observer):
    """
    Observer pattern to record structure files of MC simulation
    """
    def __init__(self, subject, write_file = True, file_name = "structure.txt", pdb_file = "trajectory.pdb"):
        """
        StructureObserver constructor

        Parameters
        ----------
        subject : Subject
            Object to write structures from
        write_file : bool
            Flag of whether to structures to file
        file_name : str
            Output file name of structure files
        pdb_file : str
            Output file name of the pdb trajectory
        """
        self.structures = []
        self.subject = subject
        self.write_file = write_file
        self.file_name = file_name
        if self.write_file:
            if os.path.isfile(pdb_file):
                os.remove(pdb_file)
            if os.path.isfile(file_name):
                os.remove(file_name)
            self.traj_writer = pyrosetta.rosetta.protocols.canonical_sampling.PDBTrajectoryRecorder()
            self.traj_writer.file_name(pdb_file)
            self.traj_writer.stride(1)


    def update(self):
        """
        This function pulls structures from subject and writes to file, depending on `self.write_file` flag.
        """
        structure = self.subject.pose.clone()
        self.structures.append(structure)
        if self.write_file is True:
            self.traj_writer.apply(structure)
            if os.path.isfile(self.file_name):
                with open(self.file_name, 'a') as f:
                    f.write(str(self.subject._total_trials) + ",")
                    f.write(str(len(self.structures)) + "\n")
            else:
                with open(self.file_name, 'w') as f:
                    f.write("TimeStep,Frame\n")
                    f.write(str(self.subject._total_trials) + ",")
                    f.write(str(len(self.structures)) + "\n")

class KTObserver(Observer):
    """
    Observer object to write kT values from a simulation
    """
    def __init__(self, subject, write_file = True, file_name = "kt_schedule.txt"):
        """
        KTObserver constructor

        Parameters
        ----------
        subject : Subject
            Object to pull kT values from
        """
        self.kts = []
        self.subject = subject
        self.write_file = write_file
        self.file_name = file_name
        if os.path.isfile(file_name):
            os.remove(file_name)

    def update(self):
        """
        This function pull kT values from subject and writes to file, depending on `self.write_file` flag.
        """
        kt = self.subject.kT
        self.kts.append(kt)
        if self.write_file is True:
            if os.path.isfile(self.file_name):
                with open(self.file_name, 'a') as f:
                    f.write(str(self.subject._total_trials) + ",")
                    f.write(str(kt) + "\n")
            else:
                with open(self.file_name, 'w') as f:
                    f.write(str(self.subject._total_trials) + ",")
                    f.write(str(kt) + "\n")

class AcceptanceRatioObserver(Observer):
    """
    Observer object to calculate acceptance ratio of MC simulation objects
    """
    def __init__(self, subject, write_file = True, file_name = "acc_ratio.txt"):
        """
        AcceptanceRatioObserver constructor

        Parameters
        ----------
        subject : Subject
            Subject object (CGMonteCarlo) to pull acceptance ratios from
        write_file : bool
            Flag on whether to write acceptance ratios to a file
        file_name : str
            Output file name
        """
        self.acc_ratios = []
        self.subject = subject
        self.write_file = write_file
        self.file_name = file_name
        if os.path.isfile(file_name):
            os.remove(file_name)

    def update(self):
        """
        This function pull acceptance ratios from the CGMonteCarlo object
        """
        acc_ratio = self.subject.get_accept_ratio(reset = False)
        if self.write_file is True:
            if os.path.isfile(self.file_name):
                with open(self.file_name, 'a') as f:
                    f.write(str(self.subject._total_trials) + ",")
                    f.write(str(acc_ratio) + "\n")
            else:
                with open(self.file_name, 'w') as f:
                    f.write(str(self.subject._total_trials) + ",")
                    f.write(str(acc_ratio) + "\n")


class CGMonteCarlo(Subject):
    """
    Object used to run a PyRosetta Monte Carlo simulation
    """
    def __init__(self,
                 pose: object = None,  # pyrosetta.rosetta.core.pose.Pose, Pose of starting structure
                 score: object = None,  # pyrosetta.ScoreFunction, # scorefunction from pyrosetta
                 seq_mover: object = None,  # pyrosetta.rosetta.protocols.moves.SequenceMover, # sequence of moves between MC evaluations
                 n_steps: int = 1000000,
                 kT: float = 1,
                 output: bool = True,
                 out_freq: int = 500,):
        """
        CGMonteCarlo constructor

        Parameters
        ----------
        pose : pyrosetta.rosetta.core.pose.Pose
            Pose object to perform MC simulation with
        score : pyrosetta.ScoreFunction
            Scorefunction object to score pose with
        seq_mover : pyrosetta.rosetta.protocols.moves.SequenceMover
            Sequence mover with desired moves for MC simulation
        n_steps : int
            Number of steps to run MC simulation for
        kT : float
            Temperature to evaluate energy differences in during the metropolis criteron
        output : bool
            Flag to print output to terminal
        out_freq : int
            Frequency of steps between outputs
        """
        super().__init__()

        # initialize input values
        self.pose = pose
        self._score = score
        self.seq_mover = seq_mover
        self._kT = kT
        self.n_steps = n_steps
        self._output = output
        self._out_freq = out_freq
        self._total_trials = 0 # Could implement as a property

        if self._output is False:
            self._out_freq = n_steps

        # Build MC Object
        self.mc = pyrosetta.MonteCarlo(self.pose, self._score, self._kT)
        self.mc_trial = pyrosetta.TrialMover(self.seq_mover, self.mc)

        if self._output:
            self.pymol = pyrosetta.PyMOLMover()
            print("Initial Energy :", self.get_energy())

    @property
    def out_freq(self):
        """
        Get or set output frequency
        """
        return(self._out_freq)

    @out_freq.setter
    def out_freq(self, new_out_freq):
        self._out_freq = new_out_freq
    
    @property
    def kT(self):
        """
        Get or set MC simulation kT
        """
        return(self._kT)

    @kT.setter
    def kT(self, kT):
        self._kT = kT
        self.mc.set_temperature(self.kT)
        self.mc_trial = pyrosetta.TrialMover(self.seq_mover, self.mc)

    def get_energy(self):
        """
        Evaluate energy of scorefunction with given pose
        """
        return(self._score(self.pose))

    # def __call__():
    def run(self):
        """
        Run MC simulation for specified number of steps (self.n_steps)
        """
        rep_mover = pyrosetta.RepeatMover(self.mc_trial, self._out_freq)
        for _ in range(int(self.n_steps/self._out_freq)):
            rep_mover.apply(self.pose)
            self._total_trials += self._out_freq
            if self._output is True:
                print("Step :", self._total_trials)
                print("Energy : ", self.get_energy())
                self.pymol.apply(self.pose)
            self.notifyObservers()

    def get_accept_ratio(self, reset = True):
        """
        Calculate the acceptance ratio from the MC simulation

        Parameters
        ----------
        reset : bool (True)
            Flag to reset counts each time acceptance ratio is calculated.
        """
        acc_ratio = self.mc_trial.acceptance_rate()
        if reset:
            self.mc.reset_counters()
        return acc_ratio

    def get_pose(self):
        """
        Get the Pose object in the MC simulation
        """
        return(self.pose)

    def get_minimum_energy_pose(self):
        """
        Get the minimum energy pose found in the MC simulation
        """
        return(self.mc.lowest_score_pose())


class CGMonteCarloAnnealer:
    """
    Annealer object used to make simulated annealing MC simulations
    using CGMonteCarlo objects
    """
    def __init__(self,
                 seq_mover: object = None,
                 score_function: object = None,
                 pose: object = None,
                 param_file_object: object = None,
                 ):
        """
        CGMonteCarloAnnealr constructor

        Parameters
        ----------
        seq_mover : pyrosetta.rosetta.protocols.moves.SequenceMover
            Sequence mover with desired moves for MC simulation
        score_function : pyrosetta.ScoreFunction
            Scorefunction object to score pose with
        pose : pyrosetta.rosetta.core.pose.Pose
            Pose object to perform MC simulation with
        param_file_object : CGMonteCarloAnnealerParameters
            Parameter object used to store relevant parameters for annealer
            objects
        """
        self.seq_mover = seq_mover
        self.score_function = score_function
        self.pose = pose
        self.param_file_object = param_file_object
        self.convergence_criterea = param_file_object.annealer_criteron
        self._cg_mc_sim = CGMonteCarlo(self.pose, self.score_function,
                                       self.seq_mover, self.param_file_object.n_inner,
                                       param_file_object.t_init,
                                       output=param_file_object.mc_output, 
                                       out_freq = param_file_object.out_freq,
                                       )
        self.kt_anneals = copy.deepcopy(param_file_object.t_anneals)

    def run_schedule(self):
        """
        Runs annealing cycles using parameters defined in the CGMonteCarloAnnealerParameter object
        """
        for kt in self.kt_anneals:
            print("Current kT = ", kt)
            self._cg_mc_sim.kT = kt
            convergence_criterea = self.convergence_criterea()
            while not convergence_criterea(self._cg_mc_sim):
                self._cg_mc_sim.run()
            self.kt_anneals = self.kt_anneals[1:]

    def get_mc_sim(self):
        """
        Returns the CGMonteCarlo object used within this object
        """
        return self._cg_mc_sim
            
    def _get_score_function(self):
        """
        Gets the score function object
        """
        return self.score_function

    def _get_seq_mover(self):
        """
        Gets the sequence mover object 
        """
        return self.seq_mover

    def _add_to_schedule(self, kT):
        """
        Add new temperature values to the annealing schedule

        Parameters
        ----------
        kT : float
        """
        self.kt_anneals.append(kT)

    def registerObserver(self, observer):
        """
        Add new observers to observer list

        Parameters
        ----------
        observer : Observer
        """
        self._cg_mc_sim.registerObserver(observer)

    def removeObserver(self, observer):
        """
        Remove observer from observer list

        Parameters
        ----------
        observer : Observer
        """
        self._cg_mc_sim.removeObserver(observer)

class CGMonteCarloDynamicAnnealer:
    """
    A MC annealing object that dynamically adjusts temperature according
    to specific criteria objects
    """
    def __init__(self,
                seq_mover: object = None,
                score_function: object = None,
                pose: object = None,
                dynamic_param_file_object: object = None,
                ):
        """
        CGMonteCarloDynamicAnnealer constructor

        Parameters
        ----------
        seq_mover : pyrosetta.rosetta.protocols.moves.SequenceMover
            Sequence mover with desired moves for MC simulation
        score_function : pyrosetta.ScoreFunction
            Scorefunction object to score pose with
        pose : pyrosetta.rosetta.core.pose.Pose
            Pose object to perform MC simulation with
        dynamic_param_file_object : CGMonteCarloDynamicAnnealerParameters
            Parameter object used to store relevant parameters for the dynamic 
            annealer objects
        """
        self.seq_mover = seq_mover
        self.score_function = score_function
        self.pose = pose
        self.dynamic_params = dynamic_param_file_object
        self._cg_mc_sim = CGMonteCarlo(self.pose, self.score_function,
                                       self.seq_mover, self.dynamic_params.n_inner,
                                       dynamic_param_file_object.t_init,
                                       output=dynamic_param_file_object.mc_output, 
                                       out_freq = dynamic_param_file_object.out_freq,
                                       )
    
    def estimate_starting_kt(self, base = 10, factor = 1):
        """
        Function to estimate a good starting kT value. Gives an initial
        temperature in the order of magnitude of the initial energy of the
        MC simulation. Setting `base=1` will return the absolute value
        of the initial energy.
        
        Parameters
        ----------
        base : float (10)
            Base value to use for order of magnitude estimate
        factor : float (1)
            Factor by which to multiply the base estimate
        """
        energy = self._cg_mc_sim.get_energy()
        if base == 1:
            estimate_kt = factor * np.abs(energy)
        else:
            estimate_kt = factor * base ** np.ceil(np.log(np.abs(self._cg_mc_sim.get_energy())) / np.log(base))
        print("Setting kT to", estimate_kt)
        self._cg_mc_sim.kT = estimate_kt

    def run_annealing(self):
        """
        Run the annealing schedule
        """
        for i in range(self.dynamic_params.n_cycles):
            print("Current kT: ", self._cg_mc_sim.kT)
            self._cg_mc_sim.run()
            self.adjust_kt()
        

    def adjust_kt(self):
        """
        Adjust kT value based on acceptance ratio. If the acceptance
        ratio is too large, kT will be decreased by a factor of
        `self.dynamic_params.anneal_rate`, if the acceptance ratio
        is too small, kT will be increased by a factor of 
        `self.dynamic_params.anneal_rate`.
        """
        acc_ratio = self._cg_mc_sim.get_accept_ratio(reset = True)
        print("Cycle acceptance ratio: " + str(acc_ratio))
        print("Target acceptance ratio: " + str(self.dynamic_params.target_ratio))
        if acc_ratio > self.dynamic_params.target_ratio + self.dynamic_params.ratio_tolerance:
            self._cg_mc_sim.kT *= self.dynamic_params.anneal_rate
        if acc_ratio < self.dynamic_params.target_ratio - self.dynamic_params.ratio_tolerance:
            self._cg_mc_sim.kT /= self.dynamic_params.anneal_rate
        
    def get_mc_sim(self):
        """
        Get MC simulation object
        """
        return self._cg_mc_sim
            
    def _get_score_function(self):
        """
        Get scoreterm object
        """
        return self.score_function

    def _get_seq_mover(self):
        """
        Get sequence mover object
        """
        return self.seq_mover

    def registerObserver(self, observer):
        """
        Add new observers to observer list

        Parameters
        ----------
        observer : Observer
        """
        self._cg_mc_sim.registerObserver(observer)

    def removeObserver(self, observer):
        """
        Remove observer from observer list

        Parameters
        ----------
        observer : Observer
        """
        self._cg_mc_sim.removeObserver(observer)
    

class Convergence(ABC):
    def __call__(self, mc_sim):
        pass

class Repeat10Convergence(Convergence):
    """
    Convergence object, used to evaluate the convergence of an MC simulation
    at a given condition. Here convergence is met with 10 iteration
    """
    def __init__(self):
        self.counter = 0

    def __call__(self, mc_sim):
        if self.counter < 10:
            self.counter += 1
            return False
        else:
            return True


class Repeat1Convergence(Convergence):
    """
    Convergence object, used to evaluate the convergence of an MC simulation
    at a given condition. Here convergence is met with 1 iteration
    """
    def __init__(self):
        self.counter = 0

    def __call__(self, mc_sim):
        if self.counter < 1:
            self.counter += 1
            return False
        else:
            return True




# class PoseAdapter(ABC):
#     @abstractmethod
#     def get_pose(self):
#         pass

class CGMonteCarloAnnealerParameters:
    """
    Data object for storing how a fixed schedule MC annealing simulation will be run
    """
    def __init__(self, n_inner, t_init, anneal_rate, n_anneals, annealer_criteron, mc_output, out_freq):
        self.n_inner = n_inner
        self.t_init = t_init
        self.anneal_rate = anneal_rate
        self.n_anneals = n_anneals
        self.annealer_criteron = annealer_criteron # Need a new object for this, unclear
        self.mc_output = mc_output
        self.out_freq = out_freq
        # list of annealing temps
        self.t_anneals = [t_init*(anneal_rate)**n for n in range(n_anneals)]

class CGMonteCarloDynamicAnnealerParameters:
    """
    Data object for storing how a fixed schedule MC annealing simulation will be run
    """
    def __init__(self, n_inner, t_init, anneal_rate, max_anneal_cycles, target_ratio, ratio_tolerance, mc_output, out_freq):
        self.n_inner = n_inner
        self.t_init = t_init
        self.anneal_rate = anneal_rate
        self.n_cycles = max_anneal_cycles
        self.mc_output = mc_output
        self.out_freq = out_freq
        self.target_ratio = target_ratio
        self.ratio_tolerance = ratio_tolerance
        # list of annealing temps

class SequenceMoverFactory:
    """
    Factory object for building pyrosetta.SequenceMover() objects with
    specific movers
    """
    def __init__(self, pose, extra_movers_dict = None):
        self.methods = {
            'small_dihe': cg_pyrosetta.CG_movers.CGSmallMover(pose),
            'small_angle': cg_pyrosetta.CG_movers.CGSmallAngleMover(pose),
            'shear_dihe': NotImplementedError(),
            'sc_small_dihe': NotImplementedError(),
            'sc_small_angle': cg_pyrosetta.CG_movers.CGSmallAngleMover(pose),
            'pymol' : cg_pyrosetta.pyrosetta.PyMOLMover(),
        }
        if extra_movers_dict:
            self.methods.update(extra_movers_dict)
            

    def build_seq_mover(self, mover_freq_map):
        seq_mover = pyrosetta.SequenceMover()
        movers = [pyrosetta.RepeatMover(self.methods[mover], mover_freq_map[mover]) for mover in  mover_freq_map.keys()]
        for i, mover in enumerate(mover_freq_map.keys()):
            print(mover)
            if mover in self.methods.keys():
                seq_mover.add_mover(movers[i])
            else:
                warnings.warn("Unimplemented Mover : "+mover+"\n Skipping mover", UserWarning)

        return(seq_mover)


class EnergyFunctionFactory:
    """
    Factory object for building pyrosetta.ScoreFunction()
    objects with specific score terms
    """
    def __init__(self):
        self.methods = {}
        for method in dir(pyrosetta.rosetta.core.scoring):
            self.methods[method] = eval('pyrosetta.rosetta.core.scoring.'+method)

    def build_energy_function(self, score_term_weight_map):
        score = pyrosetta.ScoreFunction()

        for term in score_term_weight_map.keys():
            if term in self.methods.keys():
                score.set_weight(eval('pyrosetta.rosetta.core.scoring.'+term), score_term_weight_map[term])
            else:
                warnings.warn("Energy Term not implemented :"+term+"\n Skipping term", UserWarning)

        return(score)
