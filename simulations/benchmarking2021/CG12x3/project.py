import signac
import cg_pyrosetta
import os
import flow
import time
from flow import FlowProject

# flow.status_parallelization = "none"

@FlowProject.operation
@flow.directives(fork=True)
@FlowProject.post.isfile("minimum.pdb")
def run_mc_simulation(job):
    start = time.time()
    os.chdir(job.ws)
    cg_pyrosetta.init()
    # Build Annealer Parameters
    # Annealing simulation with 500000
    annealer_params = cg_pyrosetta.CG_monte_carlo.\
        CGMonteCarloAnnealerParameters(n_inner = 1000,
                                       t_init = 10,
                                       anneal_rate = 0.9,
                                       n_anneals = 50,
                                       annealer_criteron = cg_pyrosetta.CG_monte_carlo.Repeat10Convergence,
                                       mc_output = True,
                                       out_freq = 500,
    )

    # Build Energy Function
    energy_function = cg_pyrosetta.CG_monte_carlo.EnergyFunctionFactory().build_energy_function(
        {
            "mm_twist" : 1,
            "mm_bend" : 1,
            "fa_atr" : 1,
            "fa_rep" : 1,
            "fa_intra_rep" : 1,
            "fa_intra_atr" : 1,
        }
    )

    # Pose to be folded
    pose = cg_pyrosetta.pyrosetta.pose_from_sequence("X[CG12x3:CGLower]" + "X[CG12x3]" * (job.sp.nmer - 2) + "X[CG12x3:CGUpper]")

    # Build Minimizer
    mini = cg_pyrosetta.pyrosetta.rosetta.protocols.minimization_packing.MinMover()
    mini.min_type('lbfgs_armijo_nonmonotone')
    mini.score_function(energy_function)

    movemap = cg_pyrosetta.pyrosetta.MoveMap()
    movemap.set_bb_true_range(1, pose.size())
    movemap.set_branches(True)
    movemap.set(cg_pyrosetta.pyrosetta.rosetta.core.id.THETA, True)
    movemap.set(cg_pyrosetta.pyrosetta.rosetta.core.id.PHI, True)

    # print("Using the following MoveMap:")
    # print(movemap)

    mini.movemap(movemap)

    # Build Sequence Mover for MC object
    sequence_mover_factory = cg_pyrosetta.CG_monte_carlo.SequenceMoverFactory(pose, {"mini" : mini})
    sequence_mover = sequence_mover_factory.build_seq_mover(
        {
            "small_dihe" : 1,
            "small_angle" : 1,
            "mini" : 1,
        }
    )

    cg_annealer = cg_pyrosetta.CG_monte_carlo.CGMonteCarloAnnealer(
        seq_mover=sequence_mover,
        score_function=energy_function,
        pose = pose,
        param_file_object = annealer_params
    )

    # Setup Configuration/Energy observer for saving minimum energy structures
    struct_obs = cg_pyrosetta.CG_monte_carlo.StructureObserver(cg_annealer.get_mc_sim())
    energy_obs = cg_pyrosetta.CG_monte_carlo.EnergyObserver(cg_annealer.get_mc_sim())
    cg_annealer.registerObserver(struct_obs)
    cg_annealer.registerObserver(energy_obs)


    # Run Annealer
    cg_annealer.run_schedule()
    min_pose = cg_annealer._cg_mc_sim.get_minimum_energy_pose()
    print("Writing structure to:")
    print(job.fn("minimum.pdb"))
    min_pose.dump_pdb(job.fn("minimum.pdb"))

    end = time.time()

    t_rep = end - start
    print(t_rep, "seconds")
    job.data['timing'] = t_rep
    

if __name__ == '__main__':
    FlowProject().main()
