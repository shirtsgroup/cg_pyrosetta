import signac
import cg_pyrosetta
import flow
from flow import FlowProject

flow.status_parallelization = "none"

@FlowProject.label
def set_parameters(job):
    # set parameters
    cg_pyrosetta.change_parameters.changeAngleParameters(
        {
         'CG1 CG1 CG1' : [1500, job.sp.bbb_angle],
         'CG2 CG1 CG1' : [1500, (360 - job.sp.bbb_angle)/2]
        },
        angle_file= "parameters/mm_atom_type_sets/mm_bond_angle_params.txt"
    )

@FlowProject.operation
@FlowProject.label
@FlowProject.pre(set_parameters)
@FlowProject.post.isfile("minimum.pdb")
def run_mc_simulation(job):
    cg_pyrosetta.pyrosetta.init(
                      "--add_atom_types fa_standard parameters/atom_properties.txt " +
                      "--add_mm_atom_type_set_parameters fa_standard parameters/mm_atom_type_sets/mm_atom_properties.txt " +
                      "--extra_mm_params_dir parameters/mm_atom_type_sets"
                      )
    # Build Annealer Parameters
    annealer_params = cg_pyrosetta.CG_monte_carlo.\
        CGMonteCarloAnnealerParameters(n_inner = 1000,
                                       t_init = 5,
                                       anneal_rate = 0.9,
                                       n_anneals = 30,
                                       annealer_criteron = cg_pyrosetta.CG_monte_carlo.Repeat10Convergence,
                                       traj_out = "mc-min_traj.pdb",
                                       mc_output = True,
                                       mc_traj = True,
                                       traj_freq = 250, 
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
    pose = cg_pyrosetta.pyrosetta.pose_from_sequence("X[CG11x3:CGLower]X[CG11x3]X[CG11x3]X[CG11x3]X[CG11x3:CGUpper]")
    change_lengths = cg_pyrosetta.CG_movers.setBondLengths(pose, {"BB1 BB2":job.sp.bb_length, "BB2 BB3":job.sp.bb_length, "BB3 BB1":job.sp.bb_length})
    change_lengths.apply(pose)

    # Build Minimizer
    mini = cg_pyrosetta.pyrosetta.rosetta.protocols.minimization_packing.MinMover()
    mini.min_type('lbfgs_armijo_nonmonotone')
    mini.score_function(energy_function)

    movemap = cg_pyrosetta.pyrosetta.MoveMap()
    movemap.set_bb_true_range(1, pose.size())
    movemap.set_branches(True)
    movemap.set(cg_pyrosetta.pyrosetta.rosetta.core.id.THETA, True)
    movemap.set(cg_pyrosetta.pyrosetta.rosetta.core.id.PHI, True)

    print("Using the following MoveMap:")
    print(movemap)

    mini.movemap(movemap)

    # Build Sequence Mover for MC object
    sequence_mover_factory = cg_pyrosetta.CG_monte_carlo.SequenceMoverFactory(pose, {"mini" : mini})
    sequence_mover = sequence_mover_factory.build_seq_mover(
        {
            "small_dihe" : 1,
            "small_angle" : 1,
            "mini" : 50,
        }
    )

    cg_annealer = cg_pyrosetta.CG_monte_carlo.CGMonteCarloAnnealer(
        seq_mover=sequence_mover,
        score_function=energy_function,
        pose = pose,
        param_file_object = annealer_params
    )

    # Run Annealer
    cg_annealer.run_schedule()
    min_pose = cg_annealer._cg_mc_sim.get_minimum_energy_pose()
    min_pose.dump_pdb("minimum.pdb")
    
@FlowProject.operation
@FlowProject.pre.isfile("minimum.pdb")
def low_temperature_mc(job):

    pose = cg_pyrosetta.pyrosetta.pose_from_file("minimum.pdb")

    # temporary output
    print(pose)


if __name__ == '__main__':
    FlowProject().main()