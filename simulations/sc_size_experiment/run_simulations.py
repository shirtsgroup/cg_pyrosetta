import signac
import cg_pyrosetta
import os
import flow
from flow import FlowProject

# flow.status_parallelization = "none"

@FlowProject.label
def parameters_are_set(job):
    return job.isfile("job_status.txt")

@FlowProject.operation
@FlowProject.post(parameters_are_set)
def set_parameters(job):

    os.chdir(job.ws)
    # set parameters
    print("Changing parameters in", os.path.abspath(""))
    cg_pyrosetta.change_parameters.changeAtomParameters(
        {
         'CG2' : ['X', job.sp.sc_size, 0.2],
        },
        atom_types_path = "parameters/atom_properties.txt",
        mm_atom_types_path = "parameters/mm_atom_type_sets/mm_atom_properties.txt"
    )

    with open(job.fn("job_status.txt"), "w") as f:
        f.write("parameters set\n")


@FlowProject.operation
@FlowProject.pre(parameters_are_set)
@FlowProject.post.isfile("minimum.pdb")
def run_mc_simulation(job):
    os.chdir(job.ws)
    cg_pyrosetta.init()
    # Build Annealer Parameters
    annealer_params = cg_pyrosetta.CG_monte_carlo.\
        CGMonteCarloAnnealerParameters(n_inner = 1000,
                                       t_init = 5,
                                       anneal_rate = 0.9,
                                       n_anneals = 30,
                                       annealer_criteron = cg_pyrosetta.CG_monte_carlo.Repeat10Convergence,
                                       traj_out = job.fn("mc-min_traj.pdb"),
                                       mc_output = True,
                                       mc_traj = True,
                                       out_freq = 250, 
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
    change_lengths = cg_pyrosetta.CG_movers.setBondLengths(pose, {"BB1 SC1":job.sp.sc_size, "BB2 SC2":job.sp.sc_size, "BB3 SC3":job.sp.sc_size})
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

    # print("Using the following MoveMap:")
    # print(movemap)

    mini.movemap(movemap)

    # Build Sequence Mover for MC object
    sequence_mover_factory = cg_pyrosetta.CG_monte_carlo.SequenceMoverFactory(pose, {"mini" : mini})
    sequence_mover = sequence_mover_factory.build_seq_mover(
        {
            "small_dihe" : 1,
            "small_angle" : 1,
            "mini" : 10,
        }
    )

    cg_annealer = cg_pyrosetta.CG_monte_carlo.CGMonteCarloAnnealer(
        seq_mover=sequence_mover,
        score_function=energy_function,
        pose = pose,
        param_file_object = annealer_params
    )

    # Setup Configuration/Energy observer for saving minimum energy structures
    min_energy_confs = cg_pyrosetta.CG_monte_carlo.MinEnergyConfigObserver(cg_annealer.get_mc_sim())
    cg_annealer.registerObserver(min_energy_confs)

    # Run Annealer
    cg_annealer.run_schedule()
    min_pose = cg_annealer._cg_mc_sim.get_minimum_energy_pose()
    print("Writing structure to:")
    print(job.fn("minimum.pdb"))
    min_pose.dump_pdb(job.fn("minimum.pdb"))
    
@FlowProject.operation
@FlowProject.pre.isfile("minimum.pdb")
def low_temperature_mc(job):
    os.chdir(job.ws)
    cg_pyrosetta.pyrosetta.init(
                      "--add_atom_types fa_standard parameters/atom_properties.txt " +
                      "--add_mm_atom_type_set_parameters fa_standard parameters/mm_atom_type_sets/mm_atom_properties.txt " +
                      "--extra_mm_params_dir parameters/mm_atom_type_sets"
                      )
    pose = cg_pyrosetta.pyrosetta.pose_from_file(job.fn("minimum.pdb"))

    # temporary output
    print(pose)

@FlowProject.operation
@FlowProject.pre.isfile("minimum.pdb")
@FlowProject.post(lambda job: 'minimum_energy' in job.document)
def store_minimum_energy(job):
    
    os.chdir(job.ws)
    cg_pyrosetta.pyrosetta.init(
                      "--add_atom_types fa_standard parameters/atom_properties.txt " +
                      "--add_mm_atom_type_set_parameters fa_standard parameters/mm_atom_type_sets/mm_atom_properties.txt " +
                      "--extra_mm_params_dir parameters/mm_atom_type_sets"
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

    pose = cg_pyrosetta.pyrosetta.pose_from_file("minimum.pdb")


    # Pose to be folded
    pose = cg_pyrosetta.pyrosetta.pose_from_sequence("X[CG11x3:CGLower]X[CG11x3]X[CG11x3]X[CG11x3]X[CG11x3:CGUpper]")
    change_lengths = cg_pyrosetta.CG_movers.setBondLengths(pose, {"BB1 BB2":job.sp.bb_length, "BB2 BB3":job.sp.bb_length, "BB3 BB1":job.sp.bb_length})
    change_lengths.apply(pose)

    # Build Minimizer
    mini = cg_pyrosetta.pyrosetta.rosetta.protocols.minimization_packing.MinMover()
    mini.min_type('lbfgs_armijo_nonmonotone')
    mini.score_function(energy_function)

    for _ in range(10):
        mini.apply(pose)

    with open("minimum.pdb", "r") as fn:
        for line in fn.readlines():
            print(line)
    
    min_energy = energy_function(pose)

    print(min_energy)

    pose.dump_pdb("new_minimum.pdb")

    with open("new_minimum.pdb", "r") as fn:
        for line in fn.readlines():
            print(line)

    job.document.minimum_energy = min_energy


# @FlowProject.operation
# def run_helical_analysis(job):
#     os.system("run stuff")


if __name__ == '__main__':
    FlowProject().main()