import cg_pyrosetta as cgpy

def main():
    import pdb
    pdb.set_trace()
    # Initialize CG_PyRosetta with extra mm_types
    cgpy.pyrosetta.init("--add_atom_types fa_standard atom_properties.txt --add_mm_atom_type_set_parameters fa_standard mm_atom_type_sets/mm_atom_properties.txt " +
                        "--extra_mm_params_dir mm_atom_type_sets")
    # cgpy.pyrosetta.init()
    # CG MC Annealer Parameters
    params = cgpy.CG_monte_carlo.\
        CGMonteCarloAnnealerParameters(n_inner=500,
                                    t_init=5,
                                    anneal_rate=0.9,
                                    n_anneals=20,
                                    annealer_criteron=cgpy.CG_monte_carlo.Repeat10Convergence,
                                    traj_out = "testing.pdb",
                                    mc_output = True,
                                    mc_traj = False,
                                    traj_freq = 500,
                                    )
    # Energy Function
    energy_function = cgpy.CG_monte_carlo.\
        EnergyFunctionFactory().build_energy_function(
            {
                "mm_twist": 1,
                "mm_bend": 1,
                "fa_atr": 1,
                "fa_rep": 1,
                "fa_intra_rep":1,
                "fa_intra_atr":1

            }
        )

    # Pose to be folded
    pose = cgpy.pyrosetta.pose_from_sequence("X[CG11x3:CGLower]X[CG11x3]X[CG11x3]X[CG11x3]X[CG11x3:CGUpper]")
    change_lengths = cgpy.CG_movers.setBondLengths(pose, {"BB1 BB2":1, "BB2 BB3":1, "BB3 BB1":1})
    change_lengths.apply(pose)
    print(energy_function(pose))
    # exit()
    # Build Minimizer

    mini = cgpy.pyrosetta.rosetta.protocols.minimization_packing.MinMover()
    mini.min_type('lbfgs_armijo_nonmonotone')
    mini.score_function(energy_function)

    # Build Cartesian Minimizer

    # cart_mini = cgpy.pyrosetta.rosetta.core.optimization.CartesianMinimizer()
    # cart_mmap = cgpy.pyrosetta.rosetta.core.optimization.CartesianMap()

    # MoveMap for Minimizer
    movemap = cgpy.pyrosetta.MoveMap()
    movemap.set_bb_true_range(1, pose.size())
    movemap.set_branches(True)
    movemap.set(cgpy.pyrosetta.rosetta.core.id.THETA, True)
    movemap.set(cgpy.pyrosetta.rosetta.core.id.PHI, True)

    # Insert MoveMap into Minimizer
    mini.movemap(movemap)

    # import pdb
    # pdb.set_trace()

    mini.apply(pose)


    # Build Sequence Mover for MC object
    sequence_mover_fct = cgpy.CG_monte_carlo.SequenceMoverFactory(pose, {"mini":mini})
    sequence_mover = sequence_mover_fct.build_seq_mover(
                                            {
                                                "small_dihe": 1,
                                                "small_angle": 1,
                                                # "sc_small_angle": 5,
                                                "mini": 50
                                            }
                                            )

    # Initialize Annealer
    cg_annealer = cgpy.CG_monte_carlo.CGMonteCarloAnnealer(
        seq_mover=sequence_mover,
        score_function=energy_function,
        pose=pose,
        param_file_object=params
    )

    # Run Annealer
    cg_annealer.run_schedule()
    min_pose = cg_annealer._cg_mc_sim.get_minimum_energy_pose()
    min_pose.dump_pdb("min.pdb")
    cg_annealer._cg_mc_sim.pymol.apply(min_pose)


if __name__ == "__main__":
    main()
