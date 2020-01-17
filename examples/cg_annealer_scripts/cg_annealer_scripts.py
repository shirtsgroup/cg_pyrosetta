import cg_pyrosetta as cgpy

cgpy.pyrosetta.init("--add_mm_atom_type_set_parameters fa_standard mm_atom_type_sets/mm_atom_properties.txt " +
                    "--extra_mm_params_dir mm_atom_type_sets")

params = cgpy.CG_monte_carlo.\
    CGMonteCarloAnnealerParameters(n_inner=50000,
                                   t_init=100,
                                   anneal_rate=0.9,
                                   n_anneals=50,
                                   annealer_criteron=cgpy.CG_monte_carlo.Repeat10Convergence())

energy_function = cgpy.CG_monte_carlo.\
    EnergyFunctionFactory().build_energy_function(
        {
            "mm_twist": 1,
            "mm_bend": 1,
            "mm_lj_inter_rep": 1,
            "mm_lj_inter_atr": 1,
        }
    )

pose = cgpy.pyrosetta.pose_from_sequence("X[CG11x3:CGLower]X[CG11x3]X[CG11x3]X[CG11x3:CGUpper]")

sequence_mover_fct = cgpy.CG_monte_carlo.SequenceMoverFactory(pose)

sequence_mover = sequence_mover_fct.build_seq_mover(
                                           {
                                               "small_dihe": 10,
                                               "small_angle": 5,
                                               "sc_small_angle": 5
                                           }
                                           )
print(sequence_mover)

sequence_mover.apply(pose)

cg_annealer = cgpy.CG_monte_carlo.CGMonteCarloAnnealer(
    seq_mover=sequence_mover,
    score_function=energy_function,
    pose=pose,
    param_file_object=params
)

cg_annealer.run_schedule()
