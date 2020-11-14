import cg_pyrosetta as cgpy
import os
import pickle

def main():

    cgpy.change_parameters.changeAtomParameters({"CG1" : ['X', 1.5, 0.2]},
                                    atom_types_path="parameters/atom_properties.txt",
                                    mm_atom_types_path = "parameters/mm_atom_type_sets/mm_atom_properties.txt")

    # Initialize CG_PyRosetta with extra mm_types
    cgpy.init()
    # cgpy.pyrosetta.init()
    # CG MC Annealer Parameters
    params = cgpy.CG_monte_carlo.\
        CGMonteCarloAnnealerParameters(n_inner=100,
                                    t_init=5,
                                    anneal_rate=0.9,
                                    n_anneals=3,
                                    annealer_criteron=cgpy.CG_monte_carlo.Repeat10Convergence,
                                    traj_out = "testing.pdb",
                                    mc_output = True,
                                    mc_traj = False,
                                    out_freq = 10,
                                    )
    # Energy Function
    energy_function = cgpy.CG_monte_carlo.\
        EnergyFunctionFactory().build_energy_function(
            {
                "mm_twist": 1,
                "mm_bend": 1,
                "fa_atr": 1,
                "fa_rep": 1,
                "fa_intra_rep": 1,
                "fa_intra_atr": 1
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
                                                "mini": 20
                                            }
                                            )

    # Initialize Annealer
    cg_annealer = cgpy.CG_monte_carlo.CGMonteCarloAnnealer(
        seq_mover=sequence_mover,
        score_function=energy_function,
        pose=pose,
        param_file_object=params
    )

    min_energy_confs = cgpy.CG_monte_carlo.MinEnergyConfigObserver(cg_annealer.get_mc_sim())

    cg_annealer.registerObserver(min_energy_confs)

    # Run Annealer
    cg_annealer.run_schedule()
    min_pose = cg_annealer._cg_mc_sim.get_minimum_energy_pose()
    cg_annealer._cg_mc_sim.pymol.apply(min_pose)

    
    # Extract all structures
    print("There are a total of", len(min_energy_confs.structures), "structures")
    print("Energies:")


    if not os.path.exists("structures"):
        os.mkdir("structures")

    for i in range(len(min_energy_confs.energies)):
        print("Structure:", i, "Energy:", min_energy_confs.energies[i])
        min_energy_confs.structures[i].dump_pdb("structures/structure_" + str(i) + ".pdb")

    # Extract bottom 20 structures
    sorted_structures = [pose for _, pose in sorted(zip(min_energy_confs.energies, min_energy_confs.structures))]
    if not os.path.exists("50_min_structures"):
        os.mkdir("50_min_structures")
    for i, p in enumerate(sorted_structures[:50]):
        p.dump_pdb("50_min_structures/structure_" + str(i) + ".pdb")

    # Save output config object
    with open("min_energy_configs.pkl", "wb") as f:
        pickle.dump(min_energy_confs, f, pickle.HIGHEST_PROTOCOL)




if __name__ == "__main__":
    main()
