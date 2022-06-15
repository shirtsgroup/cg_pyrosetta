import cg_pyrosetta
import time
import numpy as np
from itertools import compress



def main():
    cg_pyrosetta.init(add_atom_types = "fa_standard atom_properties.txt", 
                      add_mm_atom_type_set_parameters = "fa_standard mm_atom_properties.txt",
                      extra_res_fa = "OCT_trappe.params"
    )
    pose = cg_pyrosetta.pyrosetta.pose_from_sequence("X[OCT]X[OCT]")
    pose.dump_pdb("OCT.pdb")

    cg_pyrosetta.pyrosetta

    # Score with LJ spheres
    energy_function = cg_pyrosetta.CG_monte_carlo.EnergyFunctionFactory().build_energy_function(
        {
            "fa_atr" : 1,
            "fa_rep" : 1,
            "fa_intra_rep" : 1,
            "fa_intra_atr" : 1,
            "fa_elec" : 1,
        }
    )

    print(energy_function(pose))
    pose.dump_pdb("OCT_scored.pdb")

    # Try out movers
    small = cg_pyrosetta.CG_movers.CGSmallMover(pose)
    print(len(small.torsions))

    pymol = cg_pyrosetta.pyrosetta.PyMOLMover()
    
    torsion_atom_types = []
    keep_torsions = []
    for i in range(len(small.torsions)):
        torsion = small.torsions[i]
        atom_types = [
            pose.residue(torsion[0].rsd()).atom_type(torsion[0].atomno()).atom_type_name(),
            pose.residue(torsion[1].rsd()).atom_type(torsion[1].atomno()).atom_type_name(),
            pose.residue(torsion[2].rsd()).atom_type(torsion[2].atomno()).atom_type_name(),
            pose.residue(torsion[3].rsd()).atom_type(torsion[3].atomno()).atom_type_name(),
            ]
        if "Ch1" in atom_types[1] and "Ch1" in atom_types[2]:
            print("Removing:", " ".join(atom_types))
            keep_torsions.append(False)
            continue
        if all([a == "C1" for a in atom_types]):
            print("Removing:", " ".join(atom_types))
            keep_torsions.append(False)
            continue
        if sum([a == "Ch1" for a in atom_types]) > 2:
            print("Removing:", " ".join(atom_types))
            keep_torsions.append(False)
            continue
        if sum([a == "Ch1" for a in atom_types[:2]]) == 2 or sum([a == "Ch1" for a in atom_types[2:]]) == 2:
            print("LOOK HERE! Removing:", " ".join(atom_types))
            keep_torsions.append(False)
            continue
        print("Moving torsion", " ".join(atom_types))
        torsion_atom_types.append(atom_types)
        pymol.apply(pose)
        pose.conformation().set_torsion_angle(*torsion, pose.conformation().torsion_angle(*torsion) + np.pi/6)
        pymol.apply(pose)
        pose.conformation().set_torsion_angle(*torsion, pose.conformation().torsion_angle(*torsion) - np.pi/6)
        pymol.apply(pose)
        keep_torsions.append(True)
    
    small.torsions = list(compress(small.torsions, keep_torsions))

    for _ in range(5000):
        small.apply(pose)
        pymol.apply(pose)
    
    energy_function(pose)

    

    

if __name__ == "__main__":
    main()