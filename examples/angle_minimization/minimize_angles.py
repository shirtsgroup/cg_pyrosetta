# This is an executable script for debugging the angle
# minimization issue we've been dealing with


import cg_pyrosetta
import pyrosetta


# Build models
polyA = pyrosetta.pose_from_sequence('AAAAAAAA')
polyCG = pyrosetta.pose_from_sequence('X[CG31:CGLower]X[CG31]X[CG31]X[CG31]X[CG31:CGUpper]')

# Build PyMol visualizer
pymol = pyrosetta.PyMOLMover()
# pymol.apply(polyA)
pymol.apply(polyCG)

# Build Score Function
sf = pyrosetta.ScoreFunction()
sf.set_weight(pyrosetta.rosetta.core.scoring.mm_bend, 1)
sf.set_weight(pyrosetta.rosetta.core.scoring.mm_twist, 1)
sf.set_weight(pyrosetta.rosetta.core.scoring.fa_atr, 1)
sf.set_weight(pyrosetta.rosetta.core.scoring.fa_rep, 1)

# Build MoveMap
mm = pyrosetta.MoveMap()
mm.set(pyrosetta.rosetta.core.id.DOF_ID(pyrosetta.AtomID(3, 2),
                                        pyrosetta.rosetta.core.id.THETA),
       True)

# Build Minimizer
mini = pyrosetta.rosetta.protocols.minimization_packing.MinMover()
mini.score_function(sf)
mini.movemap(mm)

# Perturb angle
polyA.conformation().set_bond_angle(
    pyrosetta.AtomID(3, 2),
    pyrosetta.AtomID(2, 2),
    pyrosetta.AtomID(1, 2),
    0.5
)

polyCG.conformation().set_bond_angle(
    pyrosetta.AtomID(3, 2),
    pyrosetta.AtomID(2, 2),
    pyrosetta.AtomID(1, 2),
    0.5
)
# pymol.apply(polyA)
pymol.apply(polyCG)

for _ in range(20):
    mini.apply(polyA)
    pymol.apply(polyA)

for _ in range(20):
    mini.apply(polyCG)
    pymol.apply(polyCG)
