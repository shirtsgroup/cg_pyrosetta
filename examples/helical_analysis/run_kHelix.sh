#!/bin/bash

cat > input << EOF
inputhelix $1
helixout_name kHelix.out
coord_type 1
num_grid 360
natoms 15
nframes 1
grid_phi_beg 0
grid_phi_end 180
grid_theta_beg 0
grid_theta_end 180
helix_atom_names BB
print_to_plot 1
EOF
/home/tfobe/Research/Software/kHelios/bin/helios.o input
