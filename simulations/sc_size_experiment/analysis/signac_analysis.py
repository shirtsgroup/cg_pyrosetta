import signac
import numpy as np
import pandas as pd
import analyze_foldamers
import argparse
import flow
import matplotlib.pyplot as plt
from flow import FlowProject
import os


def parse_args():
    pass

@flow.directives(fork=True)
@FlowProject.operation
@FlowProject.post.isfile("energy_trajectory.pdf")
def get_energy_trajectory(job):
    old_project = signac.get_project(root=job.sp.production_path)

    # kT values for tempeature plot
    kts = np.array([10 * ( 0.9 ) ** i for i in range(50)])
    out_steps = np.array([10000 * i for i in range(50)])

    # Get trajectory and energies of each simulation
    traj_file_list = []
    all_energies = []
    energy_traj = []
    for old_job in old_project.find_jobs({"sc_size": job.sp.sc_size}):
        structure_file = old_job.fn("trajectory.pdb")
        traj_file_list.append(structure_file)
        energy_file = old_job.fn("energies.txt")
        energies = pd.read_csv(energy_file, header=None)
        energy_traj.append(energies.values[:, 1])
        all_energies.extend(energies.values[:, 1])
    all_energies = np.array(all_energies)
    
    # Energy trajectory plot
    fig, ax1 = plt.subplots(1,1,figsize = [20,10])
    ax2 = ax1.twinx()
    ax2.plot(out_steps, kts, 'r')
    ax2.set_xlim([0, 500000])
    ax2.set_ylabel("Simulated Temperature", color = 'r')
    for traj in energy_traj:
        ax1.plot(energies.values[:,0], traj, alpha = 0.4, lw=2)
        ax1.set_xlabel("Steps")
        ax1.set_ylabel("Energy (A.E.U)")
        ax1.set_title("SC size = " + str(round(job.sp['sc_size'],4)))

    fig.savefig(job.fn("energy_trajectory.jpg"), bbox_inches="tight")
    fig.savefig(job.fn("energy_trajectory.pdf"), bbox_inches="tight")

    plt.close("all")




if __name__ == "__main__":
    FlowProject().main()