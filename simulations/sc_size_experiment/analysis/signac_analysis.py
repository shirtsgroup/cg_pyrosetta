import signac
import numpy as np
import pandas as pd
import analyze_foldamers
import argparse
import flow
import matplotlib.pyplot as plt
from flow import FlowProject
import os
import pickle
import parser
import yaml

@FlowProject.label
def energy_trajectory(job):
    return job.isfile("energy_trajectory.pdf")

# @FlowProject.label
# def clustering(job):
#     if len(job.stores['clustering']['labels']) > 0:
#         return True
#     else:
#         return False

@flow.directives(fork=False)
@FlowProject.operation
@FlowProject.post.isfile("energy_trajectory.pdf")
@FlowProject.post(energy_trajectory)
def get_energy_trajectory(job):
    param_dict = get_clustering_parameters(job)
    old_project = signac.get_project(root=param_dict["signac_dir"])

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
    
    with open(job.fn('traj_files.pkl'),'wb') as fp:
        pickle.dump(traj_file_list, fp)
    job.stores['traj']['energy_traj'] = np.array(energy_traj)
    job.stores['traj']['all_energies'] = all_energies
    

    plt.close("all")

def get_default_parameters():
    default_params = {
        "eps": 0.15,
        "min_samples": 30,
        "frame_start": 700,
        "frame_stride": 3,
        "frame_end": -1,
        "output_format": "pdb",
        "plot_silhouette": True,
        "filter": True,
        "filter_ratio": 0.5,
        "output_cluster_traj": True,
        "core_points_only": False,
        "parallel" : -1,
    }
    return default_params

def get_clustering_parameters(job):
    if os.path.isfile(job.fn("cluster_parameters.yml")):
        param_dict = read_yaml_params(job.fn("cluster_parameters.yml"))
    else:
        if os.path.isfile("cluster_parameters.yml"):
            param_dict = read_yaml_params("cluster_parameters.yml")
        else:
            param_dict = get_default_parameters()
    return param_dict


def read_yaml_params(filename):
    stream = open(filename, "r")
    params = yaml.safe_load(stream)
    return params

@flow.directives(fork=False)
@FlowProject.operation
@FlowProject.pre(energy_trajectory)
@FlowProject.post.isfile("clustering.h5")
def run_clustering(job):

    # Run clustering and store output in job.data

    param_dict = get_clustering_parameters(job)
    param_dict.pop("signac_dir", None)

    with open(job.fn('traj_files.pkl'), 'rb') as fp:
        traj_file_list = pickle.load(fp)
    (
        medoid_positions,
        cluster_sizes,
        cluster_rmsd,
        n_noise,
        silhouette_avg,
        labels,
        original_indices,
    ) = analyze_foldamers.cluster.get_cluster_medoid_positions_DBSCAN(
        traj_file_list,
        None,
        output_dir=job.fn("cluster_output"),
        **param_dict
    )

    job.stores['clustering']['medoid_positions'] = medoid_positions
    job.stores['clustering']['cluster_sizes'] = cluster_sizes
    job.stores['clustering']['cluster_rmsd'] = cluster_rmsd
    job.stores['clustering']['n_noise'] = n_noise
    job.stores['clustering']['silhouette_avg'] = silhouette_avg
    job.stores['clustering']['labels'] = labels
    job.stores['clustering']['original_indices'] = original_indices


@flow.directives(fork=False)
@FlowProject.operation
@FlowProject.pre.after(run_clustering)
@FlowProject.post.isfile("cluster_energy_dist_cluster_0.pdf")
def cluster_energy_distributions(job):
    if len(job.stores['clustering']['cluster_sizes']) > 0:
        # Output cluster energy distributions
        clusters = list(np.unique(job.stores['clustering']['labels']))
        all_cluster_energies = []
        if -1 in clusters:
            clusters.remove(-1)
        for i in clusters:
            plt.figure(figsize=[5, 5], dpi=100)
            cluster_indices = np.where(job.stores['clustering']['labels'] == i)[0]
            cluster_energies = all_energies[job.stores['clustering']['original_indices'][cluster_indices]]
            all_cluster_energies.append(cluster_energies)
            mean_energy = np.mean(cluster_energies)
            std_energy = np.std(cluster_energies)
            color = cm.nipy_spectral(float(i) / len(clusters))
            plt.hist(cluster_energies, bins=50, color="red")
            plt.xlabel("Energy (A.E.U)")
            plt.ylabel("Counts")
            plt.title("Cluster " + str(i) + " Energy Distribution")
            plt.savefig(
                os.path.join(
                    sub_dir, "cluster_energy_dist_cluster_" + str(i) + ".pdf"
                ), bbox_inches="tight"
            )
            plt.savefig(
                os.path.join(
                    sub_dir, "cluster_energy_dist_cluster_" + str(i) + ".jpg"
                ), bbox_inches="tight"
            )
            plt.close("all")



    




if __name__ == "__main__":
    FlowProject().main()