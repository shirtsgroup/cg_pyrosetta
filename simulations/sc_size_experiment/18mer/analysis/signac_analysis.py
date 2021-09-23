import signac
from signac import H5Store
import numpy as np
import pandas as pd
import analyze_foldamers
import argparse
import flow
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib import ticker
from flow import FlowProject
import os
import pickle
import parser
import yaml
import mdtraj as md

# Helper functions

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

def compute_2_population_z_score(X1, X2):
    mu_1 = np.mean(X1)
    mu_2 = np.mean(X2)
    sigma_1 = np.std(X1)
    sigma_2 = np.std(X2)

    z_score = (mu_1 - mu_2) / np.sqrt(np.square(sigma_1) + np.square(sigma_2)) 
    return z_score

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

# Flow Operations

@FlowProject.label
def energy_trajectory(job):
    return job.isfile("energy_trajectory.pdf")

@flow.directives(fork=False)
@FlowProject.operation
@FlowProject.post(lambda job: matching_paramfile_statepoint(job))
def set_eps_parameter(job):
    param_dict = get_clustering_parameters(job)
    param_dict["eps"] = job.sp.eps
    stream =  open(job.fn("cluster_parameters.yml"), "w")
    yaml.dump(param_dict, stream)

def matching_paramfile_statepoint(job):
    param_dict = get_clustering_parameters(job)
    return job.sp.eps == param_dict["eps"] 

# @FlowProject.label
# def clustering(job):
#     if len(job.stores['clustering']['labels']) > 0:
#         return True
#     else:
#         return False

@flow.directives(fork=False)
@FlowProject.operation
@FlowProject.post.isfile("energy_trajectory.pdf")
@FlowProject.post.isfile("traj_files.pkl")
@FlowProject.post.isfile("traj.h5")
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

@flow.directives(fork=False)
@FlowProject.operation
@FlowProject.pre(energy_trajectory)
@FlowProject.pre(lambda job: matching_paramfile_statepoint(job))
@FlowProject.post.isfile("clustering.h5")
def run_clustering(job):

    # Run clustering and store output in job.data

    param_dict = get_clustering_parameters(job)
    param_dict.pop("signac_dir", None)
    print(param_dict)

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

    # job.stores['clustering']['medoid_positions'] = medoid_positions datatype not supported
    job.stores['clustering']['cluster_sizes'] = np.array(cluster_sizes)
    job.stores['clustering']['cluster_rmsd'] = cluster_rmsd
    job.stores['clustering']['n_noise'] = n_noise
    job.stores['clustering']['silhouette_avg'] = silhouette_avg
    job.stores['clustering']['labels'] = labels
    job.stores['clustering']['original_indices'] = original_indices


@flow.directives(fork=False)
@FlowProject.operation
@FlowProject.pre.after(run_clustering)
@FlowProject.post.isfile("cluster_energy_dist_cluster_0.pdf")
@FlowProject.post.isfile("cluster_energies.pkl")
def cluster_energy_distributions(job):
    with H5Store(job.fn('traj.h5')).open(mode='r') as energy_data:
        with H5Store(job.fn('clustering.h5')).open(mode='r') as data:
            all_energies = energy_data['all_energies'][:]
            labels = data['labels'][:]
            cluster_sizes = data['cluster_sizes'][:]
            original_indices = data['original_indices'][:]
    if len(cluster_sizes) > 0:
        # Output cluster energy distributions
        clusters = list(np.unique(labels))
        all_cluster_energies = {}
        # if -1 in clusters:
        #    clusters.remove(-1)
        for i in clusters:
            plt.figure(figsize=[5, 5], dpi=100)
            cluster_indices = np.where(labels == i)[0]
            cluster_energies = all_energies[original_indices[cluster_indices]]
            all_cluster_energies[i] = cluster_energies
            mean_energy = np.mean(cluster_energies)
            std_energy = np.std(cluster_energies)
            color = cm.nipy_spectral(float(i) / len(clusters))
            plt.hist(cluster_energies, bins=50, color="red")
            plt.xlabel("Energy (A.E.U)")
            plt.ylabel("Counts")
            plt.title("Cluster " + str(i) + " Energy Distribution")
            plt.savefig(
                job.fn("cluster_energy_dist_cluster_" + str(i) + ".pdf"
                ), bbox_inches="tight"
            )
            plt.savefig(
                job.fn("cluster_energy_dist_cluster_" + str(i) + ".jpg"
                ), bbox_inches="tight"
            )
            plt.close("all")
        
        with open(job.fn('cluster_energies.pkl'),'wb') as fp:
            pickle.dump(all_cluster_energies, fp)
        


@flow.directives(fork=False)
@FlowProject.pre.after(run_clustering)
@FlowProject.pre.after(get_energy_trajectory)
@FlowProject.post.isfile("rmsd_energy_scatter_plot.pdf")
@FlowProject.operation
def cluster_energy_rmsd_scatter(job):
    with H5Store(job.fn('traj.h5')).open(mode='r') as energy_data:
        with H5Store(job.fn('clustering.h5')).open(mode='r') as data:
            all_energies = energy_data['all_energies'][:]
            labels = data['labels'][:]
            cluster_sizes = data['cluster_sizes'][:]
            original_indices = data['original_indices'][:]
            cluster_rmsd = data['cluster_rmsd'][:]
    param_dict = get_clustering_parameters(job)
    old_project = signac.get_project(root=param_dict["signac_dir"])
    # Getting a trajectory with all frames
    sim_traj = []
    all_traj = None
    for old_job in old_project.find_jobs({'sc_size':job.sp.sc_size}):
        job_traj = md.load(old_job.fn("trajectory.pdb"))
        sim_traj.append(job_traj)
        if all_traj is None:
            all_traj = job_traj
        else:
            all_traj = all_traj.join(job_traj)
    clusters = list(np.unique(labels))
    
    with open(job.fn('cluster_energies.pkl'), 'rb') as fp:
        all_cluster_energies = pickle.load(fp)    
    
    if -1 in clusters:
       clusters.remove(-1)

    plt.figure(figsize = [15, 10], dpi=500)
    energy_stdevs = []
    print(clusters)
    for i in clusters:
        i_str = str(i)
        medoid_mdtraj = md.load(job.fn("cluster_output/medoid_" + i_str + ".pdb"))
        cluster_mdtraj = md.load(job.fn("cluster_output/cluster_"+ i_str + ".pdb"))
        cluster_indices = np.where(labels == i)[0]
        print("Indices Again!")
        print(cluster_indices)
        print("original_indices[cluster_indices]")
        print(original_indices[cluster_indices])
        rmsds_medoid_i = md.rmsd(cluster_mdtraj, medoid_mdtraj)
        print(rmsds_medoid_i)
        print("Calculated:", np.mean(rmsds_medoid_i))
        print("Expected RMSD:", cluster_rmsd[i])

        test_slice = all_traj[original_indices[cluster_indices]]
        test_slice.save_pdb(job.fn("test_cluster_"+i_str+".pdb"))
        cluster_energies = all_energies[original_indices[cluster_indices]]
        color = cm.nipy_spectral(float(i) / len(clusters))
        plt.scatter(rmsds_medoid_i, cluster_energies, s=10, alpha=0.7, c=color)
        plt.xlabel("RMSD to Cluster Medoid (A.L.U)")
        plt.ylabel("Cluster Energies (A.E.U)")
        plt.xticks(fontsize=16)
        plt.yticks(fontsize=16)
    plt.legend(["Cluster "+ str(a) for a in clusters], loc = "best", prop={'size': 16})
    plt.savefig(job.fn("rmsd_energy_scatter_plot.pdf"), bbox_inches="tight")
    plt.savefig(job.fn("rmsd_energy_scatter_plot.jpg"), bbox_inches="tight")

    plt.close("all")


@flow.directives(fork=False)
@FlowProject.pre.after(run_clustering)
@FlowProject.pre.after(get_energy_trajectory)
@FlowProject.post.isfile("intermedoid_rmsd_matrix.pdf")
@FlowProject.operation
def inter_medoid_matrix(job):
    # Distance matrix of each medoid structure to one another
    with H5Store(job.fn('traj.h5')).open(mode='r') as energy_data:
        with H5Store(job.fn('clustering.h5')).open(mode='r') as data:
            all_energies = energy_data['all_energies'][:]
            labels = data['labels'][:]
            cluster_sizes = data['cluster_sizes'][:]
            original_indices = data['original_indices'][:]

    clusters = list(np.unique(labels))
    if -1 in clusters:
        clusters.remove(-1)
    
    medoid_rmsd_matrix = 1

    cluster_output_files = os.listdir(job.fn("cluster_output"))
    cluster_output_files.sort()
    medoid_traj = None
    for medoid_file in cluster_output_files:
        if "medoid" in medoid_file:
            if medoid_traj is None:
                medoid_traj = md.load(job.fn("cluster_output/"+medoid_file))
            else:
                medoid_traj = medoid_traj.join(md.load(job.fn("cluster_output/"+medoid_file)))

    medoid_traj.superpose(medoid_traj)
    medoid_rmsd_matrix = np.zeros((medoid_traj.n_frames, medoid_traj.n_frames))

    for i in range(medoid_traj.n_frames):
        medoid_rmsd_matrix[i, :] = md.rmsd(medoid_traj, medoid_traj[i])

    job.stores['clustering']['medoid_rmsd_matrix'] = medoid_rmsd_matrix

    plt.figure(figsize = [10,10])
    plt.matshow(medoid_rmsd_matrix, cmap = "viridis")
    plt.colorbar()
    ax = plt.gca()
    ax.set_title("Medoid RMSD Matrix")


    labels = ["Medoid " + str(i) for i in clusters]
    labels.insert(0, "")
    
    ax.set_xticklabels(labels)
    ax.set_yticklabels(labels)
    
    ax.xaxis.set_major_locator(ticker.MultipleLocator(1))
    ax.yaxis.set_major_locator(ticker.MultipleLocator(1))

    plt.xticks(rotation=45)

    plt.savefig(job.fn("intermedoid_rmsd_matrix.pdf"), bbox_inches="tight")
    plt.savefig(job.fn("intermedoid_rmsd_matrix.pdf"), bbox_inches="tight")
    plt.close("all")

@flow.directives(fork=False)
@FlowProject.pre.after(cluster_energy_distributions)
@FlowProject.pre.after(get_energy_trajectory)
@FlowProject.post.isfile("energy_zscore_matrix.pdf")
@FlowProject.operation
def energy_z_score_matrix(job):
    # Distance matrix of each medoid structure to one another
    with H5Store(job.fn('traj.h5')).open(mode='r') as energy_data:
        with H5Store(job.fn('clustering.h5')).open(mode='r') as data:
            all_energies = energy_data['all_energies'][:]
            labels = data['labels'][:]
            cluster_sizes = data['cluster_sizes'][:]
            original_indices = data['original_indices'][:]

    with open(job.fn('cluster_energies.pkl'), 'rb') as fp:
        all_cluster_energies = pickle.load(fp)

    clusters = list(np.unique(labels))
    # if -1 in clusters:
    #    clusters.remove(-1)
    

    z_score_matrix = np.zeros((len(clusters), len(clusters)))
    for i in range(len(clusters)):
        for j in range(len(clusters)):
            if i != j:
                z_score_matrix[i,j] = np.abs(compute_2_population_z_score(all_cluster_energies[clusters[i]], all_cluster_energies[clusters[j]]))

    job.stores['clustering']['z_score_matrix'] = z_score_matrix
    
    plt.figure(figsize = [10,10])
    plt.matshow(z_score_matrix, cmap = "viridis")
    plt.colorbar()
    ax = plt.gca()
    ax.set_title("Energy Z-score Matrix")


    if -1 in clusters:
        labels = ["Medoid " + str(i) for i in range(0,len(clusters)-1)]
        labels.insert(0, "Noise")
    else:
        labels = ["Medoid " + str(i) for i in range(len(clusters))]

    labels.insert(0, "")

    ax.xaxis.set_major_locator(ticker.MultipleLocator(1))
    ax.yaxis.set_major_locator(ticker.MultipleLocator(1))
    
    ax.set_xticklabels(labels)
    ax.set_yticklabels(labels)

    ax.xaxis.set_major_locator(ticker.MultipleLocator(1))
    ax.yaxis.set_major_locator(ticker.MultipleLocator(1))

    plt.xticks(rotation=90)

    plt.savefig(job.fn("energy_zscore_matrix.pdf"), bbox_inches="tight")
    plt.savefig(job.fn("energy_zscore_matrix.jpg"), bbox_inches="tight")
    plt.close("all")




if __name__ == "__main__":
    FlowProject().main()