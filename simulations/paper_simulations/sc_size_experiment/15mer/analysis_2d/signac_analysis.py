import signac
from signac import H5Store
import numpy as np
import pandas as pd
import analyze_foldamers
import argparse
import flow
import matplotlib.pyplot as plt
from matplotlib import cm, ticker, rc
from flow import FlowProject
from sklearn.metrics import silhouette_samples
import os
import copy
import pickle
import parser
import yaml
import shutil as sh
import mdtraj as md

# rc('text', usetex=True)

# Helper functions

def get_default_parameters():
    default_params = {
        "eps": 0.15,
        "min_samples": 30,
        "frame_start": 1400,
        "frame_stride": 6,
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

@FlowProject.label
def has_clusters(job):
    return job.isfile("cluster_output/medoid_0.pdb")

@flow.directives(fork=False)
@FlowProject.operation
@FlowProject.post(lambda job: matching_paramfile_statepoint(job))
def set_cluster_parameter(job):
    param_dict = get_clustering_parameters(job)
    param_dict["eps"] = job.sp.eps
    param_dict["min_samples"] = job.sp.min_samples
    stream =  open(job.fn("cluster_parameters.yml"), "w")
    yaml.dump(param_dict, stream)

def matching_paramfile_statepoint(job):
    param_dict = get_clustering_parameters(job)
    return job.sp.eps == param_dict["eps"] and job.sp.min_samples == param_dict["min_samples"]

# @FlowProject.label
# def clustering(job):
#     if len(job.stores['clustering']['labels']) > 0:
#         return True
#     else:
#         return False

# @flow.directives(fork=False)
# @FlowProject.operation
# @FlowProject.post.isfile("acc_ratio_traj.pdf")
# def get_acc_ratio_traj(job):
#     print("Working on job ID:", job)
#     param_dict = get_clustering_parameters(job)
#     old_project = signac.get_project(root=param_dict["signac_dir"])
#     plt.figure()
#     for old_job in old_project.find_jobs({"sc_size": job.sp.sc_size}):
#         acc_ratio_traj = np.genfromtxt(old_job.fn("acc_ratio.txt"), delimiter = ",")
#         plt.plot(acc_ratio_traj[:, 0], acc_ratio_traj[:, 1])
#     plt.xlabel("Steps")
#     plt.ylabel("Acceptance Ratio (over 10000 steps)")
# 
#     plt.savefig(job.fn("acc_ratio_traj.pdf"), dpi=500)
#     plt.savefig(job.fn("acc_ratio_traj.png"), dpi=500)


@flow.directives(fork=False)
@FlowProject.operation
@FlowProject.post.isfile("energy_trajectory_10_10.pdf")
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
    
    aspect_ratios = [[10,10], [15,10]]
    for ar in aspect_ratios:
        # Energy trajectory plot
        fig, ax1 = plt.subplots(1,1,figsize = ar, dpi=600)
        ax2 = ax1.twinx()
        ax2.plot(out_steps, kts/0.2, 'r')
        # ax2.set_xlim([0, 500000])
        ax2.set_ylabel(r'Simulated Temperature ($\epsilon_{B}$)', color = 'r', fontsize=20)
        plt.rcParams.update({'font.size' : 18})
        for traj in energy_traj:
            ax1.plot(energies.values[:,0], traj/0.2, alpha = 0.4, lw=2)
            ax1.set_xlabel("Steps", fontsize=20)
            ax1.set_ylabel(r'Energy ($\epsilon_{B}$)', fontsize=20)
            ax1.set_title(r'$SC$ = ' + str(round(job.sp['sc_size'],4)) + r'$R^{min}_{B}$', fontsize=20)
            ax1.tick_params(axis = "both", labelsize=16)
            ax2.tick_params(axis = "both", labelsize=16)
        fig.savefig(job.fn("energy_trajectory_" + "_".join([str(i) for i in ar]) + ".jpg"), bbox_inches="tight")
        fig.savefig(job.fn("energy_trajectory_" + "_".join([str(i) for i in ar]) + ".pdf"), bbox_inches="tight")
        
    with open(job.fn('traj_files.pkl'),'wb') as fp:
        pickle.dump(traj_file_list, fp)
    job.stores['traj']['energy_traj'] = np.array(energy_traj)
    job.stores['traj']['all_energies'] = all_energies
    
    plt.close("all")

    for a in energy_traj:
        print("Sim. steps:", len(a))
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
@FlowProject.pre.isfile("clustering.h5")
@FlowProject.pre(has_clusters)
@FlowProject.post.isfile("cluster_silhouette_scores.pkl")
def get_cluster_silhouette_scores(job):
    # Saved data values we need to calculate individual cluster SSs
    with H5Store(job.fn('clustering.h5')).open(mode='r') as data:
        labels = data['labels'][:]
        cluster_sizes = data['cluster_sizes'][:]
        silhouette_avg = data['silhouette_avg']

    # If there is only 1 cluster we can't calculate SSs
    if len(cluster_sizes) > 1:
        print("There are", len(np.unique(labels)), "clusters")
       
        # Load cluster files of all structures
        new_labels = []
        clusters_traj = None
        for i in np.unique(labels):
            cluster_file = job.fn("cluster_output/cluster_"+str(i)+".pdb")
            if i == -1:
                cluster_file = job.fn("cluster_output/cluster_noise.pdb")
            traj_i = md.load(cluster_file)
            if clusters_traj  is None:
                clusters_traj = traj_i
            else:
                clusters_traj = clusters_traj.join(traj_i)
                # Reconstruct labels to match current order of files
            for _ in range(traj_i.n_frames):
                new_labels.append(i)
        new_labels = np.array(new_labels)

        # Reconstruct distance matrix
        distance_matrix = np.zeros((len(new_labels), len(new_labels)))
        for i in range(1, clusters_traj.n_frames):
            md.Trajectory.superpose(clusters_traj[i], clusters_traj[0])
        
        for i in range(clusters_traj.n_frames):
            distance_matrix[i] = md.rmsd(clusters_traj, clusters_traj, i)

        silhouette_scores = silhouette_samples(distance_matrix, new_labels)
        avg_cluster_silhouette_scores = []
        for i in np.unique(labels):
            if i == -1:
                continue
            ss_i = silhouette_scores[np.where(new_labels == i)]
            avg_cluster_silhouette_scores.append(np.mean(ss_i))

        with open(job.fn('cluster_silhouette_scores.pkl'),'wb') as fp:
            pickle.dump(avg_cluster_silhouette_scores, fp)
    else:
        print("There is", len(np.unique(labels)), "cluster")
        with open(job.fn('cluster_silhouette_scores.pkl'),'wb') as fp:
            pickle.dump([0], fp)



@flow.directives(fork=False)
@FlowProject.operation
@FlowProject.pre.after(run_clustering)
@FlowProject.pre.after(get_energy_trajectory)
@FlowProject.post.isfile("cluster_energy_dist_cluster_0.pdf")
@FlowProject.post.isfile("cluster_energies.pkl")
@FlowProject.pre(has_clusters)
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
            plt.hist(cluster_energies/0.2, bins=50, color="red")
            plt.xlabel(r'Energy ($\epsilon_{B}$)')
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
@FlowProject.operation
@FlowProject.pre.after(cluster_energy_distributions)
@FlowProject.pre.after(get_energy_trajectory)
@FlowProject.post.isfile("cluster_output/cluster_0_minimum.pdb")
@FlowProject.pre(has_clusters)
def write_cluster_minimum_structures(job):
    print("Working in", job.fn(""))
    with open(job.fn('cluster_energies.pkl'), 'rb') as fp:
        all_cluster_energies = pickle.load(fp)

    for cluster_id in all_cluster_energies.keys():
        i = cluster_id
        if cluster_id == -1:
           cluster_id = "noise"
        cluster = md.load(job.fn("cluster_output/cluster_"+ str(cluster_id) + ".pdb"))
        i_min = np.argmin(all_cluster_energies[i])

        min_structure = cluster[i_min]
        min_structure.save_pdb(job.fn("cluster_output/cluster_"+ str(cluster_id) + "_minimum.pdb"))

@flow.directives(fork=False)
@FlowProject.pre.after(run_clustering)
@FlowProject.pre.after(get_energy_trajectory)
@FlowProject.pre.after(write_cluster_minimum_structures)
@FlowProject.post.isfile("rmsd_energy_scatter_plot.pdf")
@FlowProject.pre(has_clusters)
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
    for i in clusters:
        i_str = str(i)
        medoid_mdtraj = md.load(job.fn("cluster_output/cluster_" + i_str + "_minimum.pdb"))
        cluster_mdtraj = md.load(job.fn("cluster_output/cluster_"+ i_str + ".pdb"))
        cluster_indices = np.where(labels == i)[0]
        rmsds_medoid_i = md.rmsd(cluster_mdtraj, medoid_mdtraj)
        cluster_energies = all_energies[original_indices[cluster_indices]]
        # color = cm.tab10(float(i) / len(clusters))
        plt.scatter(rmsds_medoid_i*10, cluster_energies/0.2, s=10, alpha=0.7)
        plt.xlabel(r'RMSD to Cluster Minimum ($R^{min}_{B}$)', fontsize=20)
        plt.ylabel(r'Cluster Energies ($\epsilon_{B}$)', fontsize=20)
        plt.xticks(fontsize=16)
        plt.yticks(fontsize=16)
    plt.title("SC =" + str(round(job.sp.sc_size,3)) + r'$R^{min}_{B}$', fontsize=20)
    plt.legend(["Cluster "+ str(a) for a in clusters], loc = "lower right", fontsize=16)
    plt.savefig(job.fn("rmsd_energy_scatter_plot.pdf"), bbox_inches="tight")
    plt.savefig(job.fn("rmsd_energy_scatter_plot.jpg"), bbox_inches="tight")

    plt.close("all")

@flow.directives(fork=False)
@FlowProject.pre.after(run_clustering)
@FlowProject.pre.after(get_energy_trajectory)
@FlowProject.pre.after(write_cluster_minimum_structures)
@FlowProject.post.isfile("intermedoid_rmsd_matrix.pdf")
@FlowProject.pre(has_clusters)
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

    medoid_traj = None
    for i in range(len(cluster_sizes)):
        if medoid_traj is None:
            medoid_traj = md.load(job.fn("cluster_output/cluster_"+str(i)+"_minimum.pdb"))
        else:
            medoid_traj = medoid_traj.join(md.load(job.fn("cluster_output/cluster_"+str(i)+"_minimum.pdb")))
    
    medoid_traj.superpose(medoid_traj)
    medoid_rmsd_matrix = np.zeros((medoid_traj.n_frames, medoid_traj.n_frames))
    medoid_rmsd_min_mirror_matrix = np.zeros((medoid_traj.n_frames, medoid_traj.n_frames))

    for i in range(medoid_traj.n_frames):
        medoid_rmsd_matrix[i, :] = md.rmsd(medoid_traj, medoid_traj[i])
        for j in range(medoid_traj.n_frames):
            mirror = copy.deepcopy(medoid_traj[j])
            mirror.xyz[0, :, 0] = -mirror.xyz[0, :, 0]
            medoid_rmsd_min_mirror_matrix[i, j] = np.min([
                md.rmsd(medoid_traj[j], medoid_traj[i]),
                md.rmsd(mirror, medoid_traj[i])
            ])


    job.stores['clustering']['medoid_rmsd_matrix'] = medoid_rmsd_matrix
    job.stores['clustering']['mirrored_minmium_rmsd_matrix'] = medoid_rmsd_min_mirror_matrix

    plt.figure(figsize = [10,10])
    plt.matshow(medoid_rmsd_matrix*10, cmap = "viridis")
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
@FlowProject.pre(has_clusters)
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
        labels = ["Medoid " + str(i) for i in range(len(clusters)-1)]
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
    plt.savefig(job.fn("energy_zscore_matrix.png"), bbox_inches="tight")
    plt.close("all")


@flow.directives(fork=False)
@FlowProject.pre.after(run_clustering)
@FlowProject.post.isfile("helix_fitting.pkl")
@FlowProject.pre(has_clusters)
@FlowProject.operation
def fit_medoids_to_helices(job):
    if os.path.exists(job.fn("helix_fitting")):
        sh.rmtree(job.fn("helix_fitting"))
    os.mkdir(job.fn("helix_fitting"))
    medoid_files = [f for f in os.listdir(job.fn("cluster_output")) if "minimum" in f]
    helix_fitting_info = {}
    for medoid in medoid_files:
        print("Working on", medoid)
        medoid_id = medoid.split(".")[0]
        helix_fitting_info[medoid_id] = {}
        structure = md.load(job.fn(os.path.join("cluster_output",medoid)))
        top = structure.topology
        bb_helix = structure.atom_slice(top.select("name BB1 BB2 BB3"))
        # Scale coordinates
        bb_helix_whole = 100*bb_helix.xyz[0]
        bb_helix = 100*bb_helix.xyz[0][2:-2] # fit just the internal residues
        n_residues = bb_helix.shape[0]
        print(bb_helix.shape)
        print(bb_helix)
        entries = []
        RMSEs = []        
        for i in range(20):
            x0 = np.array([1, 0, 0, 0, 0, 0, 0])
            radius, w, phi, z_tot, rotation, center, normal, sse_helix, sse_cylinder = analyze_foldamers.parameters.helical_fitting_2.fit_helix_to_points(bb_helix, x0)
            RMSE_tot = np.sqrt(sse_cylinder/bb_helix.shape[0]) + np.sqrt(sse_helix/bb_helix.shape[0])
            entries.append([radius, w, phi, z_tot, rotation, center, normal, sse_helix, sse_cylinder])
            RMSEs.append(RMSE_tot)        

        # Pick helix fit that gives the lowest cylinder and helix RMSE
        i_min = RMSEs.index(np.min(RMSEs))
        radius, w, phi, z_tot, rotation, center, normal, sse_helix, sse_cylinder = entries[i_min]
        
        # Print helix parameters
        print("Helix Prarameters")
        print("Radius:", radius/100)
        print("Omega:", w)
        print("Phi:", phi)
        print("Rotation Matrix:\n", rotation)
        print("Center:\n", center)
        print("Normal:\n", normal)
        print("Helix SSE:", sse_helix)
        print("Cylinder SSE:", sse_cylinder)


        # Calculate residues per turn
        total_rotation = z_tot*np.abs(w)/2/np.pi
        print("Total Rotation:", total_rotation)
        print("Residues per turn", (n_residues-1)/total_rotation)

        # Save fitting output to dictionary
        helix_fitting_info[medoid_id]["radius"] = radius/100
        helix_fitting_info[medoid_id]["w"] = w
        helix_fitting_info[medoid_id]["phi"] = phi
        helix_fitting_info[medoid_id]["z_tot"] = z_tot
        helix_fitting_info[medoid_id]["rotation"] = rotation
        helix_fitting_info[medoid_id]["center"] = center
        helix_fitting_info[medoid_id]["normal"] = normal
        helix_fitting_info[medoid_id]["residues_per_turn"] = (n_residues-1)/total_rotation
        helix_fitting_info[medoid_id]["rmse_cylinder"] = np.sqrt(sse_cylinder/bb_helix.shape[0])
        helix_fitting_info[medoid_id]["rmse_helix"] = np.sqrt(sse_helix/bb_helix.shape[0])
        
        # Plot figure
        fig = plt.figure()
        ax = fig.gca(projection='3d')
        ax.scatter3D(bb_helix_whole[:, 0], bb_helix_whole[:, 1], bb_helix_whole[:, 2], s=125)
        ax.plot3D(bb_helix_whole[:, 0], bb_helix_whole[:, 1], bb_helix_whole[:, 2])

        # Plot fitted helix
        t = np.linspace(0, z_tot, 100)
        fitted_helix = np.zeros([len(t), 3])
        fitted_helix[:, 0] = radius * np.cos(w*t + phi)
        fitted_helix[:, 1] = radius * np.sin(w*t + phi)
        fitted_helix[:, 2] = t
        fitted_helix = np.dot(fitted_helix, rotation.transpose())
        fitted_helix = fitted_helix -  np.mean(fitted_helix, axis = 0)
        ax.plot3D(fitted_helix[:, 0], fitted_helix[:, 1], fitted_helix[:, 2])

        # Show projection of helix on identified plane
        point = center
        ax.scatter3D(point[0], point[1], point[2], c="black")
        normal = normal / np.sqrt(np.dot(normal, normal))
        ax.plot3D(np.array([0, normal[0]])+point[0], np.array([0, normal[1]])+point[1], np.array([0, normal[2]])+point[2], "black")

        # Plot projection of helix onto plane
        centered = bb_helix - point
        dist = np.dot(centered, normal)
        projected_points = bb_helix - dist[:, np.newaxis] * normal[np.newaxis, :]
        ax.scatter3D(projected_points[:, 0], projected_points[:, 1], projected_points[:, 2], c="black")
        ax.plot3D(projected_points[:, 0], projected_points[:, 1], projected_points[:, 2], "black")
        plt.savefig(job.fn("helix_fitting/helix_fit_" + medoid_id +".pdf"), dpi=500)
        plt.savefig(job.fn("helix_fitting/helix_fit_" + medoid_id +".png"), dpi=500)
        
        plot_values = [bb_helix_whole[:, 0], bb_helix_whole[:, 1], bb_helix_whole[:, 2],
                       fitted_helix[:, 0], fitted_helix[:, 1], fitted_helix[:, 2],
                       point[0], point[1], point[2],
                       np.array([0, normal[0]])+point[0], np.array([0, normal[1]])+point[1], np.array([0, normal[2]])+point[2],
                       projected_points[:, 0], projected_points[:, 1], projected_points[:, 2]]
        pickle.dump(plot_values, open(job.fn("helix_fitting/helix_fit_" + medoid_id +".pkl"), 'wb'))

        #Save info dictionary for later use

    with open(job.fn('helix_fitting.pkl'),'wb') as fp:
        pickle.dump(helix_fitting_info, fp)

@flow.directives(fork=False)
@FlowProject.operation
@FlowProject.pre.after(cluster_energy_distributions)
@FlowProject.pre.after(get_energy_trajectory)
@FlowProject.pre(has_clusters)
@FlowProject.post.isfile("cluster_sampling.pkl")
def minimum_cluster_simulation_count(job):
    with H5Store(job.fn('clustering.h5')).open(mode='r') as data:
        original_indices = data['original_indices'][:]
        labels = data['labels'][:]
    with open(job.fn('traj_files.pkl'), 'rb') as fp:
        traj_file_list = pickle.load(fp)

    sim_delimiters = [0]
    for f in traj_file_list:
        traj = md.load(f)
        n_steps = len(traj)
        sim_delimiters.append(n_steps + sim_delimiters[-1])
    
    clusters = np.unique(labels)

    cluster_sampling_dict = {}
    for cluster_id in clusters:
        cluster_original_indices = original_indices[np.where(labels == cluster_id)[0]]
        sim_ids = np.digitize(cluster_original_indices, sim_delimiters)
        print("Cluster", cluster_id, "has structures from simulations:\n", np.unique(sim_ids))
        cluster_sampling_dict[cluster_id] = len(np.unique(sim_ids))
        
    with open(job.fn('cluster_sampling.pkl'),'wb') as fp:
        pickle.dump(cluster_sampling_dict, fp)


if __name__ == "__main__":
    FlowProject().main()
