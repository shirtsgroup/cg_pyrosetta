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
import csv
from analyze_foldamers.parameters.helical_fitting_2 import * 
import mdtraj as md

# Script to get summary stats from each sidechain size
# I want to identify clusters

# SC Size,EPS,Number of Clusters, Average Silhouette Score, Z-score between lowest clusters (max Z-score?), minimum energy cluster average RMSD to medoid,

def main():
    project = signac.get_project()
    schema = project.detect_schema()
    sc_sizes, reps = schema.items()
    sc_sizes = list(sc_sizes[1][float])
    sc_sizes.sort()
    csv_table = []
    csv_table.append(["sc_size", "eps", "n_clusters", "avg_silhouette_score", "min_z_score*", "minimum_energy", "min_energy_cluster_rmsd", "min_energy_avg_intermedoid_rmsd", "noise_cluster", "rmse_helix", "rmse_cyl", "min_medoid"])
    for sc in sc_sizes:
        silhouette_scores = []
        eps_values = []
        n_clusters = []
        for job in project.find_jobs({'sc_size':sc}):
            # Populate eps_values with values with more than 1 or 2 clusters
            # choose the smallest one
            with H5Store(job.fn('clustering.h5')).open(mode='r') as data:
                cluster_sizes = data['cluster_sizes'][:]
                avg_silhouette_score = data['silhouette_avg']
            if len(cluster_sizes) > 2:
                silhouette_scores.append(avg_silhouette_score)
                eps_values.append(job.sp.eps)
                n_clusters.append(len(cluster_sizes))
        silhouette_scores = np.array([x for _, x in sorted(zip(eps_values, silhouette_scores))])
        n_clusters = np.array([x for _, x in sorted(zip(eps_values, n_clusters))])
        eps_values.sort()
        plot_n_clusters = False
        if plot_n_clusters:
            fig, ax1 = plt.subplots()
            ax2 = ax1.twinx()
            ax1.plot(eps_values, silhouette_scores)
            ax2.plot(eps_values, n_clusters, 'r')
            ax1.set_ylabel("Average Silhouette Score")
            ax2.set_ylabel("Number of Clusters", color='r')
            ax1.set_xlabel("EPS value")
            plt.savefig("eps_silhouette_plot_sc_" + str(sc) + ".jpg")
        
        print(silhouette_scores)
        if np.any(silhouette_scores > 0.6):
            small_n_silhouette_scores = silhouette_scores[silhouette_scores > 0.6]
            small_n_n_clusters = n_clusters[silhouette_scores > 0.6]
            target_n_clusters = np.min(small_n_n_clusters)
            target_eps = eps_values[n_clusters.tolist().index(target_n_clusters)]
        else:
            target_eps = eps_values[silhouette_scores.tolist().index(np.max(silhouette_scores))]

        csv_entry = []
        # Open specific job with max avg silhouette
        for job in project.find_jobs({'sc_size':sc, 'eps':target_eps}):
            # Load saved values from analysis
            with H5Store(job.fn('clustering.h5')).open(mode='r') as data:
                cluster_sizes = data['cluster_sizes'][:]
                cluster_rmsd = data['cluster_rmsd'][:]
                avg_silhouette_score = data['silhouette_avg']
                z_score_matrix = data['z_score_matrix'][:]
                medoid_rmsd_matrix = data['medoid_rmsd_matrix'][:]
                mirror_cluster_rmsd = data['mirrored_minmium_rmsd_matrix'][:]
            csv_entry.append(sc)
            csv_entry.append(target_eps)
            csv_entry.append(len(cluster_sizes))
            csv_entry.append(avg_silhouette_score)
            with open(job.fn('cluster_energies.pkl'), 'rb') as fp:
                all_cluster_energies = pickle.load(fp)
            
            # Find minimum energy cluster
            mean_cluster_energies = []
            keys = list(all_cluster_energies.keys())
            for cluster_i in all_cluster_energies.keys():
                mean_cluster_energies.append(np.mean(all_cluster_energies[cluster_i]))
            min_energy = np.min(mean_cluster_energies)
            i_min = mean_cluster_energies.index(min_energy)

            # Get Z values from minimum energy cluster
            z_values = z_score_matrix[i_min, :]
            medoid_min = i_min
            if -1 in all_cluster_energies.keys(): # get rid of noise cluster if present
                z_values = z_score_matrix[i_min, 1:]
                noise_cluster = True
                medoid_min -= 1

            # Save minimum Z value after removing mirrored clusters
            mirror_min_rmsd = mirror_cluster_rmsd[medoid_min, :]
            if len(mirror_min_rmsd)-1 != len(mirror_min_rmsd[mirror_min_rmsd > cluster_rmsd[medoid_min]]):
                print("Removed mirrored cluster!")
                print(medoid_rmsd_matrix[keys[i_min]])
                print(mirror_min_rmsd)
            min_energy_intermedoid_rmsds = medoid_rmsd_matrix[keys[i_min]]
            min_intermedoid_rmsd = np.min(min_energy_intermedoid_rmsds[np.invert(np.isclose(min_energy_intermedoid_rmsds, 0, atol=1e-3))])
            csv_entry.append(np.min(z_values[mirror_min_rmsd >= cluster_rmsd[keys[i_min]]]))
            csv_entry.append(min_energy)
            csv_entry.append(cluster_rmsd[keys[i_min]])
            min_energy_intermedoid_rmsds = medoid_rmsd_matrix[keys[i_min]]
            csv_entry.append(min_intermedoid_rmsd)
            noise_cluster = False
            medoid_min = i_min
            if -1 in all_cluster_energies.keys():
                noise_cluster = True
                medoid_min -= 1
            csv_entry.append(noise_cluster)

            # Get helix fitting data
            with open(job.fn('helix_fitting.pkl'), 'rb') as fp:
                helix_info_dict = pickle.load(fp)
            
            csv_entry.append(helix_info_dict["medoid_"+str(medoid_min)]["rmse_helix"])
            csv_entry.append(helix_info_dict["medoid_"+str(medoid_min)]["rmse_cylinder"])
            csv_entry.append("medoid_"+str(medoid_min))
            csv_table.append(csv_entry)

    with open("summary_stats.csv", "w") as f:
        writer = csv.writer(f)
        writer.writerows(csv_table)

    
if __name__ == "__main__":
    main()
