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
    csv_table.append(["sc_size", "eps", "n_clusters", "avg_silhouette_score", "max_z_score", "minimum_energy", "min_energy_cluster_rmsd", "min_energy_avg_intermedoid_rmsd", "noise_cluster"])
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
        print(eps_values)
        print(silhouette_scores)
        print(n_clusters)
        silhoutte_scores = [x for _, x in sorted(zip(eps_values, silhouette_scores))]
        n_clusters = [x for _, x in sorted(zip(eps_values, n_clusters))]
        eps_values.sort()
        fig, ax1 = plt.subplots()
        ax2 = ax1.twinx()
        ax1.plot(eps_values, silhouette_scores)
        ax2.plot(eps_values, n_clusters, 'r')
        ax1.set_ylabel("Average Silhouette Score")
        ax2.set_ylabel("Number of Clusters", color='r')
        ax1.set_xlabel("EPS value")
        plt.savefig("eps_silhouette_plot_sc_" + str(sc) + ".jpg")
        
        max_silhouette = np.max(silhouette_scores)
        target_eps = eps_values[silhouette_scores.index(max_silhouette)]
        csv_entry = []
        print(sc)
        print(target_eps)
        for job in project.find_jobs({'sc_size':sc, 'eps':target_eps}):
            with H5Store(job.fn('clustering.h5')).open(mode='r') as data:
                cluster_sizes = data['cluster_sizes'][:]
                cluster_rmsd = data['cluster_rmsd'][:]
                avg_silhouette_score = data['silhouette_avg']
                z_score_matrix = data['z_score_matrix'][:]
                medoid_rmsd_matrix = data['medoid_rmsd_matrix'][:]
            csv_entry.append(sc)
            csv_entry.append(target_eps)
            csv_entry.append(len(cluster_sizes))
            csv_entry.append(avg_silhouette_score)

            with open(job.fn('cluster_energies.pkl'), 'rb') as fp:
                all_cluster_energies = pickle.load(fp)
            print(all_cluster_energies)
            mean_cluster_energies = []
            keys = list(all_cluster_energies.keys())
            for cluster_i in all_cluster_energies.keys():
                mean_cluster_energies.append(np.mean(all_cluster_energies[cluster_i]))
            min_energy = np.min(mean_cluster_energies)

            csv_entry.append(np.max(z_score_matrix[mean_cluster_energies.index(min_energy), :])) # Need to identify minium Z-score after degenerate distribution
            csv_entry.append(min_energy)
            csv_entry.append(cluster_rmsd[keys[mean_cluster_energies.index(min_energy)]])
            min_energy_intermedoid_rmsds = medoid_rmsd_matrix[keys[mean_cluster_energies.index(min_energy)]]
            print("Cluster:", min_energy_intermedoid_rmsds[min_energy_intermedoid_rmsds != 0])
            csv_entry.append(np.min(min_energy_intermedoid_rmsds[min_energy_intermedoid_rmsds != 0]))
            noise_cluster = False
            if -1 in all_cluster_energies.keys():
                noise_cluster = True
            csv_entry.append(noise_cluster)
            csv_table.append(csv_entry)


    with open("summary_stats.csv", "w") as f:
        writer = csv.writer(f)
        writer.writerows(csv_table)

    
if __name__ == "__main__":
    main()