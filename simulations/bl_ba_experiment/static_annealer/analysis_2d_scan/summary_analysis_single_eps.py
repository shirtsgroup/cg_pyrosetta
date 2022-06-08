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

def main():
    project = signac.get_project()
    schema = project.detect_schema()
    eps_values, bond_angles = schema.items()
    eps_values = list(eps_values[1][float])
    eps_values.sort()
    bond_angles = list(bond_angles[1][float])
    bond_angles.sort()
    for eps in eps_values:
        csv_table = []
        csv_table.append(["bond_angle", "eps", "n_clusters", "min_e_clusters_silhouette_score", "min_z_score*", "minimum_energy", "min_energy_cluster_rmsd", "min_energy_avg_intermedoid_rmsd", "noise_cluster", "rmse_helix", "rmse_cyl", "min_medoid"])
        for ba in bond_angles:
            csv_entry = []
            print("bond_angle:", ba)
            print("eps:", eps)
            for job in project.find_jobs({'bond_angle':ba, 'eps':eps}):
                with H5Store(job.fn('clustering.h5')).open(mode='r') as data:
                    cluster_sizes = data['cluster_sizes'][:]
                if len(cluster_sizes) > 2:
                    with H5Store(job.fn('clustering.h5')).open(mode='r') as data:
                        cluster_rmsd = data['cluster_rmsd'][:]
                        z_score_matrix = data['z_score_matrix'][:]
                        medoid_rmsd_matrix = data['medoid_rmsd_matrix'][:]
                        mirror_cluster_rmsd = data['mirrored_minmium_rmsd_matrix'][:]
                    csv_entry.append(ba)
                    csv_entry.append(eps)
                    csv_entry.append(len(cluster_sizes))
                    with open(job.fn('cluster_silhouette_scores.pkl'), 'rb') as fp:
                        cluster_silhouettes = pickle.load(fp)
                    cluster_silhouettes = np.array(cluster_silhouettes)
                    with open(job.fn('cluster_energies.pkl'), 'rb') as fp:
                        all_cluster_energies = pickle.load(fp)
                    
                    # Find minimum energy cluster
                    mean_cluster_energies = []
                    std_cluster_energies = []
                    keys = list(all_cluster_energies.keys())
                    for cluster_i in all_cluster_energies.keys():
                        mean_cluster_energies.append(np.mean(all_cluster_energies[cluster_i]))
                        std_cluster_energies.append(np.std(all_cluster_energies[cluster_i]))
                    mean_cluster_energies = np.array(mean_cluster_energies)
                    std_cluster_energies = np.array(std_cluster_energies)
                    min_energy = np.min(mean_cluster_energies)
                    i_min = np.argmin(mean_cluster_energies)

                    if -1 in all_cluster_energies.keys():
                        included_silhouettes = cluster_silhouettes[mean_cluster_energies[1:] < min_energy + std_cluster_energies[i_min-1]]
                        included_weights = cluster_sizes[mean_cluster_energies[1:] < min_energy + std_cluster_energies[i_min-1]]
                    else:
                        included_silhouettes = cluster_silhouettes[mean_cluster_energies < min_energy + std_cluster_energies[i_min]]
                        included_weights = cluster_sizes[mean_cluster_energies < min_energy + std_cluster_energies[i_min]]
                    csv_entry.append(np.average(included_silhouettes, weights = included_weights))

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
                    min_energy_intermedoid_rmsds = medoid_rmsd_matrix[keys[i_min]]
                    min_intermedoid_rmsd = np.min(min_energy_intermedoid_rmsds[np.invert(np.isclose(min_energy_intermedoid_rmsds, 0, atol=1e-3))])
                    print(z_values)
                    if len(z_values[mirror_min_rmsd >= cluster_rmsd[keys[i_min]]]) != 0:
                        csv_entry.append(np.min(z_values[mirror_min_rmsd >= cluster_rmsd[keys[i_min]]]))
                    else:
                        csv_entry.append("-")
                    csv_entry.append(min_energy)
                    csv_entry.append(cluster_rmsd[keys[i_min]])
                    min_energy_intermedoid_rmsds = medoid_rmsd_matrix[keys[i_min]]
                    csv_entry.append(min_intermedoid_rmsd)
                    noise_cluster = False
                    if -1 in all_cluster_energies.keys():
                        noise_cluster = True
                        medoid_min -= 1
                    csv_entry.append(noise_cluster)

                    # Get helix fitting data
                    with open(job.fn('helix_fitting.pkl'), 'rb') as fp:
                        helix_info_dict = pickle.load(fp)
                    
                    if medoid_min != -1:
                        csv_entry.append(helix_info_dict["medoid_"+str(medoid_min)]["rmse_helix"])
                        csv_entry.append(helix_info_dict["medoid_"+str(medoid_min)]["rmse_cylinder"])
                    else:
                        csv_entry.append("-")
                        csv_entry.append("-")

                    csv_entry.append("medoid_"+str(medoid_min))
                else:
                    csv_entry = [ba,"-","-","-","-","-","-","-","-","-","-","-"]
                csv_table.append(csv_entry)
        with open("summary_stats_" + str(eps) + ".csv", "w") as f:
            writer = csv.writer(f)
            writer.writerows(csv_table)

if __name__ == "__main__":
    main()

            
                


        
