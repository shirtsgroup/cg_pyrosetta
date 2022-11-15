import signac
from signac import H5Store
import numpy as np
import sklearn
from sklearn.preprocessing import normalize
import os
import pickle
import csv

def main():
    project = signac.get_project()
    schema = project.detect_schema()
    min_samples, eps, sc_sizes = schema.items()
    sc_sizes = list(sc_sizes[1][float])
    sc_sizes.sort()
    summary_csv = []
    summary_csv.append(["sc_size", "eps", "min_samples", "n_clusters",
                        "avg_silhouette_score",
                        "min_z_score*", "minimum_energy",
                        "min_energy_cluster_rmsd",
                        "min_energy_min_intermedoid_rmsd",
                        "noise_cluster", "rmse_helix", "rmse_cyl",
                        "min_medoid"])
    for sc in sc_sizes:
        print("SC:", sc)
        sil_scores = []
        n_clusters = []
        identifiers = []
        for job in project.find_jobs({"sc_size":sc}):
            with H5Store(job.fn('clustering.h5')).open(mode='r') as data:
                avg_silhouette_score = data['silhouette_avg']
                cluster_sizes = data['cluster_sizes'][:]
            if avg_silhouette_score is None or len(cluster_sizes) < 4 or job.sp.eps == 0.05 or job.sp.eps < sc/2/10:
                pass
            else:
                sil_scores.append(avg_silhouette_score)
                n_clusters.append(len(cluster_sizes))
                identifiers.append({"sc_size":sc, "min_samples":job.sp.min_samples, "eps":job.sp.eps})
        if len(identifiers) == 0:
            print("No valid clustering")
            summary_csv.append([sc, "-", "-", "-", "-", "-", "-", "-", "-", "-", "-", "-", "-"])
            continue
        sil_metric = 1 - np.array(sil_scores)
        n_clusters = np.array(n_clusters)
        sil_metric = sil_metric/np.sqrt(np.dot(sil_metric, sil_metric))
        n_clust_metric = n_clusters/np.sqrt(np.dot(n_clusters, n_clusters))
        obj_func = sil_metric**2 + n_clust_metric**2
        target_ids = identifiers[np.argmin(obj_func)]
        for job in project.find_jobs(target_ids):
            print("min_samples:", job.sp.min_samples)
            print("eps:", job.sp.eps)
            
            # Load clustering info from clusternig.h5 file
            with H5Store(job.fn('clustering.h5')).open(mode='r') as data:
                cluster_sizes = data['cluster_sizes'][:]
                cluster_rmsd = data['cluster_rmsd'][:]
                z_score_matrix = data['z_score_matrix'][:]
                medoid_rmsd_matrix = data['medoid_rmsd_matrix'][:]
                mirror_cluster_rmsd = data['mirrored_minmium_rmsd_matrix'][:]
                avg_silhouette_score = data['silhouette_avg']

            
            # Load clusters energy values
            with open(job.fn('cluster_energies.pkl'), 'rb') as fp:
                all_cluster_energies = pickle.load(fp)
            
            # Get minimum energy cluster
            mean_cluster_energies = []
            for cluster_i in all_cluster_energies.keys():
                mean_cluster_energies.append(np.mean(all_cluster_energies[cluster_i]))
            mean_cluster_energies = np.array(mean_cluster_energies)
            min_energy = np.min(mean_cluster_energies)
            i_min = np.argmin(mean_cluster_energies)

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
            keys = list(all_cluster_energies.keys())
            min_energy_intermedoid_rmsds = medoid_rmsd_matrix[keys[i_min]]
            min_intermedoid_rmsd = np.min(min_energy_intermedoid_rmsds[np.invert(np.isclose(min_energy_intermedoid_rmsds, 0, atol=1e-3))])
            min_z_score = np.min(z_values[mirror_min_rmsd >= cluster_rmsd[keys[i_min]]])
            min_cluster_rmsd = cluster_rmsd[keys[i_min]]
            noise_cluster = False
            medoid_min = i_min

            if -1 in all_cluster_energies.keys():
                noise_cluster = True
                medoid_min -= 1
            min_medoid_str = "medoid_"+str(medoid_min)

            if medoid_min == -1:
                summary_csv.append([sc, job.sp.eps, job.sp.min_samples,
                                    len(cluster_sizes), avg_silhouette_score,
                                    min_z_score, min_energy, min_cluster_rmsd,
                                    min_intermedoid_rmsd, noise_cluster,
                                    "-", "-", min_medoid_str
                ])
                continue 
            
            # Get helix fitting data
            with open(job.fn('helix_fitting.pkl'), 'rb') as fp:
                helix_info_dict = pickle.load(fp)

            min_helix_rmse = helix_info_dict["cluster_"+str(medoid_min)+"_minimum"]["rmse_helix"]
            min_cyl_rmse = helix_info_dict["cluster_"+str(medoid_min)+"_minimum"]["rmse_cylinder"]
            
            # summary_csv.append(["sc_size", "eps", "min_samples", "n_clusters",
            #                     "avg_silhouette_score",
            #                     "min_z_score*", "minimum_energy",
            #                     "min_energy_cluster_rmsd",
            #                     "min_energy_avg_intermedoid_rmsd",
            #                     "noise_cluster", "rmse_helix", "rmse_cyl",
            #                     "min_medoid"])
            
            summary_csv.append([sc, job.sp.eps, job.sp.min_samples,
                                len(cluster_sizes), avg_silhouette_score,
                                min_z_score, min_energy, min_cluster_rmsd,
                                min_intermedoid_rmsd, noise_cluster,
                                min_helix_rmse, min_cyl_rmse, min_medoid_str
            ])

    with open("summary_stats.csv", "w") as f:
        writer = csv.writer(f)
        writer.writerows(summary_csv)



        
if __name__ == "__main__":
    main()