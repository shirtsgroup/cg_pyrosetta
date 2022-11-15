import signac
from signac import H5Store
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
import shutil as sh
import os
import csv

def main():
    plot_intermediates = False
    project = signac.get_project()
    schema = project.detect_schema()
    eps, min_samples, bond_angles = schema.items()
    bond_angles = list(bond_angles[1][float])
    min_samples = list(min_samples[1][int])
    bond_angles.sort()
    min_samples.sort()
    colors = np.linspace(0, 1, len(min_samples))
    viridis = cm.get_cmap("viridis", len(min_samples))


    # Make directory to store plots + look up tables
    if not os.path.isdir("pareto_plots"):
        os.mkdir("pareto_plots")
    else:
        sh.rmtree("pareto_plots")
        os.mkdir("pareto_plots")

    for bond_angle in bond_angles:
        plt.figure(figsize = [5, 5], dpi = 300)
        plt.fill_between([0, 2], -1.2, 0.2, color = "red", alpha = 0.3)
        csv_table = []
        csv_table.append(["eps", "min_samples", "avg_silhouette", "n_clusters"])
        print("BA:", bond_angle)
        for min_sample in min_samples:
            for job in project.find_jobs({'bond_angle':bond_angle, 'min_samples':min_sample}):
                with H5Store(job.fn('clustering.h5')).open(mode='r') as data:
                    avg_silhouette_score = data['silhouette_avg']
                    cluster_sizes = data['cluster_sizes'][:]
                if avg_silhouette_score is None:
                    continue
                plt.scatter(len(cluster_sizes), -avg_silhouette_score, color = viridis(colors[min_samples.index(job.sp.min_samples)]), s = 5)
                csv_entry = [job.sp.eps, job.sp.min_samples, avg_silhouette_score, len(cluster_sizes)]
                csv_table.append(csv_entry)
            with open(os.path.join("pareto_plots", "ba_" + str(bond_angle) + "_stats.csv"), "w") as f:
                writer = csv.writer(f)
                writer.writerows(csv_table)
            plt.xlabel("N Clusters")
            plt.ylabel("Avg. Silhouette Score")
            plt.ylim([-1, 0])
            plt.title("Bond Angle: " + str(bond_angle))
            h = lambda c: plt.Line2D([],[],color=c, ls="",marker="o")
            plt.legend(handles=[h(viridis(i)) for i in colors],
                    labels=list(min_samples))
            if plot_intermediates:
                plt.savefig(os.path.join("pareto_plots", "ba_" + str(bond_angle) + "_ms_" + str(min_sample) + "_plot.png"))

        
        plt.xlabel("N Clusters")
        plt.ylabel("Avg. Silhouette Score")
        plt.ylim([-1, 0])
        plt.title("Bond Angle: " + str(bond_angle))
        h = lambda c: plt.Line2D([],[],color=c, ls="",marker="o")
        plt.legend(handles=[h(viridis(i)) for i in colors],
                labels=list(min_samples))        
        plt.savefig(os.path.join("pareto_plots", "ba_" + str(bond_angle) + "_plot.pdf"))
        plt.savefig(os.path.join("pareto_plots", "ba_" + str(bond_angle) + "_plot.png"))
        plt.close()

if __name__ == "__main__":
    main()