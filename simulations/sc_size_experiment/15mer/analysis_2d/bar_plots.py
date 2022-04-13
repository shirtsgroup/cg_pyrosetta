import numpy as np
import matplotlib.pyplot as plt
import matplotlib.pylab as pl
import csv
import seaborn as sns # improves plot aesthetics
import os
import shutil as sh

def main():
    with open("summary_stats.csv", newline="") as csvfile:
        reader = csv.reader(csvfile)
        all_data = list(reader)

    # Extract metrics from list of lists
    all_data.pop(0)
    labels = ["$SC$ = " + str(round(float(a[0]),3)) + "$\sigma_{B}$" for a in all_data][:11]
    avg_sil_score = np.array([float(a[4]) for a in all_data])[:11]
    min_z_scores = np.array([float(a[5]) for a in all_data])[:11]
    min_cluster_rmsd = np.array([float(a[7]) for a in all_data])[:11]
    min_cluster_rmsd = - min_cluster_rmsd[:11]
    inter_min_rmsd = np.array([float(a[8]) for a in all_data])[:11]

    # Normalize all metrics to be between 0 and 1
    norm_sil_scores = (avg_sil_score - np.min(avg_sil_score)) / (np.max(avg_sil_score) - np.min(avg_sil_score))
    norm_z_scores = (min_z_scores - np.min(min_z_scores)) / (np.max(min_z_scores) - np.min(min_z_scores))
    norm_cluster_rmsd = (min_cluster_rmsd - np.min(min_cluster_rmsd)) / (np.max(min_cluster_rmsd) - np.min(min_cluster_rmsd))
    norm_inter_rmsd = (inter_min_rmsd - np.min(inter_min_rmsd)) / (np.max(inter_min_rmsd) - np.min(inter_min_rmsd))
    
    # Plot normalized data set
    fig = plt.figure(figsize=[9, 3], dpi=300)

    ax = fig.add_subplot()
    x = np.arange(len(labels))
    width = 0.2
    rects1 = ax.bar(x - 1.5*width, norm_sil_scores, width, label='SS')
    rects2 = ax.bar(x - 0.5*width, norm_z_scores, width, label='EZG')
    rects3 = ax.bar(x + 0.5*width, norm_cluster_rmsd, width, label='CRMSD')
    rects4 = ax.bar(x + 1.5*width, norm_inter_rmsd, width, label='MIMR')
    
    # Plot text
    ax.set_ylabel('Normalized Metrics')
    ax.set_title('$SC$ Parameter Scan Metric Comparison')
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=45, ha="right")
    ax.legend(bbox_to_anchor=(1, 1))

    plt.grid(visible=True, axis="y")
    ax.set_axisbelow(True)
    plt.tight_layout()


    plt.savefig("bar_plot_test.png")
    
if __name__ == "__main__":
    main()