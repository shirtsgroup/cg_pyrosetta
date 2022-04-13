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
    all_data = all_data[:27]
    all_data.pop(0)
    labels = [r'$\theta_{B} = $' + str(int(float(a[0]))) for a in all_data]
    avg_sil_score = np.array([float(a[4]) for a in all_data])
    min_z_scores = np.array([float(a[5]) for a in all_data])
    min_cluster_rmsd = np.array([float(a[7]) for a in all_data])
    min_cluster_rmsd = np.max(min_cluster_rmsd) - min_cluster_rmsd
    inter_min_rmsd = np.array([float(a[8]) for a in all_data])

    # Normalize all metrics to be between 0 and 1
    norm_sil_scores = (avg_sil_score - np.min(avg_sil_score)) / (np.max(avg_sil_score) - np.min(avg_sil_score))
    norm_z_scores = (min_z_scores - np.min(min_z_scores)) / (np.max(min_z_scores) - np.min(min_z_scores))
    norm_cluster_rmsd = (min_cluster_rmsd - np.min(min_cluster_rmsd)) / (np.max(min_cluster_rmsd) - np.min(min_cluster_rmsd))
    norm_inter_rmsd = (inter_min_rmsd - np.min(inter_min_rmsd)) / (np.max(inter_min_rmsd) - np.min(inter_min_rmsd))
    
    # Plot normalized data set
    fig = plt.figure(figsize=[10, 6], dpi=300)
    ax1 = fig.add_subplot(2,1,1)
    ax2 = fig.add_subplot(2,1,2)
    x1 = np.arange(len(labels[0:14]))
    x2 = np.arange(len(labels[14:]))
    width = 0.2
    rects1 = ax1.bar(x1 - 1.5*width, norm_sil_scores[0:14], width, label='SS')
    rects2 = ax1.bar(x1 - 0.5*width, norm_z_scores[0:14], width, label='EZG')
    rects3 = ax1.bar(x1 + 0.5*width, norm_cluster_rmsd[0:14], width, label='CRMSD')
    rects4 = ax1.bar(x1 + 1.5*width, norm_inter_rmsd[0:14], width, label='MIMR')

    rects5 = ax2.bar(x2 - 1.5*width, norm_sil_scores[14:], width, label='SS')
    rects6 = ax2.bar(x2 - 0.5*width, norm_z_scores[14:], width, label='EZG')
    rects7 = ax2.bar(x2 + 0.5*width, norm_cluster_rmsd[14:], width, label='CRMSD')
    rects8 = ax2.bar(x2 + 1.5*width, norm_inter_rmsd[14:], width, label='MIMR')
    
    # Plot text
    ax1.set_ylabel('Normalized Metrics')
    ax1.set_title(r'$\theta_{B}$ Scan Metric Comparison')
    ax1.set_xticks(x1)
    ax1.set_xticklabels(labels[0:14], rotation=45, ha="right")
    ax1.legend(bbox_to_anchor=(1.15, 0))
    ax1.grid(visible=True, axis="y")


    ax2.set_ylabel('Normalized Metrics')
    ax2.set_xticks(x2)
    ax2.set_xticklabels(labels[14:], rotation=45, ha="right")
    ax2.grid(visible=True, axis="y")

    ax1.set_axisbelow(True)
    ax2.set_axisbelow(True)
    plt.tight_layout()

    plt.savefig("bar_plot_test.png")
    plt.savefig("bar_plot_test.pdf")
    
    
if __name__ == "__main__":
    main()