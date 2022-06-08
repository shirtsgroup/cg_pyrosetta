import numpy as np
import matplotlib.pyplot as plt
import matplotlib.pylab as pl
import os
import shutil as sh
import csv
import seaborn as sns # improves plot aesthetics

sns.set_theme()

def _invert(x, limits):
    """inverts a value x on a scale from
    limits[0] to limits[1]"""
    return limits[1] - (x - limits[0])

def _scale_data(data, ranges):
    """scales data[1:] to ranges[0],
    inverts if the scale is reversed"""
    for d, (y1, y2) in zip(data[1:], ranges[1:]):
        assert (y1 <= d <= y2) or (y2 <= d <= y1)
    x1, x2 = ranges[0]
    d = data[0]
    if x1 > x2:
        d = _invert(d, (x1, x2))
        x1, x2 = x2, x1
    sdata = [d]
    for d, (y1, y2) in zip(data[1:], ranges[1:]):
        if y1 > y2:
            d = _invert(d, (y1, y2))
            y1, y2 = y2, y1
        sdata.append((d-y1) / (y2-y1) 
                     * (x2 - x1) + x1)
    return sdata

class ComplexRadar():
    def __init__(self, fig, variables, ranges,
                 n_ordinate_levels=6):
        angles = np.arange(0, 360, 360./len(variables))

        axes = [fig.add_axes([0.1,0.1,0.9,0.9],polar=True,
                label = "axes{}".format(i)) 
                for i in range(len(variables))]
        l, text = axes[0].set_thetagrids(angles, 
                                         labels=variables)
        [txt.set_rotation(angle-90) for txt, angle 
             in zip(text, angles)]
        for ax in axes[1:]:
            ax.patch.set_visible(False)
            ax.grid("off")
            ax.xaxis.set_visible(False)
        for i, ax in enumerate(axes):
            grid = np.linspace(*ranges[i], 
                               num=n_ordinate_levels)
            gridlabel = ["{}".format(round(x,2)) 
                         for x in grid]
            if ranges[i][0] > ranges[i][1]:
                grid = grid[::1] # hack to invert grid
                          # gridlabels aren't reversed
            gridlabel[0] = "" # clean up origin
            ax.set_rgrids(grid, labels=gridlabel,
                         angle=angles[i])
            #ax.spines["polar"].set_visible(False)
            ax.set_ylim(*ranges[i])
        # variables for plotting
        self.angle = np.deg2rad(np.r_[angles, angles[0]])
        self.ranges = ranges
        self.ax = axes[0]
    def plot(self, data, *args, **kw):
        sdata = _scale_data(data, self.ranges)
        self.ax.plot(self.angle, np.r_[sdata, sdata[0]], *args, **kw)
    def fill(self, data, *args, **kw):
        sdata = _scale_data(data, self.ranges)
        self.ax.fill(self.angle, np.r_[sdata, sdata[0]], *args, **kw)


# example data
with open("summary_stats.csv", newline="") as csvfile:
    reader = csv.reader(csvfile)
    all_data = list(reader)

spoke_labels = all_data.pop(0)
# spoke_labels = [spoke_labels[3], spoke_labels[4], spoke_labels[6], spoke_labels[7]]
spoke_labels = ["", "", "", ""]

legend_labels = [str(round(float(a[0]),3)) for a in all_data]
columns = [4, 5, 7, 8]
no_characters = [x for x in all_data if "-" not in x]
ranges = [[0, float(max(no_characters, key=lambda x: x[i])[i])] for i in columns]

# Flip ranges for min_energy_cluster_rmsd and multipy rmsds by 10 to be consistent with sigma
ranges[2][0], ranges[2][1] = 10*ranges[2][1], 10*ranges[2][0]
ranges[3][0], ranges[3][1] = 10*ranges[3][0], 10*ranges[3][1]

colors = pl.cm.tab20(np.linspace(0,1,len(all_data)))

# Make directory to store plots + look up tables
if not os.path.isdir("radar_plots"):
    os.mkdir("radar_plots")
else:
    sh.rmtree("radar_plots")
    os.mkdir("radar_plots")

for i, data in enumerate(all_data):
    label = str(round(float(data[0]), 3))
    if "-" in data:
        print("Skipping", label)
        continue
    fig1 = plt.figure(figsize=(6, 6))
    data = [data[4], data[5], data[7], data[8]]
    data = [float(a) for a in data]
    # Scale RMSD values so that they're consistent with sigma
    data[2] *= 10
    data[3] *= 10
         
    # plotting
    radar = ComplexRadar(fig1, spoke_labels, ranges, n_ordinate_levels=4)
    radar.plot(data, color = "black", label=label)
    radar.ax.set_title(label)
    radar.fill(data, alpha=0.4, color = colors[i])
    plt.savefig(os.path.join("radar_plots", "ba_" + label + "_radar_plot.pdf"))
    plt.savefig(os.path.join("radar_plots", "ba_" + label + "_radar_plot.jpg"))
    plt.close()
# radar.ax.legend(loc='center left', bbox_to_anchor=(1.2, 0.5))