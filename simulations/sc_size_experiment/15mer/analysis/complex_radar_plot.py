import numpy as np
import matplotlib.pyplot as plt
import matplotlib.pylab as pl
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

ranges = [[0,  1], [0, 1], [0.12, 0], [0, 0.3]]   

colors = pl.cm.gist_rainbow(np.linspace(0,1,len(all_data)))
for i, data in enumerate(all_data):
    fig1 = plt.figure(figsize=(6, 6))
    label = str(round(float(data[0]), 3))
    data = [data[3], data[4], data[6], data[7]]
    data = [float(a) for a in data]
         
    # plotting
    radar = ComplexRadar(fig1, spoke_labels, ranges, n_ordinate_levels=4)
    radar.plot(data, color = colors[i], label=label)
    radar.ax.set_title(label)
    radar.fill(data, alpha=0.1, color = colors[i])
    plt.savefig("sc_" + label + "radar_plot.pdf")
    plt.savefig("sc_" + label + "radar_plot.jpg")
# radar.ax.legend(loc='center left', bbox_to_anchor=(1.2, 0.5))
plt.show()
