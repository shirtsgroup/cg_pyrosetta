import mdtraj as md
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import sys


def plot_dihedral_histogram(filename, bins):
    res = md.load_pdb(filename)
    model_num = int(filename[2])
    print(model_num)
    selector = 'name '+' '.join(['BB'+str(i) for i in range(1, model_num+1)])
    indexes = res.topology.select(selector)
    dihes_id = {i:[] for i in range(3)}
    for i in range(len(indexes)-3):
        dihes_id[i%model_num].append(indexes[i:i+4])
    

    dihes = {}
    for j in range(model_num):
        dihes_id[j] = np.array(dihes_id[j])
        dihes[j] = md.compute_dihedrals(res, dihes_id[j]).flatten()
        print(filename, 'Median :', np.median(dihes[j]))
        print(filename, 'Variance :', np.var(dihes[j]))
        plt.hist(dihes[j], bins = bins, density=True, range=[-np.pi, np.pi])
        plt.ylabel('Probability Density')
        plt.xlabel('Dihedral Angle (Radians)')
        plt.xlim([-np.pi, np.pi])
        plt.title(filename+' Dihedral Angle '+str(j+1))
        plt.savefig(filename.split('.')[0]+'_'+str(j)+'.png')
        plt.show()
    

def main():
    for name in sys.argv[2:]:
        plot_dihedral_histogram(filename=name, bins=int(sys.argv[1]))
    # plot_dihedral_histogram(filename='CG11_example_angles_0.pdb', bins=7)




if __name__ == '__main__':
    main()
