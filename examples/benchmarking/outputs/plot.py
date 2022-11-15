# I this file will take in the time_data.yml files and plot

import yaml
import matplotlib.pyplot as plt
import sys
import numpy as np

def read_yaml_file(filename):
    with open(filename, 'r') as f:
        data = yaml.load(f)
    return(data)

def array_to_dist(array):
    return(np.mean(array), np.std(array))

def main():
    fn = sys.argv[1]
    data = read_yaml_file(fn)
    names = ['CG11', 'CG21', 'CG31']
    keys = dct.keys()

    for name in names:
        spec_keys = list(filter(lambda ns: name in ns, names))
        for sk in spec_keys:
            dist = data[sk]
            mean, std = array_to_dist(dist)
            
    

if __name__ == "__main__":
    main()