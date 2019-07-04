#!/bin/bash
#SBATCH -t 1-12:00:00
#SBATCH --output=logfile
#SBATCH --nodes=1
#SBATCH --ntasks-per-node 10
#SBATCH --partition RM-shared

module load anaconda5
`
python ../../setup.py install --user
