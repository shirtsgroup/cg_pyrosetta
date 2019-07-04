#!/bin/sh
#SBATCH --time=1-13:00:00
#SBATCH --partition RM-shared
#SBATCH --job-name=ParameterScan
#SBATCH --nodes=1
#SBATCH --ntasks-per-node 10

# load in conda environment needed to run script
module load anaconda5

# build conda environment from file
conda env create -f environment.yml
source activate python36

echo "SUBMITTED" > jobstatus.txt
echo `date` > jobstatus.txt
# Run Python Script

python parameter_scan.py

echo "FINISHED" > jobstatus.txt
echo `date` > jobstatus.txt
