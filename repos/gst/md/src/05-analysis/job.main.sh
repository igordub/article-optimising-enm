#!/bin/bash
# Clean, analyse and plot AMBER trajectory

#SBATCH --job-name=clean_traj
#SBATCH --time=00:05:00
#SBATCH --mem=75gb
#SBATCH --ntasks=1
#SBATCH --output=clean_traj.%j.log         
#SBATCH --account=phys-bioagg-2019


echo My working directory is `pwd`
echo Running job on host:
echo -e '\t'`hostname` at `date`
echo

conda init bash
source activate comp-biophys

python src/05-analysis/main.py

echo
echo Job completed at `date`