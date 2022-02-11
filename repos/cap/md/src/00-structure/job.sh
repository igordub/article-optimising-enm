#!/bin/bash

#SBATCH --job-name=config
#SBATCH --time=00:05:00
#SBATCH --mem=1gb
#SBATCH --ntasks=1
#SBATCH --output=config.%j.log
#SBATCH --mail-type=BEGIN,END,FAIL              
#SBATCH --mail-user=igors.dubanevics@york.ac.uk           
#SBATCH --account=phys-covid19-2020

module load chem/Amber/16-intel-2018b-AmberTools-17-patchlevel-10-15

echo My working directory is `pwd`
echo Running job on host:
echo -e '\t'`hostname` at `date`
echo

time bash src/00-structure/01.download_pdb.sh
time bash src/00-structure/02.prepare_pdb.sh
time bash src/00-structure/03.configure_ligand.sh
time tleap -f src/00-structure/tleap.04.leapin

echo
echo Job completed at `date`
