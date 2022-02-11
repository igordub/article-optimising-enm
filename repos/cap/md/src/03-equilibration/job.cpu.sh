#!/bin/bash
# All-atom MD simulation: density equilibration

#SBATCH --job-name=eq
#SBATCH --time=02:00:00
#SBATCH --mem=10gb
#SBATCH --ntasks-per-node=32
#SBATCH --output=eq.cpu.%j.log
#SBATCH --mail-type=ALL             
#SBATCH --mail-user=igors.dubanevics@york.ac.uk           
#SBATCH --account=phys-bioagg-2019

module load chem/Amber/16-intel-2018b-AmberTools-17-patchlevel-10-15

SRC_DIR="src/03-equilibration"
LOGS_DIR="logs"
STRUC_DIR="data/00-structure"
INPUT_DIR="data/02-heating"
OUTPUT_DIR="data/03-equilibration"

echo My working directory is `pwd`
echo Running job on host:
echo -e '\t'`hostname` at `date`
echo

time mpirun -np $SLURM_NTASKS pmemd.MPI -O \
    -i ${SRC_DIR}/eq.mdin \
    -c ${INPUT_DIR}/heat.ncrst \
    -p ${STRUC_DIR}/complex.parm7 \
    -ref ${INPUT_DIR}/heat.ncrst \
    -o ${OUTPUT_DIR}/eq.mdout \
    -r ${OUTPUT_DIR}/eq.ncrst \
    -x ${OUTPUT_DIR}/eq.nc \
    -inf ${LOGS_DIR}/eq.mdinfo

echo
echo Job completed at `date`
