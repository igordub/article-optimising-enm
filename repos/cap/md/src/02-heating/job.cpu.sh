#!/bin/bash
# All-atom MD simulation: heating

#SBATCH --job-name=heat
#SBATCH --time=00:30:00
#SBATCH --mem=5gb
#SBATCH --ntasks-per-node=32
#SBATCH --output=heat.cpu.%j.log
#SBATCH --mail-type=ALL             
#SBATCH --mail-user=igors.dubanevics@york.ac.uk           
#SBATCH --account=phys-bioagg-2019

module load chem/Amber/16-intel-2018b-AmberTools-17-patchlevel-10-15

SRC_MIN="src/02-heating"
LOGS_DIR="logs"
STRUC_DIR="data/00-structure"
INPUT_DIR="data/01-minimisation"
OUTPUT_DIR="data/02-heating"

echo My working directory is `pwd`
echo Running job on host:
echo -e '\t'`hostname` at `date`
echo

time mpirun -np $SLURM_NTASKS pmemd.MPI -O \
    -i ${SRC_MIN}/heat.mdin \
    -c ${INPUT_DIR}/min.02.ncrst \
    -p ${STRUC_DIR}/complex.parm7 \
    -ref ${INPUT_DIR}/min.02.ncrst \
    -o ${OUTPUT_DIR}/heat.mdout \
    -r ${OUTPUT_DIR}/heat.ncrst \
    -x ${OUTPUT_DIR}/heat.nc
    -inf ${LOGS_DIR}/heat.mdinfo

echo
echo Job completed at `date`
