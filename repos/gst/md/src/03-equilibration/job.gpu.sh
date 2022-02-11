#!/bin/bash
# All-atom MD simulation: density equilibration

#SBATCH --job-name=eq
#SBATCH --time=00:01:00
#SBATCH --mem=1gb
#SBATCH --ntasks=1
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --output=eq.gpu.%j.log
#SBATCH --mail-type=ALL              
#SBATCH --mail-user=igors.dubanevics@york.ac.uk           
#SBATCH --account=phys-bioagg-2019

module load chem/Amber/16-foss-2018a-AmberTools-17-CUDA

SRC_DIR="src/03-equilibration"
LOGS_DIR="logs"
STRUC_DIR="data/00-structure"
INPUT_DIR="data/02-heating"
OUTPUT_DIR="data/03-equilibration"

echo My working directory is `pwd`
echo Running job on host:
echo -e '\t'`hostname` at `date`
echo

time pmemd.cuda -O \
    -i ${SRC_DIR}/eq.mdin \
    -c ${INPUT_DIR}/heat.ncrst \
    -p ${STRUC_DIR}/complex.parm7 \
    -ref ${INPUT_DIR}/heat.ncrst \
    -o ${OUTPUT_DIR}/eq.mdout \
    -r ${OUTPUT_DIR}/eq.ncrst \
    -x ${OUTPUT_DIR}/eq.nc
    -inf ${LOGS_DIR}/eq.mdinfo

echo
echo Job completed at `date`
