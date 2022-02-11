#!/bin/bash
# All-atom MD simulation: minimisation

#SBATCH --job-name=min
#SBATCH --time=00:03:00
#SBATCH --mem=1gb
#SBATCH --ntasks=1
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --output=min.gpu.%j.log
#SBATCH --mail-type=ALL              
#SBATCH --mail-user=igors.dubanevics@york.ac.uk           
#SBATCH --account=phys-bioagg-2019

module load chem/Amber/16-foss-2018a-AmberTools-17-CUDA

SRC_DIR="src/01-minimisation"
LOGS_DIR="logs"
INPUT_DIR="data/00-structure"
OUTPUT_DIR="data/01-minimisation"

echo My working directory is `pwd`
echo Running job on host:
echo -e '\t'`hostname` at `date`
echo

# Minimise the solvent
time pmemd.cuda -O \
    -i ${SRC_DIR}/min.01.mdin \
    -p ${INPUT_DIR}/complex.parm7 \
    -c ${INPUT_DIR}/complex.ncrst \
    -ref ${INPUT_DIR}/complex.ncrst \
    -inf ${LOGS_DIR}/min.01.mdinfo \
    -o ${OUTPUT_DIR}/min.01.mdout \
    -r ${OUTPUT_DIR}/min.01.ncrst \
    -x ${OUTPUT_DIR}/min.01.nc
   
# Minimise the whole system
time pmemd.cuda -O \
    -i ${SRC_DIR}/min.02.mdin \
    -c ${OUTPUT_DIR}/min.01.ncrst \
    -p ${INPUT_DIR}/complex.parm7 \
    -ref ${OUTPUT_DIR}/min.01.ncrst \
    -inf ${LOGS_DIR}/min.02.mdinfo \
    -o ${OUTPUT_DIR}/min.02.mdout \
    -r ${OUTPUT_DIR}/min.02.ncrst \
    -x ${OUTPUT_DIR}/min.02.nc 

echo
echo Job completed at `date`
