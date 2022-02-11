#!/bin/bash
# Splits production run into several chunks
# Times can be changed by varying `nstlim` in `prod.mdin`
# and changing `MD_END_JOB`
#
# Good practice is to test the script by setting `MD_END_JOB` = `MD_START_JOB`

if [ "$#" -ne 2 ]; then
   echo "Supply start and end run numbers " 
   echo "Usage: $0 <MD_START_JOB> <MD_END_JOB>"
   echo "Example: $0 5 10"
   exit 1
fi

MD_START_JOB=$1
MD_END_JOB=$2

SRC_DIR="src/04-production"
INPUT_DIR="data/03-equilibration"
OUTPUT_DIR="data/04-production"
STRUC_PATH="data/00-structure"
LOGS_DIR="logs"

# Make the density equilibration trajectory file
# the first input file for the production
if [ "${MD_START_JOB}" -eq 1 ]; then
   printf -v INPUT_FILENAME "prod.%03d" "0"
   cp ${INPUT_DIR}/eq.ncrst ${OUTPUT_DIR}/${INPUT_FILENAME}.ncrst
fi

echo -n "Starting script at: "
date
echo ""

for MD_CURRENT_JOB in $(seq $MD_START_JOB $MD_END_JOB)
do
   echo -n "Job $MD_CURRENT_JOB started at: "
   date
   MD_INPUT=$((${MD_CURRENT_JOB} - 1))
   printf -v INPUT_FILENAME "prod.%03d" ${MD_INPUT}
   printf -v OUTPUT_FILENAME "prod.%03d" ${MD_CURRENT_JOB}
   pmemd.cuda -O \
      -i ${SRC_DIR}/prod.mdin \
      -o ${OUTPUT_DIR}/${OUTPUT_FILENAME}.mdout \
      -p ${STRUC_PATH}/complex.parm7 \
      -c ${OUTPUT_DIR}/${INPUT_FILENAME}.ncrst \
      -r ${OUTPUT_DIR}/${OUTPUT_FILENAME}.ncrst \
      -x ${OUTPUT_DIR}/${OUTPUT_FILENAME}.nc \
      -inf ${LOGS_DIR}/${OUTPUT_FILENAME}.mdinfo

   echo -n "Job $MD_CURRENT_JOB finished at: "
   date
done

echo "ALL DONE"
