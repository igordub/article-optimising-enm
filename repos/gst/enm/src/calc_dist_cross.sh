#! /bin/bash

# Calculate alpha-carbon - alpha-carbon distance matrix
# and cross-correlation for given mode interval.
# 
# Usage: bash src/calc_dist_cross.sh

if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <dist_cutoff> <mode_start> <mode_end>"
    echo "Example: $0 "8.00" "11" "13""
    exit 1
fi

# if [ $2 > $3 ]; then
#     echo "Starting mode must be <= end mode."
#     exit 1
# fi

CUTOFF=$1
MODE_START=$2
MODE_END=$3

printf -v CUTOFF_PAD "%05.2f" ${CUTOFF}

INPUT_DIR="data/raw/scan-dc/c${CUTOFF_PAD}/0"

if ! [ -d "${INPUT_DIR}" ]; then
  # Take action if ${INPUT_DIR} doesn't exist.
  echo "Check if specified cutoff was simualted in\ndistance cutoff scan."
  exit 1
fi

WORK_DIR="tmp/dist_cross"
OUTPUT_DIR="data/raw/dist_cross"

mkdir -p ${WORK_DIR} ${OUTPUT_DIR}

cp ${INPUT_DIR}/CAonly.pdb ${INPUT_DIR}/matrix.eigenfacs -t ${WORK_DIR}

pushd ${WORK_DIR}

echo "Running: SPACING"
SPACING -pdb CAonly.pdb

echo "Running: CROSSCOR"
for MODE_NUM in $(seq "${MODE_START}" 1 "${MODE_END}")
do
    # Use non-trivial modes
    printf -v MODE_NUM_PAD '%04d' ${MODE_NUM}

    TRIV_MODE_NUM=$((${MODE_NUM} + 6))
    CROSCOR -i matrix.eigenfacs -s ${TRIV_MODE_NUM} -e ${TRIV_MODE_NUM}

    mv crosscor.dat crosscor.m${MODE_NUM_PAD}.dat
done

popd

mv  ${WORK_DIR}/crosscor.m*.dat \
    ${WORK_DIR}/dist.dat \
    -t ${OUTPUT_DIR}

rm -rf ${WORK_DIR}

echo "Distance and cross-correlation calculation is done!"
