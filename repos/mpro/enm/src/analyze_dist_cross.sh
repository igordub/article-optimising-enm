#! /bin/bash

# Calculate alpha-carbon - alpha-carbon distance matrix
# and cross-correlation for given mode interval.
# 
# Usage: bash src/analyze_dist_cross.sh

if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <mode_start> <mode_end>"
    echo "Example: $0 "11" "13""
    exit 1
fi

MODE_START=$1
MODE_END=$2

INPUT_DIR="data/raw/dist_cross"
OUTPUT_DIR="data/interim/dist_cross"

mkdir -p ${OUTPUT_DIR}

for MODE_NUM in $(seq "${MODE_START}" 1 "${MODE_END}")
do
    # Use non-trivial modes
    printf -v MODE_NUM_PAD '%04d' ${MODE_NUM}

    python src/plot_dist_div.py -cross ${INPUT_DIR}/crosscor.m${MODE_NUM_PAD}.dat \
        -dist ${INPUT_DIR}/dist.dat \
        -suffix "m${MODE_NUM_PAD}" \
        -outdir "${OUTPUT_DIR}" \
        -title "Distance and divergance | CAP | Mode ${MODE_NUM}"    

done

echo "Distance and divergance analysis is done!"
