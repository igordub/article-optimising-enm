#! /bin/bash

# Generate BENMs with varying backbone stiffening
# coefficients
# 
# Usage: bash src/extract_benm.sh


# for BACK_COEF in 1 10 20 30 40 50 100 150 200 250 300
# do
#     for FORM_IDX in 0 1 2
#     do
#         printf -v BACK_COEF_PAD "%04d" ${BACK_COEF}

#         INPUT_DIR="data/raw/scan-benm/b${BACK_COEF_PAD}/${FORM_IDX}"
#         OUTPUT_DIR="data/interim/scan-benm/b${BACK_COEF_PAD}/${FORM_IDX}"
#         mkdir -p ${OUTPUT_DIR}

#         #  Extract eigenvalues
#         grep -v "^#" ${INPUT_DIR}/mode.frequencies > ${OUTPUT_DIR}/eigvals.csv

#     done
# done

# Glob-way of data processing
for INPUT_DIR in data/raw/scan-benm/b????/?
do
    OUTPUT_DIR=$(echo "${INPUT_DIR}" | sed 's/data\/raw\//data\/interim\//' )
    mkdir -p ${OUTPUT_DIR}

    # Extract eigenvalues
    grep -v "^#" ${INPUT_DIR}/mode.frequencies > ${OUTPUT_DIR}/eigvals.csv
done

echo "BENMs scan eigenvalue extraction complete!"
