#! /bin/bash

# Finds protein B-factor optimised distance cutoff
# 
# Usage: bash src/extract_dc.sh


INPUT_DIR="data/raw/scan-dc"
OUTPUT_DIR="data/interim/scan-dc"

mkdir -p ${OUTPUT_DIR}

BFACTOR_CORR_FILE="${OUTPUT_DIR}/bfactor_corr.csv"
echo "form_idx,dist_cutoff,bfactor_corr" > ${BFACTOR_CORR_FILE}

for CUTOFF in $(seq 7.5 0.5 12.0)
do
    printf -v CUTOFF_PAD "%05.2f" ${CUTOFF}

    for FORM_IDX in 0 1 2
    do
        INPUT_SUBDIR="${INPUT_DIR}/c${CUTOFF_PAD}/${FORM_IDX}"
        OUTPUT_SUBDIR="${OUTPUT_DIR}/c${CUTOFF_PAD}/${FORM_IDX}"

        mkdir -p ${OUTPUT_SUBDIR}


        # B-factors
        BFACTORS="${OUTPUT_SUBDIR}/bfactors.csv"
        echo "chain_id,residue_number,bfactor" > ${BFACTORS}
        # `awk` flied separator option is not working. Using `sed`, instead.
        grep -v '\#' ${INPUT_SUBDIR}/mode.bfactors | awk '{ print $5,$6,$9 }' | sed -e "s/ /,/g" >> ${BFACTORS}

        BFACTOR_CORR=$(head -n 1 ${INPUT_SUBDIR}/mode.bfactors | grep -o "[0-9]\.[0-9][0-9][0-9]")
        echo "${FORM_IDX},${CUTOFF},${BFACTOR_CORR}" >> ${BFACTOR_CORR_FILE}


        # Eigenvalues
        grep -v "^#" ${INPUT_SUBDIR}/mode.frequencies > ${OUTPUT_SUBDIR}/eigvals.csv
    done
done

echo "Distance cutoff scan value extraction complete!"
