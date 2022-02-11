#! /bin/bash
# Performs ANM distance cutoff scan.
# 
# Usage: bash src/scan_dc.sh


OUTPUT_DIR="data/raw/scan-dc"

mkdir -p ${OUTPUT_DIR}

for CUTOFF in $(seq 5.0 0.5 15.0)
do
    for FORM_IDX in 0 1 2
    do
        printf -v CUTOFF_PAD "%05.2f" ${CUTOFF}
        CUTOFF_DIR="${OUTPUT_DIR}/c${CUTOFF_PAD}/${FORM_IDX}"

        WORK_DIR="tmp/scan-dc/c${CUTOFF_PAD}/${FORM_IDX}"
        PDB_PATH="pdb/processed/0${FORM_IDX}.pdb"
        GENENMM_FLAGS="-f 1 -ca -c ${CUTOFF} -het -mass -res -lig1"

        bash src/run_enm.sh "${WORK_DIR}" "${CUTOFF_DIR}" "${PDB_PATH}" "${GENENMM_FLAGS}"
    done    
done

echo "Distance cutoff scan complete!"
