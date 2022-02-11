#! /bin/bash

# Generate BENMs with varying backbone stiffening
# coefficients
# 
# Usage: bash src/run_benm.sh "8.00"

if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <dist-cutoff>"
    echo "Example: $0 "8.00""
    exit 1
fi

CUTOFF=$1

# for BACK_COEF in 1 10 50 100 150 200 250 300
for BACK_COEF in 20 30 40
do
    for FORM_IDX in 0 1 2
    do
        printf -v BACK_COEF_PAD "%04d" ${BACK_COEF}

        WORK_DIR="tmp/benm/b${BACK_COEF_PAD}/${FORM_IDX}"
        OUTPUT_DIR="data/raw/scan-benm/b${BACK_COEF_PAD}/${FORM_IDX}"
        PDB_PATH="pdb/processed/0${FORM_IDX}.pdb"
        GENENMM_FLAGS="-f 1 -ca -c ${CUTOFF} -b ${BACK_COEF} -het -res -mass -lig1"

        bash src/run_enm.sh "${WORK_DIR}" "${OUTPUT_DIR}" "${PDB_PATH}" "${GENENMM_FLAGS}"

    done
done

rm -rf tmp/benm

echo "BENM simulation done!"
