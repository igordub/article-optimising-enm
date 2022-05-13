#! /bin/bash

# Generate BENMs with varying backbone stiffening
# coefficients
# 
# Usage: bash src/scan_benm.sh "8.00"

if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <dist-cutoff>"
    echo "Example: $0 "8.00""
    exit 1
fi

CUTOFF=$1

for BACK_COEF in 1 5 10 20 30 40 50 100 150 200
do
    for FORM_IDX in 0 1 2
    do
        printf -v BACK_COEF_PAD "%04d" ${BACK_COEF}

        WORK_DIR="tmp/benm/b${BACK_COEF_PAD}/${FORM_IDX}"
        OUTPUT_DIR="data/raw/scan-benm/b${BACK_COEF_PAD}/${FORM_IDX}"
        PDB_PATH="pdb/processed/0${FORM_IDX}.pdb"
        GENENMM_FLAGS="-f 1 -ca -het -b ${BACK_COEF} -mass -res -ccust ddpt.cfile -spcust ddpt.spfile"

        mkdir -p ${WORK_DIR}

        # Copy auxilary files, if any are present,
        # for -mass -ca -res, -ccust, -spcust and -fcust flags. 
        cp -f misc/{resmass.dat,ddpt.cfile,ddpt.ffile,ddpt.spfile} -t ${WORK_DIR} 2> /dev/null
        cp ${PDB_PATH} ${WORK_DIR}/struct.pdb

        pushd ${WORK_DIR}
        # Generate an EN for single-bead mass-weigthed ligands
        GENENMM -pdb struct.pdb -res -mass -ca -het -lig1
        sed -i 's/\(HETATM.\{6\}\) CA /\1 LG /' CAonly.pdb
        mv CAonly.pdb prot_lig_en.pdb

        # Create DDPT cfile
        printf '%4s %7.3f\n' " CA " "${CUTOFF}" > ddpt.cfile
        printf '%4s %7.3f\n' " LG " "0.0" >> ddpt.cfile

        popd

        bash src/run_enm.sh "${WORK_DIR}" "${OUTPUT_DIR}" "${WORK_DIR}/prot_lig_en.pdb" "${GENENMM_FLAGS}"


    done
done

rm -rf tmp/benm

echo "BENM simulation done!"
