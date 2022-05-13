#! /bin/bash
# Performs ANM distance cutoff scan.
# 
# Usage: bash src/scan_dc.sh

OUTPUT_DIR="data/raw/scan-dc"

mkdir -p ${OUTPUT_DIR}

# for CUTOFF in $(seq 7.5 0.5 15.0)
# for CUTOFF in $(seq 7.5 0.5 7.5)
for CUTOFF in $(seq 7.5 0.5 12.0)
do
    for FORM_IDX in 0 1 2
    # for FORM_IDX in 2
    do
        printf -v CUTOFF_PAD "%05.2f" ${CUTOFF}
        CUTOFF_DIR="${OUTPUT_DIR}/c${CUTOFF_PAD}/${FORM_IDX}"

        WORK_DIR="tmp/scan-dc/c${CUTOFF_PAD}/${FORM_IDX}"
        PDB_PATH="pdb/processed/0${FORM_IDX}.pdb"
        GENENMM_FLAGS="-f 1 -ca -het -mass -res -ccust ddpt.cfile -spcust ddpt.spfile"

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

        bash src/run_enm.sh "${WORK_DIR}" "${CUTOFF_DIR}" "${WORK_DIR}/prot_lig_en.pdb" "${GENENMM_FLAGS}"
    done    
done

echo "Distance cutoff scan complete!"
