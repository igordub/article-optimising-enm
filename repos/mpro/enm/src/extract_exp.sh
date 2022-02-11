#! /bin/bash

# Extracts inforamtion from PDB file.
# 
# Usage: bash src/extract_exp.sh

INT_DATA_DIR="data/interim"
PDB_DIR="pdb/processed/"

# Get experimental B-factors for CA atoms
BFACTORS_EXP="${INT_DATA_DIR}/exp.bfactors.csv"
echo "form_idx,chain_id,residue_number,bfactor" > ${BFACTORS_EXP}

BFACTORS_TMP="tmp/bfactors.csv"
for FORM_IDX in 0 1 2
do
    # Add white space between ocupancy and B-factor columns in case if B-factor has 5 digits
    # `awk` flied separator option is not working. Using `sed`, instead.
    grep '^ATOM.\{8\} CA ' ${PDB_DIR}/0${FORM_IDX}.pdb | sed 's/\(.\{60\}\)\(.*\)/\1 \2/' \
        | awk '{ print $5,$6,$11 }' | sed -e "s/ /,/g" > ${BFACTORS_TMP}
    
    LINE_NUM=$(cat ${BFACTORS_TMP} | wc -l)
    yes ${FORM_IDX} | head -n ${LINE_NUM} > tmp/form_idx.csv

    paste -d "," tmp/form_idx.csv ${BFACTORS_TMP} >> ${BFACTORS_EXP}
done

rm tmp/form_idx.csv ${BFACTORS_TMP}

echo "PDB info extracted!"
