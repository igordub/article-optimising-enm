#!/bin/bash

# Prepares PDB file for AMBER
# AmebrTtools must be loaded and/or added to PATH

PDB_RAW_DIR="pdb/00-raw"
PDB_INT_DIR="pdb/01-interim"
PDB_PRO_DIR="pdb/02-processed"

# Remove MODEL pointers which brake `pdb4amber`
grep '^ATOM\|^HETATM\|TER' ${PDB_RAW_DIR}/complex.pdb | grep -v '^HETATM.\{11\}HOH' > ${PDB_INT_DIR}/complex.pdb
# grep '^ATOM\|^HETATM\|TER' ${PDB_RAW_DIR}/complex.pdb > ${PDB_INT_DIR}/complex.pdb

# Paste heteroatoms to the end of a PDB file
grep -v '^.\{21\}C' ${PDB_INT_DIR}/complex.pdb > ${PDB_INT_DIR}/protein.pdb
grep '^.\{21\}C' ${PDB_INT_DIR}/complex.pdb | grep -v '^TER' > ${PDB_INT_DIR}/ligands.pdb

# sed -i 's/^....../ATOM  /g' ${PDB_INT_DIR}/ligands.pdb | \
#     sed 's/^\(.\{17\}\).../\1 N3/g' > ${PDB_INT_DIR}/ligands.pdb

split --additional-suffix=.pdb -d -n 2 ${PDB_INT_DIR}/ligands.pdb ${PDB_INT_DIR}/ligand.model_

mv ${PDB_INT_DIR}/ligand.model_01.pdb ${PDB_INT_DIR}/ligand.model_02.pdb
mv ${PDB_INT_DIR}/ligand.model_00.pdb ${PDB_INT_DIR}/ligand.model_01.pdb

# sed -i 's/^\(.\{22\}\)..../\1   1/g' ${PDB_INT_DIR}/ligand.model_00.pdb
# sed -i 's/^\(.\{22\}\)..../\1   2/g' ${PDB_INT_DIR}/ligand.model_01.pdb

cat ${PDB_INT_DIR}/protein.pdb ${PDB_INT_DIR}/ligand.model_01.pdb ${PDB_INT_DIR}/ligand.model_02.pdb > ${PDB_INT_DIR}/complex.pdb
# cat ${PDB_INT_DIR}/protein.pdb ${PDB_INT_DIR}/ligands.pdb > ${PDB_INT_DIR}/complex.pdb

# Prepare PDB for AMBER and ENM sims
pdb4amber -i ${PDB_INT_DIR}/complex.pdb \
    -o ${PDB_INT_DIR}/complex.amber.pdb \
    --dry \
    --reduce \
    --model -1

pdb4amber -i ${PDB_INT_DIR}/complex.pdb \
    -o ${PDB_INT_DIR}/complex.enm.pdb \
    --dry \
    --nohyd \
    --model -1

cp ${PDB_INT_DIR}/complex.amber.pdb ${PDB_PRO_DIR}/complex.pdb
cp ${PDB_INT_DIR}/complex.enm.pdb ${PDB_PRO_DIR}/complex.enm.pdb

# Extract protein and ligand coordinates
grep -e '^ATOM\|^TER\|^END' ${PDB_INT_DIR}/complex.amber.pdb > ${PDB_PRO_DIR}/protein.pdb
grep -e '^HETATM' ${PDB_INT_DIR}/complex.amber.pdb > ${PDB_INT_DIR}/ligands.pdb 

# # Extract ligand separately
grep '^.\{23\}603' ${PDB_INT_DIR}/complex.amber.pdb | grep -v '^TER' > ${PDB_PRO_DIR}/ligand.model_01.pdb
grep '^.\{23\}604' ${PDB_INT_DIR}/complex.amber.pdb | grep -v '^TER\|^END' > ${PDB_PRO_DIR}/ligand.model_02.pdb
