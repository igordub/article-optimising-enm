#! /bin/bash

# Process CAP (PDB ID: 4HZF) to mathc AMBER processed
# PDB file

INPUT_DIR="pdb/external/"
OUTPUT_DIR="tmp"

grep '^ATOM\|^TER' ${INPUT_DIR}/4hzf.pdb > ${OUTPUT_DIR}/4hzf.enm.pdb
grep '^HETATM.\{10\} CMP A 301 \|^HETATM.\{10\} CMP B 301 ' ${INPUT_DIR}/4hzf.pdb >> ${OUTPUT_DIR}/4hzf.enm.pdb
echo "END" >> ${OUTPUT_DIR}/4hzf.enm.pdb

# Keep only A conformation
sed -ni '/^.\{16\}[B-Z]/!p' ${OUTPUT_DIR}/4hzf.enm.pdb 

# Trim N-termini
sed -ni '/^ATOM.\{18\}   [1-9]/!p' ${OUTPUT_DIR}/4hzf.enm.pdb
sed -ni '/^ATOM.\{18\}  10/!p' ${OUTPUT_DIR}/4hzf.enm.pdb
# Trim C-termini
sed -ni '/^ATOM.\{18\} 20[8-9]/!p' ${OUTPUT_DIR}/4hzf.enm.pdb
