#!/bin/bash
# Downlaods raw PDB file to pdb/00-raw/
if [ "$#" -ne 1 ]; then
    echo "Incorrect variable number" 
    echo "Usage: $0 <pdb-id>"
    echo "Example: $0 4hzf"
    exit 1
fi

PDB_ID=$1

wget "https://files.rcsb.org/view/${PDB_ID}.pdb1" -O pdb/00-raw/complex.pdb
