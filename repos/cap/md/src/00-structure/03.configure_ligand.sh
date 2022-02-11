#!/bin/bash

# Creates force field modification for the ligand
# Each ligand must be treated separately

PDB_RAW_DIR="pdb/00-raw"
PDB_INT_DIR="pdb/01-interim"
PDB_PRO_DIR="pdb/02-processed"
ANTECHAM_DIR="data/antechamber"

# Ligand name: cAMP
# Ligand charge: -1
# Multiplicity: 2

# Run antechamber 
antechamber -i ${PDB_PRO_DIR}/ligand.chain_A.pdb -fi pdb -o ${ANTECHAM_DIR}/ligand.chainA.mol2 -fo mol2 \
    -s 2 -c bcc -at gaff2 \
    -nc -1 \
    -m 2
rm ANTECHAMBER_*.AC ANTECHAMBER_*.AC0 ATOMTYPE.INF sqm.*

antechamber -i ${PDB_PRO_DIR}/ligand.chain_B.pdb -fi pdb -o ${ANTECHAM_DIR}/ligand.chainB.mol2 -fo mol2 \
    -s 2 -c bcc -at gaff2 \
    -nc -1 \
    -m 2
rm ANTECHAMBER_*.AC ANTECHAMBER_*.AC0 ATOMTYPE.INF sqm.*

# Run parmch2
parmchk2 -i ${ANTECHAM_DIR}/ligand.chainA.mol2 -f mol2 -o ${ANTECHAM_DIR}/ligand.chainA.frcmod -s gaff2
parmchk2 -i ${ANTECHAM_DIR}/ligand.chainB.mol2 -f mol2 -o ${ANTECHAM_DIR}/ligand.chainB.frcmod -s gaff2
