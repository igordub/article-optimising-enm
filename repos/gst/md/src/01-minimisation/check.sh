#!/bin/bash
# Check equilibration MD simulation results

STRUC_DIR="data/00-structure"
MIN_DIR="data/01-minimisation"
HEAT_DIR="data/02-heating"
EQ_DIR="data/03-equilibration"
PROD_DIR="data/04-production"
TMP_DIR="tmp"

module load chem/Amber/16-intel-2018b-AmberTools-17-patchlevel-10-15

ambpdb -p ${STRUC_DIR}/complex.parm7 -c ${MIN_DIR}/min.01.ncrst > ${TMP_DIR}/min.01.pdb
ambpdb -p ${STRUC_DIR}/complex.parm7 -c ${MIN_DIR}/min.02.ncrst > ${TMP_DIR}/min.02.pdb
