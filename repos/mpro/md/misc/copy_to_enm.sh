#!/bin/bash
# Copies eigenvalues and B-factors
# to ENM repo

ENM_DIR="/users/id583/scratch/repos/enm-scov2-mpro/"

cp output/eigvals.csv ${ENM_DIR}/data/external/md.eigvals.00.csv
cp output/bfactors.csv ${ENM_DIR}/data/external/md.bfactors.00.csv
