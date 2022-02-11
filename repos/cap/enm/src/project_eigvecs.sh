#! /bin/bash

# Project eigenvectors and replace experimental B-factors 
# in PDB file with ENM computed B-factors
# 
# Usage: bash src/project_eigvecs.sh

if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <dist-cutoff>"
    echo "Example: $0 "8.00""
    exit 1
fi

CUTOFF=$1
printf -v CUTOFF_PAD "%05.2f" ${CUTOFF}

INPUT_DIR="data/raw/scan-dc/c${CUTOFF_PAD}/0"

if ! [ -d "${INPUT_DIR}" ]; then
  # Take action if ${INPUT_DIR} doesn't exist.
  echo "Check if specified cutoff was simualted in\ndistance cutoff scan."
  exit 1
fi

WORK_DIR="tmp/project_eigvecs"
OUTPUT_DIR="data/raw/project_eigvecs"

mkdir -p ${WORK_DIR} ${OUTPUT_DIR}

cp ${INPUT_DIR}/CAonly.pdb ${INPUT_DIR}/matrix.eigenfacs -t ${WORK_DIR}

pushd ${WORK_DIR}

echo "Running: PROJECT"

for MODE_NUM in {1..25}
do
    # Use non-trivial modes
    printf -v MODE_NUM_PAD '%04d' ${MODE_NUM}

    TRIV_MODE_NUM=$((${MODE_NUM} + 6))
    PROJECT -i matrix.eigenfacs -pdb CAonly.pdb -s ${TRIV_MODE_NUM} -e ${TRIV_MODE_NUM} -het -scale 0.5

    RMSCOL -i matrix.eigenfacs -pdb CAonly.pdb -s ${TRIV_MODE_NUM} -e ${TRIV_MODE_NUM} -het > /dev/null

    # Extract residue numbers and RMSD values
    grep -v '^#' rms.val | awk '{ print $6,$7 }' > resnum.rms
    # Normalize RMSD
    awk 'FNR==NR{max=($2+0>max)?$2:max;next} {print $1,$2/max}' resnum.rms resnum.rms > resnum.rms_norm

    PLOTPDB -i resnum.rms_norm -pdb CAonly.pdb
    
    mv Mode_???.pdb mode.m${MODE_NUM_PAD}.pdb
    cp plot.pdb rmsd.m${MODE_NUM_PAD}.pdb
done

popd

mv  ${WORK_DIR}/mode.m*.pdb \
    ${WORK_DIR}/rmsd.m*.pdb \
    -t ${OUTPUT_DIR}

rm -rf ${WORK_DIR}

echo "Eigenvector projection complete!"
