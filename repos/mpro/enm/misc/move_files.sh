#! /bin/bash

OUTPUT_DIR="data/raw/find-dist-cutoff"
STASH_DIR="tmp/stash"

mkdir -p ${STASH_DIR}

for CUTOFF in $(seq 5.5 0.5 15.0)
do
    printf -v CUTOFF_PAD "%05.2f" ${CUTOFF}
    OLD_DIR="${OUTPUT_DIR}/c${CUTOFF_PAD}"
    NEW_DIR="${OUTPUT_DIR}/c${CUTOFF_PAD}/0"

    mv ${OLD_DIR}/* ${STASH_DIR}

    mkdir -p ${NEW_DIR}

    mv ${STASH_DIR}/* ${NEW_DIR}
done

rm -rf ${STASH_DIR}
