#! /bin/bash

# Copies figures and tables to the article's repository

OUTPUT_DIR="/users/id583/scratch/repos/article-jmb-si-allostery"

cp scratch/*.png ${OUTPUT_DIR}/img/mpro/
cp scratch/*.tex ${OUTPUT_DIR}/tables/mpro/
