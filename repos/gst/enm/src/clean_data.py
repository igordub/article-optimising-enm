#!/usr/bin/env python
""" 
    This script cleans the data needed for analysis,
    putting it into the same format and merging it
    together.
"""
import os
import glob

import pandas as pd
import numpy as np
from os.path import join as join_paths, basename as get_basename
import subprocess
import utilities as utils


def clean_dict(data):
    """ Reads interim eigenvalues from BENM's scan 
        and writes them to a single file.
    """
    df = pd.DataFrame(data)

    mode_nums = np.arange(1, df.shape[0]+1)
    df['mode'] = mode_nums
    df = df.set_index('mode', drop=True)

    return df

def clean_exp():
    """ Clean experimental data from PDB.
    """
    subprocess.run(['bash', 'src/extract_exp.sh'])

    return None

def clean_md():
    """ Clean all-atom MD data.
    """
    config = utils.read_config()
    input_dir = config['data']['extFilePath']
    output_dir = config['data']['intFilePath']

    eigvals_dict = {}
    for form_idx in range(1):
        file_path = join_paths(input_dir, 'md.eigvals.0{}.csv'.format(form_idx))

        eigvals_dict = {form_idx: pd.read_csv(file_path, header=0, index_col='mode')['eigval']}

    # Combine all dataframes into one
    # and add another index
    eigvlas_md = pd.concat(eigvals_dict,
        ignore_index=False)
    eigvlas_md.index = eigvlas_md.index.set_names(['form_idx','mode'])

    # Save eigenvalues
    eigvlas_md.to_csv(join_paths(output_dir, 'md.eigvals.csv'),
        float_format='%.6f')

    bfactors_dict = {}
    for form_idx in range(1):
        file_path = join_paths(input_dir, 'md.bfactors.0{}.csv'.format(form_idx))

        bfactors_dict = {form_idx: pd.read_csv(file_path, header=0, index_col='residue_number')['bfactor']}

    bfactors_md = pd.concat(bfactors_dict,
        ignore_index=False)
    bfactors_md.index = bfactors_md.index.set_names(['form_idx','residue_number'])

    # Save eigenvalues
    bfactors_md.to_csv(join_paths(output_dir, 'md.bfactors.csv'),
        float_format='%.2f')


def clean_dc():
    """ Clean distance cutoff scan.
    """
    # Extract eigenvalues
    subprocess.run(['bash', 'src/extract_dc.sh'])

    config = utils.read_config()
    input_dir = join_paths(config['data']['intFilePath'], 'scan-dc')
    output_dir = config['data']['intFilePath']

    # Eigenvalues
    eigvlas_dc = {}
    for form_idx in range(3):
        file_paths = sorted(glob.glob(
            join_paths(input_dir,
                       'c??.??',
                       str(form_idx),
                       'eigvals.csv')))


        data_dict = {
            float(filename.split('/')[-3].replace('c', '')): np.loadtxt(filename)
            for filename in file_paths}

        eigvlas_dc[form_idx] = clean_dict(data_dict)

    # Save eigenvalues
    eigvlas_dc[0].to_csv(join_paths(
        output_dir, 'dc.eigvals.00.csv'), float_format='%.6f')
    eigvlas_dc[1].to_csv(join_paths(
        output_dir, 'dc.eigvals.01.csv'), float_format='%.6f')
    eigvlas_dc[2].to_csv(join_paths(
        output_dir, 'dc.eigvals.02.csv'), float_format='%.6f')

    # B-factors
    bfactors_dc = {}
    for form_idx in range(3):
        file_paths = sorted(glob.glob(
            join_paths(input_dir,
                       'c??.??',
                       str(form_idx),
                       'bfactors.csv')))


        data_dict = {
            float(filename.split('/')[-3].replace('c', '')): pd.read_csv(filename, header=0, index_col='residue_number')
                for filename in file_paths}

        # Combine all dataframes into one
        # and add another index
        bfactors_cut = pd.concat(data_dict,
            ignore_index=False)
        bfactors_cut.index = bfactors_cut.index.set_names(['dist_cutoff','residue_number'])

        bfactors_dc[form_idx] = bfactors_cut

    # Concatenate with form index
    bfactors_dc = pd.concat(bfactors_dc,
            ignore_index=False)
    bfactors_dc.index = bfactors_dc.index.set_names(['form_idx','dist_cutoff','residue_number'])

    # Save eigenvalues
    bfactors_dc.to_csv(join_paths(
        output_dir, 'dc.bfactors.csv'), float_format='%.2f')

def clean_benms():
    """ Clean BENMs.
    """
    # Extract eigenvalues
    subprocess.run(['bash', 'src/extract_benm.sh'])

    config = utils.read_config()
    input_dir = join_paths(config['data']['intFilePath'], 'scan-benm')
    output_dir = config['data']['intFilePath']

    eigvlas_benms = {}
    for form_idx in range(3):
        file_paths = sorted(glob.glob(
            join_paths(input_dir,
                       'b????',
                       str(form_idx),
                       'eigvals.csv')))

        data_dict = {
            float(filename.split('/')[-3].replace('b', '')): np.loadtxt(filename)
            for filename in file_paths}

        eigvlas_benms[form_idx] = clean_dict(data_dict)

    # Save eigenvalues
    eigvlas_benms[0].to_csv(join_paths(
        output_dir, 'benm.eigvals.00.csv'), float_format='%.6f')
    eigvlas_benms[1].to_csv(join_paths(
        output_dir, 'benm.eigvals.01.csv'), float_format='%.6f')
    eigvlas_benms[2].to_csv(join_paths(
        output_dir, 'benm.eigvals.02.csv'), float_format='%.6f')


def main():
    """ Master script which takes all raw series, cleans them,
        and outputs to a flat file.
    """
    clean_exp()
    clean_md()
    clean_dc()
    clean_benms()

    return None


if __name__ == "__main__":
    # Master function
    main()
