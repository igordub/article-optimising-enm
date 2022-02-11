import numpy as np
import pandas as pd
from scipy.stats import chisquare
from os.path import join as join_paths
import subprocess

import glob
from matplotlib import pyplot as plt
import utilities as utils

def calc_chi2(obs, exp):
    """ Calculate chisquare sum between
        observed and expected 1D datasets.
    """
    if not (len(obs) == len(exp)):
        print("Array must match in length.")
        return None
    
    chi2 = (obs - exp)**2 / exp
    chi2_sum = chi2.sum()

    return chi2_sum


def project_eigvecs(dist_cutoff):
    """ Project eigenvectors on a PDB file.
    """
    subprocess.run(['bash', 'src/project_eigvecs.sh', '{}'.format(dist_cutoff)])

    return None

def calc_dist_cross(dist_cutoff, start_mode, end_mode):
    """ Calculate alpha-carbon - alpha-carbon distance matrix
        and cross-correlation for given mode interval.
    """
    subprocess.run(['bash', 'src/calc_dist_cross.sh', '{}'.format(dist_cutoff),
        start_mode,
        end_mode])


    return None

def analyze_dist_div(input_dir, txt_outpath, modes=[11,12,13], residues= [1,2]):
    data = {'m{:04d}'.format(mode_num):pd.read_csv(join_paths(input_dir,'div.m{:04d}.csv'.format(mode_num)), header=0, index_col=0) 
        for mode_num in modes}
    for _, value in data.items(): 
        value.columns = value.columns.astype(int)

    with open(txt_outpath, 'w') as file_txt:
        for key, value in data.items():
            file_txt.write(key)
            file_txt.write('\n')
            file_txt.write(value.loc[residues, residues].to_string())
            file_txt.write('\n')
            

def process_dist_cross(dist_inpath, cross_inpath, suffix=''):
    """ Process distances and cross-correlation data.
    """
    mi = []
    mi2 = []
    mj = []
    mj2 = []
    ol = []
    ol2 = []
    i = -1
    j = -1

    # Read distance data
    inlines = open(dist_inpath, 'r').readlines()

    if inlines[-1] == '\n':
        inlines[-1:] = []

    i = i+1
    mi.append([])
    mj.append([])
    ol.append([])

    for line in inlines:
        if line == '\n':
            i = i+1
            mi.append([])
            mj.append([])
            ol.append([])

        else:
            mi[i].append(int(line.split()[0]))
            mj[i].append(int(line.split()[1]))
            ol[i].append(float(line.split()[2]))

    mi = np.array(mi)
    mj = np.array(mj)
    ol = np.array(ol)

    # ------------------------------------------------------------------
    # Read cross-correlation data
    inlines = open(cross_inpath, 'r').readlines()

    if inlines[-1] == '\n':
        inlines[-1:] = []

    j = j+1
    mi2.append([])
    mj2.append([])
    ol2.append([])

    for line in inlines:
        if line == '\n':
            j = j+1
            mi2.append([])
            mj2.append([])
            ol2.append([])

        else:
            mi2[j].append(int(line.split()[0]))
            mj2[j].append(int(line.split()[1]))
            ol2[j].append(float(line.split()[2]))

    mi2 = np.array(mi2)
    mj2 = np.array(mj2)
    ol2 = np.array(ol2)

    residues = mi[:,0]

    # NumPy array to dataframe
    df_dist = pd.DataFrame(data=ol, index=residues, columns=residues)
    df_cross = pd.DataFrame(data=ol2, index=residues, columns=residues)
    
    # Calculate divergance
    column_maxes = df_dist.max()
    df_max = column_maxes.max()
    normalized_df = df_dist / df_max
    normalized_df = normalized_df.replace(0, np.NaN)
    df_div = df_cross / normalized_df

    # Save data
    if suffix == "":
        df_dist.to_csv('dist.csv',
                float_format='%.2f')

        df_cross.to_csv('cross.csv',
                float_format='%.2f')

        df_div.to_csv('div.csv',
                float_format='%.2f')
    else:
        df_dist.to_csv('dist.{}.csv'.format(suffix),
                float_format='%.2f')

        df_cross.to_csv('cross.{}.csv'.format(suffix),
                float_format='%.2f')

        df_div.to_csv('div.{}.csv'.format(suffix),
                float_format='%.2f')

def concat_bfactors(data_enm, data_md, data_exp):
    """ Concatenates ENM, aaMD and experimental
        B-factors into a single dataframe.
    """
    df_1 = data_enm['bfactor'].reset_index(drop=True)
    df_2 = data_md['bfactor'].reset_index(drop=True)
    df_3 = data_exp['bfactor'].reset_index(drop=True)

    df = pd.concat([df_1, df_2, df_3], axis=1)
    df['residue_number'] = data_enm['residue_number'].reset_index(drop=True)

    df.columns = ['bfactor_enm', 'bfactor_md', 'bfactor_exp', 'residue_number']
    df = df[['residue_number', 'bfactor_enm', 'bfactor_md', 'bfactor_exp']]

    return df


def chisquare_to_md(data_enm, data_md, bin_width=15):
    """ Perform chi-square test on raw and binned
        ENMs data to all-atom MD data.
    """
    # Dictionary to categorize and save data
    eigvals_chi2 = {'variable': [],
                    'sort_type': [],
                    'modes': [],
                    'chi2': []}

    for column in data_enm.columns:
        # Number of modes
        mode_list = np.arange(1,101)
        mode_list= np.append(mode_list, data_md.index.max())
        for m in mode_list:
            # SCATTER DATA
            chi2 = calc_chi2(data_enm[column].iloc[:m],
                data_md.iloc[:m])

            eigvals_chi2['variable'].append(column)
            eigvals_chi2['sort_type'].append('scatter')
            eigvals_chi2['modes'].append(m)
            eigvals_chi2['chi2'].append(chi2)

        # BINNED DATA
        for m in mode_list:
            bins = np.arange(0, data_md[:m].max(), bin_width)
            hist_kwargs = dict(density=False,
                               bins=bins)

            bin_vals_md, _, _ = plt.hist(data_md[:m],
                                         **hist_kwargs)
            bin_vals_enm, _, _ = plt.hist(data_enm[column][:m],
                                          **hist_kwargs)

            chi2 = calc_chi2(bin_vals_enm,
                bin_vals_md)

            eigvals_chi2['variable'].append(column)
            eigvals_chi2['sort_type'].append('hist')
            eigvals_chi2['modes'].append(m)
            eigvals_chi2['chi2'].append(chi2)

    eigvals_chi2 = pd.DataFrame(eigvals_chi2).astype(
        {'variable': float,
         'sort_type': str,
         'modes': int,
         'chi2': float})

    eigvals_chi2 = eigvals_chi2.sort_values(
        by=['sort_type', 'variable'])

    return eigvals_chi2


def scale_values(df, value_fit):
    """ Scales values in all dataframe
        columns  to fit a value in first row.
    """
    df_scaled = df.copy(deep=True)
    scaling_factors = value_fit / df_scaled.iloc[0, :]
    for column in df.columns:
        df_scaled[column] = df_scaled[column] * scaling_factors[column]

    return df_scaled, scaling_factors


def calc_coop(df_0, df_1, df_2):
    """ Calculates cooperativity for ENM apo, holo1
        and holo2 forms.
        Arrays or data frames must be of the same format
        for all forms.
    """
    # Calcualte dissociation constants
    # and cooperativity
    kd_1 = df_1 / df_0
    kd_2 = df_2 / df_1
    coop = kd_2 / kd_1

    # Remove missing values
    coop, kd_2, kd_1 = coop.dropna(), kd_2.dropna(), kd_1.dropna()

    # Calculate cumulative values
    coop, kd_2, kd_1 = np.cumprod(coop), np.cumprod(kd_2), np.cumprod(kd_2)

    return coop, kd_2, kd_1


def analyze_dc(dist_cutoff):
    # Get config file
    config = utils.read_config()

    # Load experimental data
    bfactors_exp = pd.read_csv(
        join_paths(config['data']['intFilePath'], 'exp.bfactors.csv'),
        header=0,
        index_col=None,
        dtype={'form_idx':int,'residue_numnber':int}).sort_values(['form_idx', 'residue_number'])

    # Load all-atom MD data
    eigvals_md = pd.read_csv(
        join_paths(config['data']['intFilePath'], 'md.eigvals.csv'),
        header=0,
        index_col=None,
        dtype={'form_idx':int,'mode':int}).sort_values(['form_idx', 'mode'])

    bfactors_md = pd.read_csv(
        join_paths(config['data']['intFilePath'], 'md.bfactors.csv'),
        header=0,
        index_col=None,
        dtype={'form_idx':int,'residue_numnber':int}).sort_values(['form_idx', 'residue_number'])


    eigvlas_dc = {}
    for form_idx in range(3):
        eigvlas_dc[form_idx] = pd.read_csv(
            join_paths(config['data']['intFilePath'],
                       'dc.eigvals.0{}.csv'.format(form_idx)),
            header=0,
            index_col='mode')


    bfactors_enm = pd.read_csv(
            join_paths(config['data']['intFilePath'],
                       'dc.bfactors.csv'.format(form_idx)),
            header=0,
            index_col=None,
            dtype={'form_idx':int,'dist_cutoff':float,'residue_numnber':int}).sort_values(['form_idx', 'dist_cutoff','residue_number'])

    # Filter ENMs with floppy modes
    df = eigvlas_dc[0]
    for column in df.columns:
        if df[column][1] > 1e-1:
            # No floppy modes
            pass
        else:
            # Floppy mode/s
            print("Distance cutoff of {} A has a floppy mode.".format(column))
            for form_idx in range(3):
                eigvlas_dc[form_idx] = eigvlas_dc[form_idx].drop(
                    columns=[column])

    for form_idx in range(3):
        eigvlas_dc[form_idx].to_csv(
            join_paths(config['data']['proFilePath'],
                       'dc.eigvals.0{}.csv'.format(form_idx)),
            float_format='%.6f')

    # Calculate cooperativity for BENM scan
    coop, kd_2, kd_1 = calc_coop(
        eigvlas_dc[0],
        eigvlas_dc[1],
        eigvlas_dc[2])

    # Save
    kd_1.to_csv(join_paths(
        config['data']['proFilePath'], 'dc.kd.01.csv'),
        float_format='%.6f')
    kd_2.to_csv(join_paths(
        config['data']['proFilePath'], 'dc.kd.02.csv'),
        float_format='%.6f')
    coop.to_csv(join_paths(
        config['data']['proFilePath'], 'dc.coop.csv'),
        float_format='%.6f')


    # Scale ENM eigenvalues to the first aaMD eigenvalue
    eigvals_scaled = {}
    value_fit = eigvals_md[(eigvals_md['form_idx']==0) & 
        (eigvals_md['mode']==1)]['eigval'].iloc[0]

    for form_idx in range(1):
        eigvals_scaled[form_idx], scaling_factors = scale_values(
            eigvlas_dc[form_idx],
            value_fit)

        print('Form index: {}'.format(form_idx))
        k_wt = 1
        k_sc = scaling_factors.pow(2) * k_wt
        print("Scaled spring constants [kcal/mol/A^2]")
        # round to two decimal places in python pandas
        pd.options.display.float_format = '{:.2f}'.format
        print(k_sc)

    for form_idx in range(1):
        eigvals_scaled[form_idx].to_csv(
            join_paths(config['data']['proFilePath'],
                       'dc.eigvals.scaled.0{}.csv'.format(form_idx)),
            float_format='%.6f')

    # Chi2-test
    eigvals_chi2 = chisquare_to_md(eigvals_scaled[0],
        eigvals_md[eigvals_md['form_idx']==0]['eigval'],
        bin_width=15)
    eigvals_chi2.to_csv(
        join_paths(config['data']['proFilePath'],
                   'dc.eigvals.chi2.csv'),
        index=None,
        float_format='%.2f')

    # B-factors
    bfactors = concat_bfactors(
        bfactors_enm[(bfactors_enm['form_idx']==0) & 
            (bfactors_enm['dist_cutoff']==dist_cutoff)],
        bfactors_md[bfactors_md['form_idx']==0],
        bfactors_exp[bfactors_exp['form_idx']==0])

    bfactors.to_csv(
        join_paths(config['data']['proFilePath'],
                   'bfactors.csv'),
        index=None,
        float_format='%.2f')

    # Copy aaMD data to output directory
    eigvals_md.to_csv(
        join_paths(config['data']['proFilePath'],
                   'md.eigvals.csv'),
        index=None,
        float_format='%.2f')

    bfactors_md.to_csv(
        join_paths(config['data']['proFilePath'],
                   'md.bfactors.csv'),
        index=None,
        float_format='%.2f')

    # Project eigenvectors for a specific cutoff
    # project_eigvecs(dist_cutoff)

def analyze_benms(dist_cutoff):
    # Get config file
    config = utils.read_config()

    # Load all-atom MD data
    eigvals_md = pd.read_csv(
        join_paths(config['data']['intFilePath'], 'md.eigvals.csv'),
        header=0,
        index_col=None,
        dtype={'form_idx':int,'mode':int})
    
    bfactors_md = pd.read_csv(
        join_paths(config['data']['intFilePath'], 'md.bfactors.csv'),
        header=0,
        index_col=None,
        dtype={'form_idx':int,'residue_numnber':int})

    eigvlas_benms = {}
    for form_idx in range(3):
        eigvlas_benms[form_idx] = pd.read_csv(
            join_paths(config['data']['intFilePath'],
                       'benm.eigvals.0{}.csv'.format(form_idx)),
            header=0,
            index_col='mode')

    # Save eigenvalues
    for form_idx in range(3):
        eigvlas_benms[form_idx].to_csv(
            join_paths(config['data']['proFilePath'],
                       'benm.eigvals.0{}.csv'.format(form_idx)),
            float_format='%.6f')

    # Calculate cooperativity for BENM scan
    coop_benm, kd_benms_2, kd_benms_1 = calc_coop(
        eigvlas_benms[0],
        eigvlas_benms[1],
        eigvlas_benms[2])

    # Save
    kd_benms_1.to_csv(
        join_paths(config['data']['proFilePath'], 'benm.kd.01.csv'),
        float_format='%.6f')
    kd_benms_2.to_csv(
        join_paths(config['data']['proFilePath'], 'benm.kd.02.csv'),
        float_format='%.6f')
    coop_benm.to_csv(join_paths(
        config['data']['proFilePath'], 'benm.coop.csv'),
        float_format='%.6f')

    # Chi-square test
    eigvals_scaled = {}
    value_fit = eigvals_md[(eigvals_md['form_idx']==0) & (eigvals_md['mode']==1)]['eigval'].iloc[0]
    for form_idx in range(1):
        eigvals_scaled[form_idx], scaling_factors = scale_values(
            eigvlas_benms[form_idx],
            value_fit)

        print('Form index: {}'.format(form_idx))
        k_wt = 1
        k_sc = scaling_factors.pow(2) * k_wt
        print("Scaled spring constants [kcal/mol/A^2]")
        # round to two decimal places in python pandas
        pd.options.display.float_format = '{:.2f}'.format
        print(k_sc)

    for form_idx in range(1):
        eigvals_scaled[form_idx].to_csv(
            join_paths(config['data']['proFilePath'],
                       'benm.eigvals.scaled.0{}.csv'.format(form_idx)),
            float_format='%.6f')

    eigvals_chi2 = chisquare_to_md(eigvals_scaled[0],
                                   eigvals_md[eigvals_md['form_idx']==0]['eigval'],
                                   bin_width=15)

    eigvals_chi2.to_csv(
        join_paths(config['data']['proFilePath'], 'benm.eigvals.chi2.csv'),
        index=None,
        float_format='%.2f')


def main(dist_cutoff):
    """ Master script which takes all raw series, cleans them,
        and outputs to a flat file.
    """

    analyze_dc(dist_cutoff)
    analyze_benms(dist_cutoff)

    return None


if __name__ == "__main__":
    config = utils.read_config()

    dist_cutoff = 8.5
    # main(dist_cutoff)

    analyze_dist_div(join_paths(config['data']['intFilePath'], 'dist_cross'),
        txt_outpath=join_paths(config['data']['outPathScratch'], 'div_res.txt'),
        modes=[9,10,11],
        residues=[492,493])
