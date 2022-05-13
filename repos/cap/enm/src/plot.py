import numpy as np
import pandas as pd

from os.path import join as join_paths, basename as get_basename
import glob
import subprocess

from cycler import cycler
from matplotlib import pyplot as plt

import utilities as utils

def viz_eigvecs():
    """ Vizualize first 50 modes of a chosen
        distance cutoff ANM.
    """
    subprocess.run(['pymol', '-qc', 'src/viz_eigvecs.py'], stdout=None)

    return None

def viz_benms():
    """ Vizualize BENMs.
    """
    subprocess.run(['pymol', '-qc', 'src/viz_benms.py'], stdout=None)

    return None

def scatter_df(df, modes_plot=25):
    """ Plot scatter plot for Dataframe
        and label columns.
    """
    fig, ax = plt.subplots()

    for column in df.columns:
        eigvals_enm = df[column]
        data = eigvals_enm

        scatter_kwargs = dict(
            alpha=0.6,
            # color = colourWheel[j%len(colourWheel)],
            # linestyle = lineStyles_hist[j%len(lineStyles_hist)],
            label=column)
        ax = ax
        ax.scatter(data.index[:modes_plot],
                   data.iloc[:modes_plot], **scatter_kwargs)

    return fig, ax


def lineplot_df(df, modes_plot=25):
    """ Plot lineplot for Dataframe 
        and label columns.
    """
    fig, ax = plt.subplots()

    for column in df.columns:
        data = df[column]

        lineplot_kwargs = dict(
            alpha=0.8,
            # color = colourWheel[j%len(colourWheel)],
            # linestyle = lineStyles_hist[j%len(lineStyles_hist)],
            label=column)
        ax = ax
        ax.plot(data.index[:modes_plot],
                data.iloc[:modes_plot], **lineplot_kwargs)

    return fig, ax


def hist_df(df, bins=range(0, 1000, 25)):
    """ Plot histogram for Dataframe 
        and label columns.
    """
    fig, ax = plt.subplots()

    for column in df.columns:
        data = df[column]

        hist_kwargs = dict(histtype='step',
                           alpha=0.8,
                           # color = colourWheel[j%len(colourWheel)],
                           linestyle='dashed',
                           density=False,
                           label=column,
                           bins=bins)
        ax = ax
        bin_vals, _, _ = ax.hist(data, **hist_kwargs)

    return fig, ax


def compare_to_md(df, eigvals_md, output_dir, data_label='tmp', ttl='', modes_plot = 25):
    """ Plots apo form eigenvalues for differnt ENMs
        and all-atom MD eigenvalues: scatter plot, line plot
        and histogram.
    """
    # Make number of modes equal
    modes_total = min(eigvals_md.shape[0], df.shape[0])
    df = df.iloc[:modes_total, :]
    eigvals_md = eigvals_md[:modes_total]

    # FIRST 100 MODES
    # SCATTER PLOT
    fig, ax = scatter_df(df, modes_plot=modes_plot)

    scatter_kwargs = dict(
        alpha=1,
        marker='x',
        color='r',
        # linestyle = lineStyles_hist[j%len(lineStyles_hist)],
        label='MD')
    data = eigvals_md
    ax.scatter(df.index[:modes_plot],
               data[:modes_plot],
               **scatter_kwargs)

    ax.set_xlabel("mode")
    ax.set_ylabel("$\lambda$ [$cm^{-1}$]")
    ax.set_title(ttl)
    ax.legend(frameon=False,
              ncol=7,
              fontsize='small',
              columnspacing=1,
              handlelength=1)

    plt.tight_layout()
    plt.savefig(join_paths(
        output_dir, '{}.eigvals.m{:04d}.png'.format(data_label, modes_plot)))
    plt.close()

    # ALL MODES
    # LINEPLOT
    modes_plot = modes_total

    fig, ax = lineplot_df(df, modes_plot=modes_plot)

    scatter_kwargs = dict(
        alpha=1,
        color='r',
        # linestyle = lineStyles_hist[j%len(lineStyles_hist)],
        label='MD')
    data = eigvals_md
    ax.plot(df.index[:modes_plot], data[:modes_plot], **scatter_kwargs)

    ax.set_xlabel("mode")
    ax.set_ylabel("$\lambda$ [$cm^{-1}$]")
    ax.set_title(ttl)
    ax.legend(frameon=False,
              ncol=7,
              fontsize='small',
              columnspacing=1,
              handlelength=1)

    plt.tight_layout()
    plt.savefig(join_paths(output_dir,
                '{}.eigvals.m{:04d}.png'.format(data_label, modes_total)))
    plt.close()

    # Density of states
    # HISTOGRAM
    fig, ax = plt.subplots()

    bin_width = 25
    bins = np.arange(0, eigvals_md.max()+bin_width, bin_width)

    fig, ax = hist_df(df, bins=bins)

    data = eigvals_md
    hist_kwargs = dict(histtype='step',
                       alpha=1,
                       color='r',
                       density=False,
                       label='MD',
                       bins=bins)

    bin_vals_md, _, _ = ax.hist(data, **hist_kwargs)

    ax.set_ylabel("# modes per bin")
    ax.set_xlabel("$\lambda$ [$cm^{-1}$]")
    ax.set_title(
        "{} | bin width = {} $cm^{{-1}}$".format(ttl, bin_width))
    ax.legend(frameon=False,
              ncol=7,
              fontsize='small',
              columnspacing=1,
              handlelength=1)

    plt.tight_layout()
    plt.savefig(join_paths(
        output_dir, '{}.dos.m{}.png'.format(data_label, modes_total)))
    plt.close()


def plot_benms_forms(df_0, df_1, df_2, output_dir):
    """ Plots apo, holo1 and holo2 forms' eigenvalues
        for all modes and all backbone stiffening coefficients.
    """
    # apo
    modes_total = df_0.shape[0]
    _, ax = lineplot_df(df_0, modes_plot=modes_total)

    ax.set_xlabel("mode #")
    ax.set_ylabel("$\lambda$ [$cm^{-1}$]")
    ax.set_title("apo | BENM | # modes = {} ".format(modes_total))
    ax.legend(frameon=False,
              ncol=7,
              fontsize='small',
              columnspacing=1,
              handlelength=1)

    plt.tight_layout()
    plt.savefig(join_paths(output_dir,
                'benm.eigvals.00.m{}.png'.format(modes_total)))
    plt.close()

    # holo1
    modes_plot = df_1.shape[0]
    _, ax = lineplot_df(df_1, modes_plot=modes_plot)

    ax.set_xlabel("mode #")
    ax.set_ylabel("$\lambda$ [$cm^{-1}$]")
    ax.set_title("holo1 | BENM | # modes = {} ".format(modes_plot))
    ax.legend(frameon=False,
              ncol=7,
              fontsize='small',
              columnspacing=1,
              handlelength=1)

    plt.tight_layout()
    plt.savefig(join_paths(output_dir,
                'benm.eigvals.01.m{}.png'.format(modes_plot)))
    plt.close()

    # holo2
    modes_plot = df_2.shape[0]
    _, ax = lineplot_df(df_2, modes_plot=modes_plot)

    ax.set_xlabel("mode #")
    ax.set_ylabel("$\lambda$ [$cm^{-1}$]")
    ax.set_title("holo2 | BENM | # modes = {} ".format(modes_plot))
    ax.legend(frameon=False,
              ncol=7,
              fontsize='small',
              columnspacing=1,
              handlelength=1)

    plt.tight_layout()
    plt.savefig(join_paths(output_dir,
                           'benm.eigvals.02.m{}.png'.format(modes_plot)))

    plt.close()

    return None


def plot_chi2(df, output_dir, data_label='tmp', modes_plot=25, ttl=''):
    # SCATTER
    # Convert long to wide format
    df_plot = df[df['sort_type'] == 'scatter'].pivot(
        index='modes',
        columns='variable',
        values='chi2')

    _, ax = lineplot_df(df_plot, modes_plot=modes_plot)
    ax.set_xlabel("modes")
    ax.set_ylabel("$\chi^2$")
    ax.set_title(ttl)
    ax.legend(frameon=False,
              ncol=7,
              fontsize='small',
              columnspacing=1,
              handlelength=1)

    plt.tight_layout()
    plt.savefig(join_paths(
        output_dir, '{}.eigvals.chi2.m{:04d}.png'.format(data_label, modes_plot)))
    plt.close()

    return None

def plot_bfactors(df, output_dir, data_label='tmp', ttl=""):
    """ Plot B-factors.
    """
    fig, ax = plt.subplots()

    for column in df.columns:
        bfactors= df[column]
        residue_numbers = df.index

        scatter_kwargs = dict(
            alpha=1,
            # color = colourWheel[j%len(colourWheel)],
            # linestyle = lineStyles_hist[j%len(lineStyles_hist)],
            label=column)
        ax.plot(residue_numbers,
                   bfactors, **scatter_kwargs)

    ax.set_xlabel("Residue number")
    ax.set_ylabel("B-factor [$\AA{}^2$]")
    ax.set_title(ttl)

    ax.legend(frameon=False,
                ncol=1,
                fontsize='small',
                columnspacing=1,
                handlelength=1)

    plt.tight_layout()
    plt.savefig(join_paths(
        output_dir, '{}.bfactors.png'.format(data_label)))
    plt.close()

def plot_allo(kd_1, kd_2, coop, output_dir, data_label='tmp', modes_plot=25, ttl=""):
    """ Plot allostery data.
    """
    # plt.rcParams["axes.prop_cycle"] = plt.cycler("color", plt.cm.Greys.colors)
    ################# K1 #################
    data = kd_1
    total_modes = coop.shape[0]
    for m in [modes_plot, total_modes]:
        if m > 100:
            fig, ax = lineplot_df(data, modes_plot=m)
        else:
            fig, ax = scatter_df(data, modes_plot=m)
        ax.set_xlabel("mode")
        ax.set_ylabel("$K_1$")
        ax.set_title(ttl)
        ax.legend(frameon=False,
                  ncol=7,
                  fontsize='small',
                  columnspacing=1,
                  handlelength=1)

        plt.tight_layout()
        plt.savefig(join_paths(
            output_dir, '{}.kd.01.m{:04d}.png'.format(data_label, m)))
        plt.close()

    ################# K2 #################
    data = kd_2
    for m in [modes_plot, total_modes]:
        if m > 100:
            fig, ax = lineplot_df(data, modes_plot=m)
        else:
            fig, ax = scatter_df(data, modes_plot=m)
        ax.set_xlabel("mode")
        ax.set_ylabel("$K_2$")
        ax.set_title(ttl)
        ax.legend(frameon=False,
                  ncol=7,
                  fontsize='small',
                  columnspacing=1,
                  handlelength=1)

        plt.tight_layout()
        plt.savefig(join_paths(
            output_dir, '{}.kd.02.m{:04d}.png'.format(data_label, m)))
        plt.close()

    ################# K2/K1 #################
    data = coop
    for m in [modes_plot, 100, 200, 500, 750, 800, total_modes]:
        fig, ax = lineplot_df(data, modes_plot=m)

        # Non-cooperative value
        ax.axhline(y=1.0,
                   color='black',
                   linestyle=':',
                   linewidth=1)

        ax.set_xlabel("mode")
        ax.set_ylabel("$K_2/K_1$")
        ax.set_title(ttl)
        ax.legend(frameon=False,
                  ncol=7,
                  fontsize='small',
                  columnspacing=1,
                  handlelength=1)

        plt.tight_layout()
        plt.savefig(join_paths(
            output_dir, '{}.coop.m{:04d}.png'.format(data_label, m)))
        plt.close()
    #################


def plot_dc(dist_cutoff, protein_name):
    config = utils.read_config()
    plt.style.use(config['viz']['default'])

    input_dir = config['data']['proFilePath']
    output_dir = config['data']['outPathScratch']

    # LOAD DATA
    # aaMD
    eigvals_md = pd.read_csv(
        join_paths(config['data']['proFilePath'], 'md.eigvals.csv'),
        header=0,
        index_col='mode')

    # Eigenvalues
    eigvals_dc = {}
    for form_idx in range(3):
        eigvals_dc[form_idx] = pd.read_csv(
            join_paths(input_dir, 'dc.eigvals.0{}.csv'.format(form_idx)),
            header=0,
            index_col='mode')

    eigvals_scaled = {}
    for form_idx in range(1):
        eigvals_scaled[form_idx] = pd.read_csv(
            join_paths(input_dir, 'dc.eigvals.scaled.0{}.csv'.format(form_idx)),
            header=0,
            index_col='mode')

    # Allostery
    kd_dc_1 = pd.read_csv(join_paths(input_dir, 'dc.kd.01.csv'),
                          header=0,
                          index_col='mode')
    kd_dc_2 = pd.read_csv(join_paths(input_dir, 'dc.kd.02.csv'),
                          header=0,
                          index_col='mode')
    coop_dc = pd.read_csv(join_paths(input_dir, 'dc.coop.csv'),
                          header=0,
                          index_col='mode')

    # Chi-square data
    eigvals_chi2 = pd.read_csv(join_paths(input_dir, 'dc.eigvals.chi2.csv'),
                               header=0,
                               index_col=None)

    # B-factors
    bfactors = pd.read_csv(join_paths(input_dir, 'bfactors.csv'),
                               header=0,
                               index_col='residue_number')

    # PLOTTING
    plot_bfactors(bfactors,
        output_dir,
        data_label='dc',
        ttl="{} | $d_c$ = {} $\AA{{}}$".format(protein_name, dist_cutoff))

    n = len(eigvals_dc[0].columns.unique())
    plt.rcParams["axes.prop_cycle"] = plt.cycler(
        "color", plt.cm.viridis(np.linspace(0.0, 1, n)))

    plot_allo(kd_dc_1.iloc[:, :],
              kd_dc_2.iloc[:, :],
              coop_dc.iloc[:, :],
              output_dir,
              data_label='dc',
              modes_plot=25,
              ttl="{} | $d_c$ scan".format(protein_name))

    compare_to_md(eigvals_scaled[0].iloc[:, :],
                  eigvals_md[eigvals_md['form_idx']==0]['eigval'],
                  output_dir,
                  data_label='dc',
                  ttl="{} | $d_c$ scan".format(protein_name),
                  modes_plot = 25)

    # Chi-square test
    # Make Latex table
    eigvals_chi2[(eigvals_chi2['sort_type'] == 'scatter') &
                 (eigvals_chi2['modes'] == 25)].loc[:, eigvals_chi2.columns != 'sort_type'].to_latex(
        join_paths(config['data']['outPathScratch'],
                   'dc.eigvals.chi2.tex'),
        index=False,
        caption="Chi-square test for aaMD and ENMs.",
        label="tab:chi2",
        position='center')

    plot_chi2(eigvals_chi2,
              output_dir,
              data_label='dc',
              modes_plot=25,
              ttl="{} | $d_c$ scan".format(protein_name))
    
    # Vizualize eigenvectors
    # viz_eigvecs()

    return None




def plot_benms(dist_cutoff, protein_name):
    config = utils.read_config()
    plt.style.use(config['viz']['default'])

    input_dir = config['data']['proFilePath']
    output_dir = config['data']['outPathScratch']

    # LOAD DATA
    # aaMD
    eigvals_md = pd.read_csv(
        join_paths(config['data']['proFilePath'], 'md.eigvals.csv'),
        header=0,
        index_col='mode')

    # Eigenvalues
    eigvals_benms = {}
    for form_idx in range(3):
        eigvals_benms[form_idx] = pd.read_csv(
            join_paths(input_dir, 'benm.eigvals.0{}.csv'.format(form_idx)),
            header=0,
            index_col='mode')

    eigvals_scaled = {}
    for form_idx in range(1):
        eigvals_scaled[form_idx] = pd.read_csv(
            join_paths(
                input_dir, 'benm.eigvals.scaled.0{}.csv'.format(form_idx)),
            header=0,
            index_col='mode')

    # Allostery
    kd_benms_1 = pd.read_csv(join_paths(input_dir, 'benm.kd.01.csv'),
                             header=0, index_col='mode')
    kd_benms_2 = pd.read_csv(join_paths(input_dir, 'benm.kd.02.csv'),
                             header=0, index_col='mode')
    coop_benms = pd.read_csv(join_paths(input_dir, 'benm.coop.csv'),
                             header=0, index_col='mode')

    # Chi-square data
    eigvals_chi2 = pd.read_csv(join_paths(input_dir, 'benm.eigvals.chi2.csv'),
                               header=0,
                               index_col=None)

    # PLOTTING
    n = len(eigvals_benms[0].columns.unique())
    plt.rcParams["axes.prop_cycle"] = plt.cycler(
        "color", plt.cm.plasma(np.linspace(0.0, 1, n)))

    plot_allo(kd_benms_1.iloc[:, :],
              kd_benms_2.iloc[:, :],
              coop_benms.iloc[:, :],
              output_dir,
              data_label='benm',
              modes_plot=25,
              ttl="{} | BENM scan".format(protein_name))

    compare_to_md(eigvals_scaled[0],
                  eigvals_md[eigvals_md['form_idx']==0]['eigval'],
                  output_dir,
                  data_label='benm',
                  ttl="{} | BENM scan".format(protein_name),
                  modes_plot = 25)

    # Chi-square test
    # Make Latex table
    eigvals_chi2[(eigvals_chi2['sort_type'] == 'hist') &
                 (eigvals_chi2['modes'] == eigvals_chi2['modes'].max())].loc[:, eigvals_chi2.columns != 'sort_type'].to_latex(
        join_paths(config['data']['outPathScratch'],
                   'benm.eigvals.chi2.tex'),
        index=False,
        caption="Chi-square test for aaMD and ENMs.",
        label="tab:chi2",
        position='center')

    # Vizualize BENMs
    # viz_benms()
    
    return None


def main(dist_cutoff, protein_name):
    """ Master script.
    """
    plot_dc(dist_cutoff, protein_name)
    plot_benms(dist_cutoff, protein_name)


if __name__ == "__main__":
    dist_cutoff = 8.0 
    protein_name = "CAP"
    main(dist_cutoff, protein_name)
