from email import header
import pytraj as pt

from matplotlib import pyplot as plt
import os
from os.path import join as join_paths, basename as get_basename
import pandas as pd
import numpy as np

import utilities as utils

def plot_lineplot(x, y, xlabel='', ylabel='', output_path='plot.png'):
    """ Plot and save a simple lineplot.
    """
    plt.plot(x,y)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)

    plt.tight_layout()

    plt.savefig(output_path)
    plt.close()

    return None

def plot_hist(data,bins=range(0,1000,25),output_path='hist.png'):
    """ Plot histogram with density of states.
    """

    _, ax = plt.subplots()

    hist_kwargs = dict(histtype='step',
                    alpha=1,
                    bins=bins)

    ax.hist(data, **hist_kwargs, label="MD")

    ax.set_ylabel("# modes per bin")
    ax.set_xlabel("$\lambda$ [$cm^{-1}$]")
    bin_width = int(bins[1] - bins[0])
    ax.set_title("# modes = {} | bin width = {} $cm^{{-1}}$".format(len(data), bin_width))

    ax.legend(frameon=False)

    plt.tight_layout()
    plt.savefig(output_path)
    plt.close()

    return None

def plot_traj(input_dir, output_dir):
    """ Plots analyzed trajectory results. 
    """
    # ANALYSIS
    # RMSD
    rmsd_data = pd.read_csv(join_paths(input_dir, 'rmsd.csv'), 
        header=0,
        index_col='frame_number')

    plot_lineplot(rmsd_data.index,rmsd_data['rmsd'], 
        output_path=join_paths(output_dir,'rmsd.png'),
        xlabel='Frame',
        ylabel='RMSD ($\AA{}$)')

    # RMSF
    rmsf_data = pd.read_csv(join_paths(output_dir, 'rmsf.csv'),
        header=0,
        index_col='atom_number')

    plot_lineplot(rmsf_data.index, rmsf_data['rmsf'], 
        output_path=join_paths(output_dir,'rmsf.png'),
        xlabel='Atom number',
        ylabel='RMSF ($\AA{}$)')


    # B-FACTORS
    bfactors_data = pd.read_csv(join_paths(output_dir, 'bfactors.csv'),
        header=0,
        index_col='residue_number')

    plot_lineplot(bfactors_data.index, bfactors_data['bfactor'], 
        output_path=join_paths(output_dir,'bfactors.png'),
        xlabel='Residue number',
        ylabel='B-factor ($\AA{}^2$)')


    # EIGENVALUES
    eigvals = pd.read_csv(join_paths(output_dir, 'eigvals.csv'),
        header=0,
        index_col='mode')

    plot_modes = 100

    _, [ax1, ax2] = plt.subplots(2, 1, constrained_layout=False)
    ax = ax1
    ax.scatter(eigvals.index[:plot_modes], eigvals['eigval'][:plot_modes], label='MD')
    ax.set_ylabel("$\lambda$ [$cm^{-1}$]")
    ax.legend(frameon=False)

    ax = ax2
    ax.plot(eigvals.index, eigvals['eigval'], label='MD')
    ax.set_ylabel("$\lambda$ [$cm^{-1}$]")
    ax.set_xlabel("mode")
    ax.legend(frameon=False)

    plt.tight_layout()

    plt.savefig(join_paths(output_dir, 'eigvals.png'))
    plt.close()

    bin_width = 25
    bins=np.arange(0, eigvals['eigval'].max()+bin_width, bin_width)
    plot_hist(eigvals['eigval'], bins=bins, output_path=join_paths(output_dir,'dos.png'))

    return None

def main():
    """ Master script.
    """
    input_dir = 'scratch'
    output_dir = input_dir

    # Custom settings for matplotlib
    config = utils.read_config()
    plt.style.use(config['viz']['default'])

    plot_traj(input_dir, output_dir)

if __name__ == "__main__":
    main()