import os
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from os.path import join as join_paths


dist_inpath = 'dist.dat'
cross_inpath = 'crosscor.dat'
xlbl = 'Residue Number'
ylbl = 'Residue Number'
ttl = ''
suffix=''
output_dir='.'

#############################################################################
# Read arguments from terminal, and assign input files and a name that all output files will contain.
#############################################################################
for x in range(1, len(sys.argv)):
    if sys.argv[x] == '-dist':
        dist_inpath = sys.argv[x+1]

    if sys.argv[x] == '-cross':
        cross_inpath = sys.argv[x+1]

    if sys.argv[x] == '-xlabel':
        xlbl = sys.argv[x+1]

    if sys.argv[x] == '-ylabel':
        ylbl = sys.argv[x+1]

    if sys.argv[x] == '-title':
        ttl = sys.argv[x+1]

    if sys.argv[x] == '-suffix':
        suffix = sys.argv[x+1]

    if sys.argv[x] == '-outdir':
        output_dir = sys.argv[x+1]

    if sys.argv[x] == '-help':
        print('\n\nProgram to plot overlap data...\n\nOPTIONS:\n'
            '-i = Name of input file (Default=overlap.dat)\n'
            '-xlabel = Label for x axis (Default=mode i)\n'
            '-ylabel = Label for y axis (Default=mode j)\n'
            '-title = Title for plot\n')
        exit()


def plot(dist_inpath, cross_inpath, xlbl, ylbl, ttl):
    """ Plot distance cross-correlation plot.
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

    # ---------------------------------------------------------------------
    # Igor's suggestions
    # Set Seaborn
    sns.set()

    # Create masking matrix
    mask_dist = np.triu(df_dist.to_numpy())
    mask_cross = np.tril(df_cross.to_numpy())
    mask_div = np.tril(df_div.to_numpy())

    f = 1
    l = 1
    sns.set_style("ticks")
    sns.set(rc={'figure.figsize': (13, 8)})
    sns.set_context("paper", font_scale=f, rc={"lines.linewidth": l})

    nice_plot = {
        'font.size': 24,
        'xtick.labelsize': 22,
        'ytick.labelsize': 22,
        'figure.autolayout': False,
        'axes.titlesize': 22,
        'axes.labelsize': 24,
        'lines.linewidth': 0.5,
        'lines.markersize': 6,
        'legend.fontsize': 13
    }

    plt.rcParams.update(nice_plot)

    fig = plt.figure(1)
    ax = fig.add_subplot(111)

    ax1 = sns.heatmap(df_dist.to_numpy(),
                      vmin=0,
                      vmax=16,
                      mask=mask_dist,
                      square=True,
                      cmap=plt.cm.gist_yarg_r,
                      cbar_kws={'ticks': np.arange(0, 18, 2)}, cbar=True)

    column_maxes = df_div.max()
    div_max = column_maxes.max()
    column_mins = df_div.min()
    div_min = column_mins.min()    

    ax2 = sns.heatmap(df_div.to_numpy(),
                      vmin=div_min,
                      vmax=div_max,
                      mask=mask_div,
                      square=True,
                      cmap=plt.cm.RdBu_r,
                      cbar_kws={},
                      cbar=True)

    # Draw spines
    for _, spine in ax.spines.items():
        spine.set_visible(True)

    # ax1.set_xticks(residues)
    # ax1.set_yticks(residues)

    ax.locator_params(axis="x", nbins=6)
    ax.locator_params(axis="y", nbins=6)

    ax.set_xlabel(xlbl)
    ax.set_ylabel(ylbl)
    ax.set_title(ttl)
    plt.gca()

    # -------------------------------------------------------------------
    # Igor's suggestions

    # Invert y-axis
    plt.gca().invert_yaxis()
    plt.xticks(rotation=0)

    residue_total, _ = ax1.get_xlim()

    return fig, ax, df_dist, df_cross, df_div

print('Plotting...')
fig, ax, df_dist, df_cross, df_div= plot(dist_inpath, cross_inpath, xlbl, ylbl, ttl)

print('Saving...')
if suffix == "":
    # Save data
    df_dist.to_csv(join_paths(output_dir,'dist.csv'),
            float_format='%.2f')

    df_cross.to_csv(join_paths(output_dir,'cross.csv'),
            float_format='%.2f')

    df_div.to_csv(join_paths(output_dir,'div.csv'),
            float_format='%.2f')

    # Save plot
    plt.tight_layout()
    plt.savefig(join_paths(output_dir,'dist_div.png'))
    plt.close()
else:
    # Save data
    df_dist.to_csv(join_paths(output_dir,'dist.{}.csv'.format(suffix)),
            float_format='%.2f')

    df_cross.to_csv(join_paths(output_dir,'cross.{}.csv'.format(suffix)),
            float_format='%.2f')

    df_div.to_csv(join_paths(output_dir,'div.{}.csv'.format(suffix)),
            float_format='%.2f')

    # Save plot
    plt.tight_layout()
    plt.savefig(join_paths(output_dir,'dist_div.{}.png'.format(suffix)))
    plt.close()

print('DONE!')