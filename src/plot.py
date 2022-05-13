from turtle import color
import numpy as np
import pandas as pd

from os.path import join as join_paths, basename as get_basename
import glob
import subprocess

from cycler import cycler
from matplotlib import pyplot as plt

import utilities as utils

def select_df_list(df, lst):
    """ Select dataframe colums based on list values.
    """
    
    df_out = df[df.columns.intersection(lst)]
    return df_out

def cm_to_inch(cm_value):
    """ Convert cm to incehs.
    """
    inch_value = cm_value / 2.54

    return inch_value

def load_data(input_dir):
    """ Loads processed data for a protein complex.
    """
    data_dic={}
    
    # aaMD
    data_dic['eigvals_md'] = pd.read_csv(
        join_paths(input_dir, 'md.eigvals.csv'),
        header=0,
        index_col='mode')

    # dc scan
    # Eigenvalues
    data_dic['eigvals_dc'] = {}
    for form_idx in range(3):
        data_dic['eigvals_dc'][form_idx] = pd.read_csv(
            join_paths(input_dir, 'dc.eigvals.0{}.csv'.format(form_idx)),
            header=0,
            index_col='mode')

    data_dic['eigvals_dc_scaled'] = {}
    for form_idx in range(1):
        data_dic['eigvals_dc_scaled'][form_idx] = pd.read_csv(
            join_paths(input_dir, 'dc.eigvals.scaled.0{}.csv'.format(form_idx)),
            header=0,
            index_col='mode')

    # Allostery
    data_dic['kd_dc_1'] = pd.read_csv(join_paths(input_dir, 'dc.kd.01.csv'),
                        header=0,
                        index_col='mode')
    data_dic['kd_dc_2'] = pd.read_csv(join_paths(input_dir, 'dc.kd.02.csv'),
                        header=0,
                        index_col='mode')
    data_dic['coop_dc'] = pd.read_csv(join_paths(input_dir, 'dc.coop.csv'),
                        header=0,
                        index_col='mode')

    # Chi-square data
    data_dic['eigvals_dc_chi2'] = pd.read_csv(join_paths(input_dir, 'dc.eigvals.chi2.csv'),
                            header=0,
                            index_col=None)

    # B-factors
    data_dic['bfactors'] = pd.read_csv(join_paths(input_dir, 'bfactors.csv'),
                            header=0,
                            index_col='residue_number')

    # BENM scan
    # Eigenvalues
    data_dic['eigvals_benm'] = {}
    for form_idx in range(3):
        data_dic['eigvals_benm'][form_idx] = pd.read_csv(
            join_paths(input_dir, 'benm.eigvals.0{}.csv'.format(form_idx)),
            header=0,
            index_col='mode')

    data_dic['eigvals_benm_scaled'] = {}
    for form_idx in range(1):
        data_dic['eigvals_benm_scaled'][form_idx] = pd.read_csv(
            join_paths(input_dir, 'benm.eigvals.scaled.0{}.csv'.format(form_idx)),
            header=0,
            index_col='mode')

    # Allostery
    data_dic['kd_benm_1'] = pd.read_csv(join_paths(input_dir, 'benm.kd.01.csv'),
                        header=0,
                        index_col='mode')
    data_dic['kd_benm_2'] = pd.read_csv(join_paths(input_dir, 'benm.kd.02.csv'),
                        header=0,
                        index_col='mode')
    data_dic['coop_benm'] = pd.read_csv(join_paths(input_dir, 'benm.coop.csv'),
                        header=0,
                        index_col='mode')

    # Chi-square data
    data_dic['eigvals_benm_chi2'] = pd.read_csv(join_paths(input_dir, 'benm.eigvals.chi2.csv'),
                            header=0,
                            index_col=None)
   
    return data_dic

def scatter_df(df, modes_plot=25, axis=None, 
    colour_map=[], marker_wheel=plt.cycler("marker",['o'])):
    """ Plot scatter plot for Dataframe
        and label columns.
    """
    if axis == None:
        fig, ax = plt.subplots()
    else:
        ax=axis

    for i, column in enumerate(df.columns):
        data = df[column]
        if len(colour_map) == 0:    
            scatter_kwargs = dict(
                alpha=1,
                label=column)
        else:
            scatter_kwargs = dict(
                alpha=1,
                color = colour_map[i%len(colour_map)],
                label=column)

        ax.plot(data.index[:modes_plot],
                   data.iloc[:modes_plot], **scatter_kwargs)

    return ax


def lineplot_df(df, modes_plot=25, axis=None, 
    colour_map=[], linestyle_wheel=plt.cycler("marker",['-'])):
    """ Plot lineplot for Dataframe 
        and label columns.
    """
    if axis == None:
        fig, ax = plt.subplots()
    else:
        ax=axis

    for i, column in enumerate(df.columns):
        data = df[column]
        if len(colour_map) == 0:    
            lineplot_kwargs = dict(
                alpha=1,
                label=column)
        else:
            lineplot_kwargs = dict(
                alpha=1,
                color = colour_map[i%len(colour_map)],
                label=column)

        ax.plot(data.index[:modes_plot],
                data.iloc[:modes_plot], **lineplot_kwargs)

    return ax


def hist_df(df, bins=range(0, 1000, 25), axis=None, alpha_list=[]):
    """ Plot histogram for Dataframe 
        and label columns.
    """
    if axis == None:
        fig, ax = plt.subplots()
    else:
        ax=axis

    for i, column in enumerate(df.columns):
        data = df[column]
        if len(alpha_list) == 0:
            hist_kwargs = dict(
                histtype='step',
                alpha=1,
                density=False,
                label=column,
                bins=bins)
        else:
            hist_kwargs = dict(
                histtype='step',
                alpha=alpha_list[i],
                density=False,
                label=column,
                bins=bins)

        bin_vals, _, _ = ax.hist(data, **hist_kwargs)

    return ax


def compare_to_md(eigvals_enm, eigvals_md, output_path, legend_kwargs=None, alpha_lists=[[]]):
    """ Plots apo form eigenvalues for differnt ENMs
        and all-atom MD eigenvalues: scatter plot, line plot
        and histogram.
    """
    subplot_width = 6.1 # cm
    subplot_width = cm_to_inch(subplot_width)
    subplot_height = subplot_width 

    fig, axs = plt.subplots(3,3, 
        figsize=(3*subplot_width, 3*subplot_height),
        constrained_layout=False,
        tight_layout=False)

    for i, axs_row in enumerate(axs):
        for j, ax in enumerate(axs_row):
            df_enm = eigvals_enm[i]
            df_md =  eigvals_md[i]
            alpha_list = alpha_lists[i]

            bin_width = 15
            bins = np.arange(0, df_md.max()+bin_width, bin_width)

            if df_enm.shape[0] != df_md.shape[0]:
                print('Check number of modes for ENM and MD.')
                return None

            # Plot
            if j == 0:
                # Low-frequency modes
                # ENM
                _ = lineplot_df(df_enm, modes_plot=25, axis=ax)

                # MD
                scatter_kwargs = dict(
                    alpha=1,
                    color='r',
                    label='MD')

                ax.plot(df_md.index[:25],
                        df_md[:25],
                        **scatter_kwargs)
                # Add ticks
                ax.set_xticks(np.arange(0, 26, 5))

                _, y_max = ax.get_ylim()
                ax.set_ylim((0, y_max))
                ax.set_yticks(np.arange(0, y_max+2.5, 2.5))

                ax.locator_params(axis="y", nbins=6)

            elif j == 1:
                # All modes
                # ENM
                modes_total = df_enm.shape[0]
                _ = lineplot_df(df_enm, modes_plot=modes_total, axis=ax)

                scatter_kwargs = dict(
                    alpha=1,
                    color='r',
                    linestyle = '-',
                    label='MD')
                ax.plot(df_md.index[:modes_total], df_md[:modes_total], **scatter_kwargs)
                # Add ticks
                x_max = df_md.index.max()
                ax.set_xticks(np.arange(0, x_max, 500))

                _,y_max = ax.get_ylim()
                ax.set_yticks(np.arange(0, y_max+200, 200))

                ax.locator_params(axis="y", nbins=8)

            else:
                # Density of states
                # Choose certain columns
                _ = hist_df(df_enm,
                    bins=bins, 
                    axis=ax,
                    alpha_list=alpha_list)

                hist_kwargs = dict(histtype='step',
                                alpha=1,
                                color='r',
                                density=False,
                                label='MD',
                                bins=bins)

                _, _, _ = ax.hist(df_md, **hist_kwargs)

                # Add ticks
                _, x_max = ax.get_xlim()
                ax.set_xticks(np.arange(0, x_max+250, 250))
                ax.locator_params(axis="x", nbins=6)

                _,y_max = ax.get_ylim()
                ax.set_yticks(np.arange(0, y_max+50, 50))
                ax.locator_params(axis="y", nbins=8)

    # Labels
    axs[1][-1].set_ylabel('# modes per bin')
    axs[-1][-1].set_xlabel('$\omega$ ($cm^{-1}$)')


    fig.supxlabel('Mode',y=0.08, x=0.4)
    fig.supylabel('$\omega$ ($cm^{-1}$)')

    plt.tight_layout()
    plt.subplots_adjust(bottom=0.15)

    # Legend
    axs[-1][1].legend(**legend_kwargs)

    # SUBPLOT LABELS
    text_kwargs = dict(
        ha='left', 
        va='top', 
        fontsize=18,
        fontfamily='sans-serif',
        fontweight='bold', 
        color='black')

    fig.text(0.0, 1, 'a', **text_kwargs)
    fig.text(0.0, 0.71, 'b', **text_kwargs)
    fig.text(0.0, 0.41, 'c', **text_kwargs)
    
    for extension in ['tiff','pdf','eps','png']:
        plt.savefig(output_path+'.'+extension)
    plt.close()


def choose_dc(chi2_data, bfactors_data, output_path, cb_colors, dc_lbls=[],legend_kwargs=None):
    """ Plots apo form eigenvalues for differnt ENMs
        and all-atom MD eigenvalues: scatter plot, line plot
        and histogram.
    """
    subplot_width = 6.1 # cm
    subplot_width = cm_to_inch(subplot_width)
    subplot_height = subplot_width 

    fig, axs = plt.subplots(3,2, 
        figsize=(3*subplot_width, 2*subplot_height),
        constrained_layout=False,
        tight_layout=False)

    
    for i, axs_row in enumerate(axs):
        for j, ax in enumerate(axs_row):
            chi2 = chi2_data[i]
            bfactors = bfactors_data[i]

            # Plot
            if j == 0:
                # Chi-square test results
                _ = lineplot_df(chi2, modes_plot=25, axis=ax)

                ax.set_xticks(np.arange(0, 26, 5))
                ax.locator_params(axis="y", nbins=5)

            elif j == 1:
                # B-factors                    
                for column in bfactors.columns:
                    data= bfactors[column]
                    residue_numbers = bfactors.index

                    if '_enm' in column:
                        lbl ='ENM'
                        color = cb_colors['blue']
                    elif '_md' in column:
                        lbl ='MD',
                        color=cb_colors['vermillion']
                    else:
                        lbl ='X-ray',
                        color=cb_colors['bluish_green']

                    scatter_kwargs = dict(
                        alpha=1,
                        label=lbl,
                        color=color)
                    
                    ax.plot(residue_numbers,
                                data, linewidth=1,**scatter_kwargs)
                ax.locator_params(axis="y", nbins=5)

            ax.autoscale_view()


    # B-factor labels
    if dc_lbls != [] and len(dc_lbls) != 0:
        for i, axs_row in enumerate(axs):
            axs_row[1].text(0.9, 0.975, dc_lbls[i],
                horizontalalignment='center', verticalalignment='center',
                transform=axs_row[1].transAxes)

    axs[1][0].set_ylabel('$\chi^2$ ($cm^{-1}$)', labelpad=5, fontsize='large')
    axs[1][1].set_ylabel('B-factor ($\AA{}^2$)', labelpad=5)

    axs[-1][0].set_xlabel('Modes')
    axs[-1][-1].set_xlabel('Residue number')

    plt.tight_layout()
    plt.subplots_adjust(bottom=0.2)

    # Legend
    axs[-1][-1].legend(['ENM', 'MD', 'X-ray'],
        loc='upper center',
        bbox_to_anchor=(0.5,-0.45), 
        fancybox=False, 
        shadow=False, 
        ncol=3,
        fontsize='medium',
        handlelength=2)

    axs[-1][0].legend(**legend_kwargs)



    # SUBPLOT LABELS
    text_kwargs = dict(
        ha='left', 
        va='top', 
        fontsize=18,
        fontfamily='sans-serif',
        fontweight='bold', 
        color='black')

    fig.text(0.0, 1, 'a', **text_kwargs)
    fig.text(0.0, 0.72, 'b', **text_kwargs)
    fig.text(0.0, 0.43, 'c', **text_kwargs)
    for extension in ['tiff','pdf','eps','png']:
        plt.savefig(output_path+'.'+extension)
    plt.close()

def plot_coop(coop_dc, coop_benm, dc_cycler, benm_cycler, output_path, lgnd_dc, lgnd_benm, dc_lbls=[]):
    """ Plot allostery data.
    """
    subplot_width = 6.1 # cm
    subplot_width = cm_to_inch(subplot_width)
    subplot_height = subplot_width 

    fig, axs = plt.subplots(3,2, 
        figsize=(3*subplot_width, 2*subplot_height),
        constrained_layout=False,
        tight_layout=False,
        sharex=True,
        sharey=True)

    # SET CYCLER
    for i, axs_row in enumerate(axs):
        for j, ax in enumerate(axs_row):
            if j == 0:
                # dc scan
                axs[i][j].set_prop_cycle(dc_cycler)
            elif j == 1:
                # BENM scan
                axs[i][j].set_prop_cycle(benm_cycler)
    # PLOT
    for i, axs_row in enumerate(axs):
        for j, ax in enumerate(axs_row):
            if j == 0:
                # dc scan
                _ = lineplot_df(coop_dc[i], 
                    modes_plot=100, 
                    axis=ax)

                
                ax.set_xticks(np.arange(0, 101, 25))
                ax.locator_params(axis="y", nbins=5)

                # Non-cooperative value
                ax.axhline(y=1.0,
                   color='black',
                   linestyle=':',
                   linewidth=1)

            elif j == 1:
                # BENM scan
                _ = lineplot_df(coop_benm[i],
                    modes_plot=100, 
                    axis=ax)

                ax.set_prop_cycle(benm_cycler)
                ax.set_xticks(np.arange(0, 101, 25))
                ax.locator_params(axis="y", nbins=5)

                # Non-cooperative value
                ax.axhline(y=1.0,
                   color='black',
                   linestyle=':',
                   linewidth=1)

            ax.autoscale_view()


    fig.supylabel('$K_2/K_1$', x=0.03)
    fig.supxlabel('Modes', y=0.1, x=0.55)

    plt.tight_layout()
    plt.subplots_adjust(bottom=0.2)

    # dc labels
    if dc_lbls != [] and len(dc_lbls) != 0:
        for i, axs_row in enumerate(axs):
            axs_row[-1].text(0.9, 0.975, dc_lbls[i],
                horizontalalignment='center', verticalalignment='center',
                transform=axs_row[1].transAxes)

    # Legend
    axs[-1][0].legend(**lgnd_dc)
    axs[-1][0].text(0.5,-0.875, '$d_c$', fontsize=14,
        va='bottom',
        transform=axs[-1][0].transAxes)
    axs[-1][-1].legend(**lgnd_benm)
    axs[-1][-1].text(0.5,-0.875, '$\epsilon$',fontsize=14,
        va='bottom',
        transform=axs[-1][-1].transAxes)

    # Subplot labels
    text_kwargs = dict(
        ha='left', 
        va='top', 
        fontsize=18,
        fontfamily='sans-serif',
        fontweight='bold', 
        color='black')

    fig.text(0.0, 1, 'a', **text_kwargs)
    fig.text(0.0, 0.725, 'b', **text_kwargs)
    fig.text(0.0, 0.45, 'c', **text_kwargs)


    for extension in ['tiff','pdf','eps','png']:
        plt.savefig(output_path+'.'+extension)
    plt.close()
    #################


def plot():
    config = utils.read_config()
    plt.style.use(config['viz']['default'])


    cb_colors = {key : "#{:02x}{:02x}{:02x}".format(*[int(x) for x in value.split(',')]) for key, value in config['colors'].items()}

    output_dir = config['data']['outPathScratch']

    # LOAD DATA
    data_cap = load_data(config['data']['cap'])
    data_gst = load_data(config['data']['gst'])
    data_mpro = load_data(config['data']['mpro'])

    data_list = [data_cap, data_gst, data_mpro]

    # PLOTTING
    eigvals_md = [
        data['eigvals_md'][data['eigvals_md'] ['form_idx']==0]['eigval'] for data in data_list]

    # dc scan
    dist_cutoffs = np.arange(7.5,12.5,0.5)
    dist_cutoffs_str = [str(elem) for elem in dist_cutoffs]

    # Select only cutoffs from 7.5 to 12.0 A
    eigvals_dc_0 = [select_df_list(data['eigvals_dc_scaled'][0].iloc[:, :], dist_cutoffs_str)
        for data in data_list]

    chi2_data = [data['eigvals_dc_chi2'] for data in data_list]

    chi2_data = [chi2[chi2['sort_type']=='scatter'].pivot(
        index='modes',
        columns='variable',
        values='chi2') for chi2 in chi2_data]

    chi2_data = [select_df_list(data, dist_cutoffs) for data in chi2_data]


    bfactors_data = [data['bfactors'] for data in data_list]


    coop_dc = [select_df_list(data['coop_dc'].iloc[:, :], dist_cutoffs_str) for data in data_list]
    coop_benm = [data['coop_benm'].drop(columns=['250.0','300.0'], errors='ignore') for data in data_list]
    

    # Color map for distance scan
    # dc_cm = plt.cm.rainbow(np.linspace(0.0, 1, len(dist_cutoffs)))
    dc_cm = plt.cm.tab20(np.arange(len(dist_cutoffs)))

    dc_cycler = (plt.cycler("color", dc_cm) +
                  plt.cycler("linestyle", ['-', '--']*5))

    rc_params = {'axes.prop_cycle':
            dc_cycler}

    alphas = np.ones(len(dist_cutoffs))

    alpha_lists = [alphas,alphas,alphas]

    # COMPARE TO MD
    with plt.rc_context(rc_params):
        legend_kwargs = dict(
            loc='upper center', 
            bbox_to_anchor=(0.35,-0.3),
            fancybox=False, 
            shadow=False, 
            ncol=6,
            fontsize='large')
        compare_to_md(eigvals_dc_0,
                    eigvals_md,
                    join_paths(config['data']['outPathScratch'], 'dc_scan'),
                    legend_kwargs=legend_kwargs,
                    alpha_lists=alpha_lists)

    # Choose distance cutoff
    rc_params = {'axes.prop_cycle':dc_cycler,
        'lines.linewidth': 1.5}

    dc_lbls = ["$d_c = {:.1f} \: \AA{{}}$".format(dc) for dc in [8,8.5,8.5]]

    # CHOOSE DC
    with plt.rc_context(rc_params):
        lgnd_dc = dict(
            loc='upper center', 
            bbox_to_anchor=(0.5,-0.45),
            fancybox=False, 
            shadow=False, 
            ncol=5,
            fontsize='medium')
        choose_dc(chi2_data,
            bfactors_data,
            join_paths(config['data']['outPathScratch'], 'choose_dc'),
            cb_colors,
            dc_lbls=dc_lbls,
            legend_kwargs=lgnd_dc)

    # BENM scan
    back_coeffs = np.array([1,5,10,20,30,40,50,100,150,200])
    eigvals_benm_0 = [data['eigvals_benm_scaled'][0].iloc[:, :].drop(columns=['250.0','300.0'], errors='ignore')
        for data in data_list]

    # BENM SCAN COLOURMAP
    cm_min = 0.3
    cm_max = 1
    values_scaled = [cm_min + ((cm_max-cm_min)/(back_coeffs.max()-back_coeffs.min())) * (value - back_coeffs.min()) 
        for value in back_coeffs]
    values_scaled = np.array(values_scaled)
    benm_cm = plt.cm.Greys(values_scaled)
    # benm_cm = plt.cm.turbo(np.arange(len(back_coeffs)))

    benm_cycler = (plt.cycler("color", benm_cm) +
                plt.cycler("linestyle", ['-', '--', ':','-', '--', ':','-', '--',':','-']))

    rc_params = {'axes.prop_cycle':
            benm_cycler}
    
    lgnd_benm = dict(
        loc='upper center', 
        bbox_to_anchor=(0.5,-0.3),
        fancybox=False, 
        shadow=False, 
        mode=None, 
        ncol=6,
        fontsize='large')
    
    alphas = np.zeros(len(back_coeffs))
    alpha_lists = [alphas.copy() for i in range(3)]
    np.put(alpha_lists[0], [0,7,9], [1,1])
    np.put(alpha_lists[1], [0,7,9], [1,1])
    np.put(alpha_lists[2], [0,7,9], [1,1])



    # COMPARE TO MD
    with plt.rc_context(rc_params):
        compare_to_md(eigvals_benm_0,
                    eigvals_md,
                    join_paths(config['data']['outPathScratch'], 'benm_scan'),
                    legend_kwargs = lgnd_benm,
                    alpha_lists=alpha_lists)
    
    # Cooperativity
    lgnd_dc = dict(
        loc='upper center', 
        bbox_to_anchor=(0.45,-0.4),
        fancybox=False, 
        shadow=False, 
        ncol=5,
        fontsize='small')

    lgnd_benm = dict(
        loc='upper center', 
        bbox_to_anchor=(0.5,-0.4),
        fancybox=False, 
        shadow=False, 
        mode=None, 
        ncol=5,
        fontsize='small') 

    plot_coop(coop_dc,
        coop_benm, 
        dc_cycler, 
        benm_cycler,
        join_paths(config['data']['outPathScratch'], 'coop'),
        lgnd_dc, 
        lgnd_benm, 
        dc_lbls=dc_lbls)

    # Chi-square test
    # Make Latex table
    # eigvals_dc_chi2[(eigvals_dc_chi2['sort_type'] == 'scatter') &
    #              (eigvals_dc_chi2['modes'] == 25)].loc[:, eigvals_dc_chi2.columns != 'sort_type'].to_latex(
    #     join_paths(config['data']['outPathScratch'],
    #                'dc.eigvals.chi2.tex'),
    #     index=False,
    #     caption="Chi-square test for aaMD and ENMs.",
    #     label="tab:chi2",
    #     position='center')

    # plot_chi2(eigvals_dc_chi2,
    #           output_dir,
    #           data_label='dc',
    #           modes_plot=25,
    #           ttl="{} | $d_c$ scan".format(protein_name))

    return None

def plot_stand_wave():
    """ Plots qualitaive plots of sinusondal homogenous
        and inhomogenous standing waves."""
    config = utils.read_config()
    plt.style.use(config['viz']['default'])

    fig, (ax1) = plt.subplots(1,1, 
        figsize=(cm_to_inch(6.1), cm_to_inch(3.05)),
        constrained_layout=False,
        tight_layout=True,
        sharex=True,
        sharey=True)


    x = np.linspace(0, 2*np.pi, 400)
    y1 = np.sin(x)
    y2 = np.sin(x + x**2 / (0.5*np.pi))

    # Constant wavelength
    ax1.plot(x, y1, lw=1, zorder=-1)
    # Varying wavelength
    ax1.plot(x, y2, lw=1, zorder=-1)

    # EN beads and springs
    ax1.hlines(0, 0, 2*np.pi, color='black', lw=1, zorder=-1)
    beads_xloc = np.linspace(0, 2, num=9) * np.pi
    ax1.scatter(beads_xloc, np.zeros(len(beads_xloc)),
        color='#aec7e8',
        edgecolors='black',
        s=50)

    ax1.get_xaxis().set_ticks([])
    ax1.get_yaxis().set_ticks([])

    for ax in [ax1]:
        for key in ['right', 'top', 'bottom', 'left']:
            ax.spines[key].set_visible(False)

    plt.xlim(-0.1*np.pi, 2.1*np.pi)

    # Annotation
    ax.annotate('EN bead', xy=(0, 0), xytext=(-0.1*np.pi, -0.5), size=6, 
        ha="center", va="center",
        arrowprops=dict(facecolor='black', 
            arrowstyle="-|>",
            connectionstyle="arc3"),
        )
    ax.annotate('EN spring', xy=(3/8 * np.pi,0.05), xytext=(2.5/8 * np.pi, -0.8), size=6, 
        ha="center", va="bottom",
        arrowprops=dict(facecolor='black', 
            arrowstyle="-|>",
            connectionstyle="arc3"),
        )

    ax.annotate('$\lambda_{fixed}$', 
        xy=(0.5 * np.pi, 1), 
        xytext=(0.55 * np.pi, 1),
        color='#bc80bd', 
        size=6, 
        ha="center", 
        va="bottom",
        arrowprops=dict(facecolor='black', 
            arrowstyle="-",
            lw=0),
        )

    ax.annotate('$\lambda_{varying}$', 
        xy=(1.5 * np.pi, 1), 
        xytext=(1.55 * np.pi, 1), 
        color='#fb8072',
        size=6, 
        ha="center", 
        va="bottom",
        arrowprops=dict(facecolor='black', 
            arrowstyle="-",
            lw=0),
        )
    for extension in ['tiff','pdf','eps','png']:
        plt.savefig(join_paths(config['data']['outPathScratch'], 'stand_waves.' + extension))

    return None


def main():
    """ Master script.
    """
    plot()
    plot_stand_wave()


if __name__ == "__main__":
    main()
