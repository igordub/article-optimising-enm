from turtle import color
import numpy as np
import pandas as pd

from os.path import join as join_paths, basename as get_basename
import glob
import subprocess

from cycler import cycler
from matplotlib import pyplot as plt

import utilities as utils

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


def hist_df(df, bins=range(0, 1000, 25), axis=None, 
    colour_map=[], linestyle_wheel=plt.cycler("marker",['-'])):
    """ Plot histogram for Dataframe 
        and label columns.
    """
    if axis == None:
        fig, ax = plt.subplots()
    else:
        ax=axis

    for i, column in enumerate(df.columns):
        data = df[column]
        if len(colour_map) == 0:
            hist_kwargs = dict(
                histtype='step',
                alpha=1,
                density=False,
                label=column,
                bins=bins)
        else:
            hist_kwargs = dict(
                histtype='step',
                alpha=1,
                color = colour_map[i%len(colour_map)],
                density=False,
                label=column,
                bins=bins,
                orientation='horizontal')

        bin_vals, _, _ = ax.hist(data, **hist_kwargs)

    return ax


def compare_to_md(eigvals_enm, eigvals_md, output_path, legend_kwargs=None):
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
                _ = hist_df(df_enm, bins=bins, axis=ax)

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



            # ax.autoscale_view()

    # Labels
    axs[1][-1].set_ylabel('# modes per bin')
    axs[-1][-1].set_xlabel('$\omega$ ($cm^{-1}$)')


    fig.supxlabel('Mode',y=0.08, x=0.4)
    fig.supylabel('$\omega$ ($cm^{-1}$)')

    plt.tight_layout()
    plt.subplots_adjust(bottom=0.15)

    # Legend
    axs[-1][1].legend(**legend_kwargs)
    
    plt.savefig(output_path)
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
                # Chi-square test restults
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
                                data, **scatter_kwargs)
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
    axs[0][1].legend(['ENM', 'MD', 'X-ray'],
        loc='upper left', 
        fancybox=False, 
        shadow=False, 
        ncol=1,
        fontsize='medium',
        handlelength=1)

    axs[-1][0].legend(**legend_kwargs)
    
    plt.savefig(output_path)
    plt.close()

def plot_coop(coop_dc, coop_benm, dc_cm, benm_cm, output_path, lgnd_dc, lgnd_benm, dc_lbls=[]):
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

    
    for i, axs_row in enumerate(axs):
        for j, ax in enumerate(axs_row):
            # Plot
            if j == 0:
                # dc scan
                _ = lineplot_df(coop_dc[i], 
                    modes_plot=25, 
                    colour_map=dc_cm,
                    axis=ax)

                ax.set_xticks(np.arange(0, 26, 5))
                ax.locator_params(axis="y", nbins=5)

                # Non-cooperative value
                ax.axhline(y=1.0,
                   color='black',
                   linestyle=':',
                   linewidth=1)

            elif j == 1:
                # BENM scan
                _ = lineplot_df(coop_benm[i],
                    modes_plot=25, 
                    colour_map=benm_cm,
                    axis=ax)

                ax.set_xticks(np.arange(0, 26, 5))
                ax.locator_params(axis="y", nbins=5)

                # Non-cooperative value
                ax.axhline(y=1.0,
                   color='black',
                   linestyle=':',
                   linewidth=1)

            ax.autoscale_view()


    fig.supylabel('$K_2/K_1$')
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
    axs[-1][-1].legend(**lgnd_benm)

    plt.savefig(output_path)
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
    dist_cutoffs = np.arange(7.5,15.5,0.5)

    eigvals_dc_0 = [data['eigvals_dc_scaled'][0].iloc[:, :].drop(columns=['7.0'], errors='ignore')
        for data in data_list]

    chi2_data = [data['eigvals_dc_chi2'] for data in data_list]

    chi2_data = [chi2[chi2['sort_type']=='scatter'].pivot(
        index='modes',
        columns='variable',
        values='chi2').drop(columns=[7.0], errors='ignore') for chi2 in chi2_data]

    bfactors_data = [data['bfactors'] for data in data_list]


    coop_dc = [data['coop_dc'].drop(columns=['7.0'], errors='ignore') for data in data_list]
    coop_benm = [data['coop_benm'].drop(columns=['250.0','300.0'], errors='ignore') for data in data_list]

    # dc_cm = plt.cm.rainbow(np.linspace(0.0, 1, len(dist_cutoffs)))
    dc_cm = plt.cm.tab20(np.arange(len(dist_cutoffs)))

    rc_params = {'axes.prop_cycle':
            plt.cycler("color", dc_cm)}


    with plt.rc_context(rc_params):
        legend_kwargs = dict(
            loc='upper center', 
            bbox_to_anchor=(0.35,-0.4),
            fancybox=False, 
            shadow=False, 
            ncol=9,
            fontsize='medium')
        compare_to_md(eigvals_dc_0,
                    eigvals_md,
                    join_paths(config['data']['outPathScratch'], 'dc_scan.png'),
                    legend_kwargs=legend_kwargs)

    # Choose distance cutoff
    rc_params = {'axes.prop_cycle':plt.cycler("color", dc_cm),
        'lines.linewidth': 1.5}

    dc_lbls = ["$d_c = {:.1f} \: \AA{{}}$".format(dc) for dc in [8,8.5,8.5]]

    with plt.rc_context(rc_params):
        lgnd_dc = dict(
            loc='upper center', 
            bbox_to_anchor=(1,-0.45),
            fancybox=False, 
            shadow=False, 
            ncol=8,
            fontsize='medium')
        choose_dc(chi2_data,
            bfactors_data,
            join_paths(config['data']['outPathScratch'], 'choose_dc.png'),
            cb_colors,
            dc_lbls=dc_lbls,
            legend_kwargs=lgnd_dc)

    # BENM scan
    back_coeffs = np.array([1,10,20,30,40,50,100,150,200])
    eigvals_benm_0 = [data['eigvals_benm_scaled'][0].iloc[:, :].drop(columns=['250.0','300.0'], errors='ignore')
        for data in data_list]


    cm_min = 0.3
    cm_max = 0.8
    values_scaled = [cm_min + ((cm_max-cm_min)/(back_coeffs.max()-back_coeffs.min())) * (value - back_coeffs.min()) 
        for value in back_coeffs]
    values_scaled = np.array(values_scaled)
    benm_cm = plt.cm.Greys(values_scaled)
    rc_params = {'axes.prop_cycle':
            plt.cycler("color", benm_cm)}
    
    lgnd_benm = dict(
        loc='upper center', 
        bbox_to_anchor=(0.5,-0.3),
        fancybox=False, 
        shadow=False, 
        mode=None, 
        ncol=6,
        fontsize='medium')

    with plt.rc_context(rc_params):
        compare_to_md(eigvals_benm_0,
                    eigvals_md,
                    join_paths(config['data']['outPathScratch'], 'benm_scan.png'),
                    legend_kwargs = lgnd_benm)
    
    # Cooperativity
    lgnd_dc = dict(
        loc='upper center', 
        bbox_to_anchor=(0.35,-0.4),
        fancybox=False, 
        shadow=False, 
        ncol=6,
        fontsize='small')

    lgnd_benm = dict(
        loc='upper center', 
        bbox_to_anchor=(0.6,-0.4),
        fancybox=False, 
        shadow=False, 
        mode=None, 
        ncol=4,
        fontsize='small') 

    plot_coop(coop_dc,
        coop_benm, 
        dc_cm, 
        benm_cm,
        join_paths(config['data']['outPathScratch'], 'coop.png'),
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



def main():
    """ Master script.
    """
    plot()


if __name__ == "__main__":
    main()
