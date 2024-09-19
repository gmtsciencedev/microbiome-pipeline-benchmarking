#!/usr/bin/python3.6

from __future__ import division
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D
import os, sys

from radar_mpl import radar_factory
from final_figures import ctxt2, figwidth, plot_height

#--------------------------------------------------------------------------
#figures settings
plt.style.use('ggplot')
plt.style.use('seaborn-whitegrid')

cmap = mpl.cm.get_cmap("Paired")
marks = ["o", "s", "^", "d", "*", "+", "x", "."]
# cols = [cmap([j])[0] for j in range(10)]
cols = [cmap([j])[0] for j in [0, 1, 6, 7]]

cmap_default = plt.rcParams['axes.prop_cycle'].by_key()['color']
# cmap_default = np.array([mpl.colors.to_rgb(x) for x in cmap_default])
# cmap_default = np.c_[cmap_default, np.ones(cmap_default.shape[0])]

mpl.rcParams.update(ctxt2)

#--------------------------------------------------------------------------
#names and labels
names = ["kraken2 + bracken", "metaphlan3", "motus3", "metaphlan4"]
fnames = ["kraken2", "metaphlan3", "motus3", "metaphlan4"]
simuls = ["refKrak", "refMet4"]
spaces = ["gtdb207", "uhgg"]

simu_labs = {"refMet4": "refMet4", "refKrak": "refKrak"}
space_labs = {"gtdb207": "GTDB", "uhgg" : 'UHGG'}
tool_labs = {"kraken2 + bracken": "Krak/Brac", 
             "kraken2" : "Krak/Brac",
             "metaphlan3": "MPA3", 
             "metaphlan4": "MPA4", 
             "motus3": "mOTUs3", 
             "biomscope": "BiomS", 
             "simulation" : "Ref",
             "reference" : "Ref"}

labs = ["{}, {}".format(simu_labs[simu], space_labs[sp]) for simu in simuls for sp in spaces]
tlabs = [tool_labs[name] for name in names]
idx =pd.MultiIndex.from_product([simuls, spaces], names=["space", "simu"])

nn_smpls = np.arange(343)
nsats = len(names)

if __name__ =="__main__":
    """"
    Radar plots for comparing tools performance across conditions (simu, space)
    """
    plt.close("all")
    indir = "../analyses/combined_metrics/"

    figdir = "../figures_radar/"
    if not os.path.exists(figdir):
        os.makedirs(figdir)

    # meas = pd.read_csv(indir + "median_basic.tsv", sep="\t")
    # meds = pd.read_csv(indir + "median_measures.tsv", sep="\t", index_col=(0,1,2))

    #alpha-diversity
    richness = pd.read_csv(indir + "richness.tsv", sep="\t", index_col=(0,1,2)).swaplevel()
    shannon = pd.read_csv(indir + "shannon.tsv", sep="\t", index_col=(0,1,2)).swaplevel()
    
    #distance to the ground truth
    dist_to_ref = pd.read_csv(indir + "dist_to_ref.tsv", sep="\t", index_col=(0,1,2)).swaplevel()
    dist_to_ref_wu = pd.read_csv(indir + "dist_to_ref_wu.tsv", sep="\t", index_col=(0,1,2)).swaplevel()
    
    #presence metrics and FPRA
    superconf = pd.read_csv(indir + "superconf.tsv", sep="\t", index_col=[0,1,2,3]).swaplevel()
    TPR = superconf.xs("TPR", level=1)
    PPV = superconf.xs("PPV", level=1)
    FPRA = superconf.xs("FPRA", level=1)


    #--------------------------------------------------------------------------
    #attempting radar plot with matplotlib solution
    # dg1 = richness[names].div(richness["simulation"], axis=0).groupby(["simu", "space"]).median()
    dg2 = shannon[names].div(shannon["simulation"], axis=0).groupby(["simu", "space"]).median()

    dg1 = richness[names].sub(richness["simulation"], axis=0).\
    div(richness[names].add(richness["simulation"], axis=0)).abs().\
    groupby(["simu", "space"]).median()

    dg2 = shannon[names].sub(shannon["simulation"], axis=0).\
    div(shannon[names].add(shannon["simulation"], axis=0)).abs().\
    groupby(["simu", "space"]).median()

    dfs = [1. - dg1, 1. - dg2, 
           (1. - dist_to_ref).groupby(["simu", "space"]).median(), 
           (1. - dist_to_ref_wu).groupby(["simu", "space"]).median(), 
           TPR.groupby(["simu", "space"]).median(), 
           PPV.groupby(["simu", "space"]).median(), 
           (1. - FPRA).groupby(["simu", "space"]).min()]
    measures = ["Richness similarity", "Shannon similarity", "Bray-Curtis\nsimilarity",
                "UniFrac similarity", "Sensitivity", "Precision", "minimal TPRA"]

    # N = 6
    # theta = radar_factory(N, frame="polygon")

    colors = ['b', 'r', 'g', 'm', 'y']
    # fig, axs = plt.subplots(2, 4, figsize=(24, 12), subplot_kw=dict(projection='radar'))

    # for ax, df, meas in zip(axs.flatten(), dfs, measures):
    #     for name, c in zip(names, colors):
    #         d = df[name]
    #         ax.plot(theta, d, color=c)
    #         ax.fill(theta, d, facecolor=c, alpha=0.25, label='_nolegend_')
    #         ax.set_title(meas)
    #         ax.set_varlabels(labs)
    # legend = axs[0,0].legend(tlabs, loc=(0.9, .95),
    #                           labelspacing=0.1, fontsize='small')
    # axs[-1, -1].axis("off")
    # plt.savefig(figdir + "radar_by_measure.png")
    # plt.close(fig)


    df = pd.concat([df.stack() for df in dfs], axis=1)
    df.columns = measures

    N = 7
    theta = radar_factory(N, frame="polygon")

    # ratio=1.; top = 0.93; bottom = 0.13; left=0.08; hspace=0.35; wspace=0.4
    # Ly, lx, ly = plot_height(figwidth, 3, 2, ratio, bottom=bottom, top=top, left=left, hspace=hspace, wspace=wspace)
    ratio=1.; top = 0.8; bottom = 0.2; left=0.06; hspace=0.35; wspace=0.5
    Ly, lx, ly = plot_height(figwidth, 4, 1, ratio, bottom=bottom, top=top, left=left, hspace=hspace, wspace=wspace)
    Ly = np.round(Ly, 2)

    lls = ["-", "--", ":", "-.", "-"]
    fig, axs = plt.subplots(1, 4, figsize=(figwidth, Ly), subplot_kw=dict(projection='radar'))
    for ax, i, lab in zip(axs.flatten(order="F"), idx, labs):
        for name, c, ls in zip(names, cmap_default, lls):
            d = df.loc[i + (name,)]
            ax.plot(theta, d, color=c)
            # ax.fill(theta, d, facecolor=c, alpha=0.25, label='_nolegend_')
            ax.set_title(lab, position=(0.5, 1.1),
                horizontalalignment='center', verticalalignment='center')
            ax.set_varlabels(measures)
            ax.spines['polar'].set_visible(False)
    # legend = axs[0,0].legend(tlabs, loc=(0.9, .95),
    #                           labelspacing=0.1, fontsize='small')
    # handles, labels = axs[0,0].get_legend_handles_labels()
    # handles = [mpatches.Patch(color=c, label=lab) for c, lab in zip(cmap_default, tlabs)]
    handles = [Line2D([0], [0], color=c, label=lab) for c, lab in zip(cmap_default, tlabs)]
    fig.legend(handles, tlabs,\
     ncol=4, bbox_to_anchor=(0.5, 0.), loc="lower center")
    plt.subplots_adjust(bottom=bottom, top=top, left=left, hspace=hspace, wspace=wspace)
    plt.savefig(figdir + "FigS3_radar_by_conditions.png")
    plt.close(fig)


    ratio=1.; top = 0.8; bottom = 0.2; left=0.06; hspace=0.35; wspace=0.5
    Ly, lx, ly = plot_height(figwidth, nsats, 1, ratio, bottom=bottom, top=top, left=left, hspace=hspace, wspace=wspace)
    Ly = np.round(Ly, 2)

    lls = ["-", "-", "--", "--", ":", ":"]
    fig, axs = plt.subplots(1, nsats, figsize=(figwidth, Ly), subplot_kw=dict(projection='radar'))
    for ax, name in zip(axs.flatten(), names):
        for i, lab, c, ls in zip(idx, labs, cols, lls):
            d = df.loc[i + (name,)]
            ax.plot(theta, d, color=c, ls=ls)
            # ax.fill(theta, d, facecolor=c, alpha=0.25, label='_nolegend_')
            ax.set_title(tool_labs[name], position=(0.5, 1.2),
                horizontalalignment='center', verticalalignment='center')
            ax.set_varlabels(measures, fontsize=4)
            ticklabels = ax.get_xticklabels()
            # ticklabels[2].set_va("top")
            # ticklabels[-2].set_va("bottom")
            ax.spines['polar'].set_visible(False)
            ax.set_ylim([0., 1.])
            ax.set_yticks([0.2, 0.4, 0.6, 0.8, 1.])
            ax.set_yticklabels([str(x) for x in [0.2, 0.4, 0.6, 0.8, 1.]], fontsize=4, va="center")
            # ax.tick_params(axis='y', which='major', labelsize=4)

    # legend = axs[0,0].legend(labs, loc=(0.9, .95),
    #                           labelspacing=0.1, fontsize='small')
    # handles = [mpatches.Patch(color=c, label=lab) for c, lab in zip(cols, labs)]
    handles = [Line2D([0], [0], color=c, ls=ls, label=lab) for c, lab, ls in zip(cols, labs, lls)]
    fig.legend(handles, labs,\
     ncol=4, bbox_to_anchor=(0.5, 0.), loc="lower center")
    # axs[-1, -1].axis("off")
    plt.subplots_adjust(bottom=bottom, top=top, left=left, hspace=hspace, wspace=wspace)
    plt.savefig(figdir + "Fig8_radar_by_tool.png")
    plt.close(fig)

    sys.exit()
    #--------------------------------------------------------------------------
    #displaying tool rankings
    metrics = ["richness_diff", 
               "shannon_diff", 
               "braycurtis", 
               "unifrac", 
               "sensitivity", 
               "precision", 
               "fpra"]
    
    metrics_names = ["Richness",
                     "Shannon diversity",
                     "Bray-Curtis distance",
                     "UniFrac distance",
                     "Sensitivity",
                     "Precision",
                     "FPRA"]

    dfs = []
    for metric in metrics:
        fpath = "../analyses/analyses_stat_tests/rank_{}.tsv".format(metric)
        df = pd.read_csv(fpath, sep="\t", index_col=(0,1))
        df.columns = np.arange(1, 6)
        df = df.stack().reset_index(level=2).set_index(0, append=True)
        df.index.set_names("tool", level=2, inplace=True)
        df.columns = ["rank"]
        dfs.append(df)
        # break
        # df["simu"] = df["simu"].apply(lambda x: simu_labs[x])
        # df["space"] = df["space"].apply(lambda x: space_labs[x])
        # df.rename(columns={"simu": "Simulation", "space": "Common space"}, inplace=True)
        # df.fillna("", inplace=True)

    df = pd.concat(dfs, axis=1)
    df.columns = metrics_names


    fig, axs = plt.subplots(2, 3, figsize=(18, 12), subplot_kw=dict(projection='radar'))

    for ax, i, lab in zip(axs.flatten(order="F"), idx, labs):
        for name, c in zip(tlabs, colors):
            d = df.loc[i + (name,)]
            ax.plot(theta, 5-d, color=c)
            ax.fill(theta, 5-d, facecolor=c, alpha=0.25, label='_nolegend_')
            ax.set_title(lab)
            ax.set_varlabels(metrics_names)
            ax.set_yticks(np.arange(5))
            ax.set_yticklabels(np.arange(5, 0, -1))
            ax.spines['polar'].set_visible(False)
    legend = axs[0,0].legend(tlabs, loc=(0.9, .95),
                              labelspacing=0.1, fontsize='small')
    plt.savefig(figdir + "ranking_by_conditions.png")
    plt.close(fig)


    fig, axs = plt.subplots(2, 3, figsize=(18, 12), subplot_kw=dict(projection='radar'))

    for ax, name in zip(axs.flatten(), tlabs):
        for i, lab, c in zip(idx, labs, cols):
            d = df.loc[i + (name,)]
            ax.plot(theta, 5-d, color=c)
            ax.fill(theta, 5-d, facecolor=c, alpha=0.25, label='_nolegend_')
            ax.set_title(name)
            ax.set_varlabels(metrics_names)
            ax.set_yticks(np.arange(5))
            ax.set_yticklabels(np.arange(5, 0, -1))
            ax.spines['polar'].set_visible(False)
    legend = axs[0,0].legend(labs, loc=(0.9, .95),
                              labelspacing=0.1, fontsize='small')
    axs[-1, -1].axis("off")
    plt.savefig(figdir + "ranking_by_tool.png")
    plt.close(fig)
