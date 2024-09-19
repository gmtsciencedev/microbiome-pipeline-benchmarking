#!/usr/bin/python3.6

from __future__ import division
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import os, sys

pd.options.mode.chained_assignment = None  # default='warn'

#--------------------------------------------------------------------------
#figures settings
plt.style.use('ggplot')
plt.style.use('seaborn-whitegrid')

cmap = mpl.cm.get_cmap("Paired")
marks = ["o", "s", "^", "d", "*", "+", "x", "."]
cols = [cmap([j])[0] for j in range(10)]
cols1 = [cmap([j])[0] for j in [0, 1, 6, 7]]

#Nature asks figures fitting a page of 210 x 276 mm (8.27 by 10.86 inch)
figwidth = 8.27
ctxt2 = {"figure.dpi" : 600,
         "font.size" : 8,
         "legend.fontsize" : 8,
         "xtick.labelsize" : 6,
         "ytick.labelsize" : 6,
         "axes.titlesize" : 8,
         "axes.labelsize": 8,

         "figure.subplot.hspace": 0.05,
         "figure.subplot.wspace": 0.05,
         "figure.subplot.left" : 0.06,
         "figure.subplot.right" : 0.95,
         
         "lines.linewidth" : 1.,
         "lines.markersize" : 1,
         "grid.linewidth" : 0.5,

         "boxplot.flierprops.markersize" : 1,
         "boxplot.boxprops.linewidth": 0.5,
         "boxplot.capprops.linewidth": 0.5,
         "boxplot.medianprops.linewidth": 0.5,
         "boxplot.meanprops.linewidth": 0.5,
         "boxplot.whiskerprops.linewidth": 0.5}

mpl.rcParams.update(ctxt2)

def plot_height(Lx, nx, ny, ratio, left=None, bottom=None, right=None, top=None, wspace=None, hspace=None):
    """
    Calculate plot height, given its width and desired axis ratio.

    Parameters
    --------- 
    Lx : plot width
    nx : number of panels in x direction
    ny : number of panels in y direction
    ratio : panel height-to-width ratio
    left, bottom, right, top, wspace, hspace : plt.subplots_adjust parameters, taken from rcParams or context

    Returns
    -------
    Ly : plot hight
    lx : panel width
    ly : panel height
    """

    if left is None:
        left = mpl.rcParams["figure.subplot.left"]
    if right is None:
        right = mpl.rcParams["figure.subplot.right"]
    if bottom is None:
        bottom = mpl.rcParams["figure.subplot.bottom"]
    if top is None:
        top = mpl.rcParams["figure.subplot.top"]
    if hspace is None:
        hspace = mpl.rcParams["figure.subplot.hspace"]
    if wspace is None:
        wspace = mpl.rcParams["figure.subplot.wspace"]

    lx = Lx * (right - left) / (nx + (nx - 1) * wspace)
    ly = lx * ratio
    Ly = ly * (ny + (ny - 1) * hspace) / (top - bottom)

    return Ly, lx, ly

def simple_beeswarm(y, nbins=None, width=1.):
    """
    Returns x coordinates for the points in ``y``, so that plotting ``x`` and
    ``y`` results in a bee swarm plot.
    """
    y = np.asarray(y)
    if nbins is None:
        nbins = np.ceil(len(y) / 6).astype(int)

    # Get upper bounds of bins
    x = np.zeros(len(y))

    nn, ybins = np.histogram(y, bins=nbins)
    nmax = nn.max()

    #Divide indices into bins
    ibs = [np.nonzero((y>=ybins[0])*(y<=ybins[1]))[0]]
    for ymin, ymax in zip(ybins[1:-1], ybins[2:]):
        i = np.nonzero((y>ymin)*(y<=ymax))[0]
        ibs.append(i)

    # Assign x indices
    dx = 0.5 * width / (nmax // 2)
    for i in ibs:
        yy = y[i]
        if len(i) > 1:
            j = len(i) % 2
            i = i[np.argsort(yy)]
            a = i[j::2]
            b = i[j+1::2]
            x[a] = (0.5 + j / 3 + np.arange(len(a))) * dx
            x[b] = (0.5 + j / 3 + np.arange(len(b))) * -dx

    return x

def colored_box_swarm_plot(ax, data, idx, cols, ms=.5, levs=(1, 2)):
    for j, (i, c) in enumerate(zip(idx, cols)):
        if len(levs) == 1:
            yy = data.xs(i, level=levs[0]).to_numpy()
        else:
            yy = data.xs(i, level=levs).to_numpy()
        ax.boxplot(yy, widths=.8, positions=[j+1.], notch=True,\
                    showfliers=False,
                    showcaps=False,
                    boxprops=dict(color=c),
                    capprops=dict(color=c),
                    whiskerprops=dict(color=c),
                    flierprops=dict(color=c, markeredgecolor=c),
                    medianprops=dict(color=c))
        xx = simple_beeswarm(yy, width=0.4) + j + 1.
        ax.plot(xx, yy, "o", color=c, markersize=ms, alpha=0.5)

#--------------------------------------------------------------------------
#names and labels
names = ["kraken2 + bracken", "metaphlan3", "motus3", "metaphlan4"]
fnames = ["kraken2", "metaphlan3", "motus3", "metaphlan4"]
simuls = ["refKrak", "refMet4"]
spaces = ["gtdb207", "uhgg"]

simu_labs = {"refMet4": "refMet4", "refKrak": "refKrak"}
space_labs = {"gtdb207": "GTDB", "uhgg": "UHGG"}
tool_labs = {"kraken2 + bracken": "Krak/Brac", 
             "kraken2" : "Krak/Brac",
             "metaphlan3": "MPA3", 
             "metaphlan4": "MPA4", 
             "motus3": "mOTUs3", 
             "biomscope": "BiomS", 
             "simulation" : "Ref",
             "reference" : "Ref"}

labs = ["{}, {}".format(simu_labs[simu], space_labs[sp]) for simu in simuls for sp in spaces]
idx =pd.MultiIndex.from_product([simuls, spaces], names=["space", "simu"])

nn_smpls = np.arange(343)
nsats = len(names)

if __name__ == "__main__":
    """
    Comparing refMet4 simulations in spaces gtdb207 and msp
    """
    plt.close("all")
    indir = "../analyses/combined_metrics/"

    figdir = "../figures/"
    if not os.path.exists(figdir):
        os.makedirs(figdir)

    #--------------------------------------------------------------------------
    #Fig 1. Lost abundance and clades            
    ratio=1.; top = 0.93; bottom = 0.16; hspace=0.2 
    Ly, lx, ly = plot_height(figwidth, nsats +1, 2, ratio, bottom=bottom, top=top, hspace=hspace)

    lost_summary = pd.read_csv(indir + "lost_summary.tsv", sep="\t", index_col=[0,1,2,3]).swaplevel()
    lost_clades_frac = lost_summary.xs("Lost_mOTUs_frac", level=1)
    lost_clades = lost_summary.xs("Lost_mOTUs", level=1)
    lost_abundance = lost_summary.xs("Lost_abundance", level=1)

    fig, axs = plt.subplots(2, nsats+1, figsize=(figwidth, np.round(Ly, 2)),\
     sharey="row", sharex=True)
    for j, (i, lab) in enumerate(zip(idx, labs)):
        lost_ab = lost_abundance.xs(i, level=[1,2])
        lost_cl = lost_clades.xs(i, level=[1,2])
        for ax, fname in zip(axs.T, ["reference"] + fnames):
            ax[1].plot(nn_smpls, lost_ab[fname].sort_values(ascending=False), color=cols1[j], label=lab)
            ax[0].plot(nn_smpls, lost_cl[fname].sort_values(ascending=False), color=cols1[j], label=lab)

    axs[1,0].set_ylabel("Lost relative abundance")
    axs[0,0].set_ylabel("Lost features")
    [ax.set_yscale("symlog", linthresh=10.**-5) for ax in axs[1]]
    [ax.set_yscale("symlog", linthresh=10.**-3) for ax in axs[0]]
    [ax.set_title(tool_labs[name]) for ax, name in zip(axs[0], ["reference"] + fnames)]
    handles, labels = axs[0,0].get_legend_handles_labels() 
    fig.legend(handles, labels,\
     ncol=nsats, bbox_to_anchor=(0.5, 0.), loc="lower center")

    fig.text(0.005, .99, "a", va="top", ha="left")
    fig.text(0.005, .57, "b", va="top", ha="left")
    fig.text(0.5, 0.09, "Samples", ha="center")
    plt.subplots_adjust(top=top, bottom=bottom, hspace=hspace)
    plt.savefig(figdir + "Fig1a_lost_abundance_clades_log.png")
    plt.close(fig)

    #--------------------------------------------------------------------------
    #Fig 2. Richness and Shannon diversity

    richness = pd.read_csv(indir + "richness.tsv", sep="\t", index_col=(0,1,2)).swaplevel()
    shannon = pd.read_csv(indir + "shannon.tsv", sep="\t", index_col=(0,1,2)).swaplevel()

    richness_diff = richness.sub(richness["simulation"], axis=0)
    shannon_diff = shannon.sub(shannon["simulation"], axis=0)


    ratio=1.; top = 0.93; bottom = 0.13; left=0.08; hspace=0.2
    Ly, lx, ly = plot_height(figwidth, nsats, 2, ratio, bottom=bottom, top=top, left=left, hspace=hspace)
    Ly = np.round(Ly, 2)

    fig, axs = plt.subplots(2, nsats, figsize=(figwidth, Ly), sharey="row", sharex=True)
    for name, ax in zip(names, axs[0]):
        colored_box_swarm_plot(ax, richness_diff[name], idx, cols1, ms=1.)
        ax.set_title(tool_labs[name])
        ax.set_xticks(np.arange(1, len(labs)+1))
        ax.set_xticklabels([])
        ax.set_yticks(np.arange(-600, 601, 300))
        ax.axhline(0., ls="--", c="k", lw=0.5)

    for name, ax in zip(names, axs[1]):
        colored_box_swarm_plot(ax, shannon_diff[name], idx, cols1, ms=1.)
        ax.set_xticks(np.arange(1, len(labs)+1))
        ax.set_xticklabels([])
        ax.set_yticks(np.arange(-2, 3, 2))
        ax.axhline(0., ls="--", c="k", lw=0.5)

    handles = [mpatches.Patch(color=c, label=lab) for c, lab in zip(cols1, labs)]
    fig.legend(handles, labs,\
     ncol=nsats, bbox_to_anchor=(0.5, 0.), loc="lower center")

    axs[0,0].set_ylabel("Richness\n(estimate - reference)")
    axs[1,0].set_ylabel("Shannon diversity\n(estimate - reference)")

    fig.text(0.5, 0.08, "Simulation + projection space", ha="center", va="bottom")
    fig.text(0.01, .99, "a", va="top", ha="right")
    fig.text(0.005, .57, "b", va="top")
    # plt.suptitle("Species richness")
    plt.subplots_adjust(bottom=bottom, top=top, left=left, hspace=hspace)
    plt.savefig(figdir + "Fig2_richness_shannon_diff.png")
    plt.close(fig)

    #--------------------------------------------------------------------------
    #FigS7 Richness and Shannon diversity on log scale
    richness_log = np.log10(richness)
    shannon_log = np.log10(shannon)

    ratio=1.; top = 0.93; bottom = 0.13; left=0.08; hspace=0.2
    Ly, lx, ly = plot_height(figwidth, nsats+1, 2, ratio, bottom=bottom, top=top, left=left, hspace=hspace)
    Ly = np.round(Ly, 2)

    fig, axs = plt.subplots(2, nsats+1, figsize=(figwidth, Ly), sharey="row", sharex=True)
    for name, ax in zip(['simulation'] + names, axs[0]):
        colored_box_swarm_plot(ax, richness_log[name], idx, cols1, ms=1.)
        ax.set(
            title=tool_labs[name],
            xticks=np.arange(1, len(labs)+1),
            xticklabels=[],
            yticks = [1., 2., 3.],
            yticklabels=[r'$10^1$', r'$10^2$', r'$10^3$']
            )

    for name, ax in zip(['simulation'] + names, axs[1]):
        colored_box_swarm_plot(ax, shannon_log[name], idx, cols1, ms=1.)
        ax.set(
            xticks=np.arange(1, len(labs)+1),
            xticklabels=[],
            yticks=[0., 1.],
            yticklabels=[r'$10^0$', r'$10^1$']
            )

    handles = [mpatches.Patch(color=c, label=lab) for c, lab in zip(cols1, labs)]
    fig.legend(handles, labs,\
     ncol=nsats, bbox_to_anchor=(0.5, 0.), loc="lower center")

    axs[0,0].set_ylabel("Richness")
    axs[1,0].set_ylabel("Shannon diversity")

    fig.text(0.5, 0.08, "Simulation + projection space", ha="center", va="bottom")
    fig.text(0.01, .99, "a", va="top", ha="right")
    fig.text(0.005, .57, "b", va="top")
    # plt.suptitle("Species richness")
    plt.subplots_adjust(bottom=bottom, top=top, left=left, hspace=hspace)
    plt.savefig(figdir + "FigS7_richness_shannon_diff_log.png")
    plt.close(fig)

    #--------------------------------------------------------------------------
    #Fig 3. Bray-Curtis and UniFrac distance
    dist_to_ref = pd.read_csv(indir + "dist_to_ref.tsv", sep="\t", index_col=(0,1,2)).swaplevel()
    dist_to_ref_wu = pd.read_csv(indir + "dist_to_ref_wu.tsv", sep="\t", index_col=(0,1,2)).swaplevel()

    ratio=1.; top = 0.95; bottom = 0.05; left=0.06
    Ly, lx, ly = plot_height(figwidth, nsats, nsats+1, ratio, bottom=bottom, top=top, left=left)
    Ly = np.round(Ly, 2)

    fig, axs = plt.subplots(idx.size, nsats, figsize=(figwidth, Ly), sharex=True, sharey=True)
    for ax, i, lab, c in zip(axs, idx, labs, cols):
        for name, a in zip(names, ax):
            a.scatter(dist_to_ref.xs(i, level=(1,2))[name],\
                      dist_to_ref_wu.xs(i, level=(1,2))[name],\
                      label=lab, color=c, s=3)
    fig.text(0.5, 0.01, "Bray-Curtis distance", ha="center")
    fig.text(0.01, 0.5, "UniFrac distance", va="center", rotation="vertical")
    for ax, lab in zip(axs[:,-1], labs):
        ax.yaxis.set_label_position("right")
        ax.set_ylabel(lab)#, fontsize=7)
    [ax.set_title(tool_labs[name]) for ax, name in zip(axs[0], names)]
    plt.subplots_adjust(bottom=bottom, top=top, left=left)
    plt.savefig(figdir + "FigS1_distance_scatter.png")
    plt.close(fig)

    #--------------------------------------------------------------------------
    #Fig 5. Precision-Recall plane

    superconf = pd.read_csv(indir + "superconf.tsv", sep="\t", index_col=[0,1,2,3]).swaplevel()
    TPR = superconf.xs("TPR", level=1)
    PPV = superconf.xs("PPV", level=1)

    # for plotting ellipse see https://matplotlib.org/stable/gallery/statistics/confidence_ellipse.html
    from matplotlib.patches import Ellipse
    import matplotlib.transforms as transforms
    def confidence_ellipse(x, y, ax, n_std=3.0, facecolor='none', **kwargs):
        """
        Create a plot of the covariance confidence ellipse of *x* and *y*.

        Parameters
        ----------
        x, y : array-like, shape (n, )
            Input data.

        ax : matplotlib.axes.Axes
            The axes object to draw the ellipse into.

        n_std : float
            The number of standard deviations to determine the ellipse's radiuses.

        **kwargs
            Forwarded to `~matplotlib.patches.Ellipse`

        Returns
        -------
        matplotlib.patches.Ellipse
        """
        if x.size != y.size:
            raise ValueError("x and y must be the same size")

        cov = np.cov(x, y)
        pearson = cov[0, 1]/np.sqrt(cov[0, 0] * cov[1, 1])
        # Using a special case to obtain the eigenvalues of this
        # two-dimensional dataset.
        ell_radius_x = np.sqrt(1 + pearson)
        ell_radius_y = np.sqrt(1 - pearson)
        ellipse = Ellipse((0, 0), width=ell_radius_x * 2, height=ell_radius_y * 2,
                          facecolor=facecolor, **kwargs)

        # Calculating the standard deviation of x from
        # the squareroot of the variance and multiplying
        # with the given number of standard deviations.
        scale_x = np.sqrt(cov[0, 0]) * n_std
        mean_x = np.mean(x)

        # calculating the standard deviation of y ...
        scale_y = np.sqrt(cov[1, 1]) * n_std
        mean_y = np.mean(y)

        transf = transforms.Affine2D() \
            .rotate_deg(45) \
            .scale(scale_x, scale_y) \
            .translate(mean_x, mean_y)

        ellipse.set_transform(transf + ax.transData)
        return ax.add_patch(ellipse)

    ratio=1.; top = 0.9; bottom = 0.25
    Ly, lx, ly = plot_height(figwidth, nsats, 1, ratio, bottom=bottom, top=top)
    Ly = np.round(Ly, 2)

    fig, axs = plt.subplots(1, nsats, figsize=(figwidth, Ly), sharex=True, sharey=True)
    for name, ax in zip(names, axs):
        for c, i, lab, m in zip(cols1, idx, labs, marks):
            xx, yy = TPR.xs(i, level=(1,2))[name], PPV.xs(i, level=(1,2))[name]
            # starplot(ax, xx, yy, color=c, edgecolor=c, label=lab, marker=m, ms_big=10, small_markers=True)
            ax.scatter(xx, yy, color=c, marker=m, s=2.)
            cface = c.copy()
            cface[-1] = 0.2
            confidence_ellipse(xx, yy, ax, n_std=3., edgecolor=c, facecolor=cface, lw=0.5, label=lab)
        ax.set_title(tool_labs[name])
    axs[0].set_ylabel("Precision")
    fig.text(0.5, 0.14, "Sensitivity", ha="center", va="bottom")

    handles, labels = ax.get_legend_handles_labels() 
    fig.legend(handles, labels,\
     ncol=nsats, bbox_to_anchor=(0.5, 0.), loc="lower center")
    plt.subplots_adjust(bottom=bottom, top=top)
    plt.savefig(figdir + "Fig4_precision_recall.png")
    plt.close(fig)

    #--------------------------------------------------------------------------
    #Fig3 : Quartile plot of Bray-Curtis and UniFrac distances

    richness_quarts = pd.read_csv(indir + "quarts/richness_diff.tsv".format(name), sep="\t", index_col=(0,1,2))
    shannon_quarts = pd.read_csv(indir + "quarts/shannon_diff.tsv".format(name), sep="\t", index_col=(0,1,2))
    bc_quarts = pd.read_csv(indir + "quarts/braycurtis.tsv".format(name), sep="\t", index_col=(0,1,2))
    unifrac_quarts = pd.read_csv(indir + "quarts/unifrac.tsv".format(name), sep="\t", index_col=(0,1,2))
    nnx = np.arange(1, len(labs)+1)

    ratio=1.; top = 0.9; bottom = 0.34; left=0.1; wspace=0.2
    Ly, lx, ly = plot_height(figwidth/2., 2, 1, ratio, bottom=bottom, top=top, left=left, wspace=wspace)
    Ly = np.round(Ly, 2)

    fig, axs = plt.subplots(1, 2, figsize=(figwidth/2., Ly), sharex=True, sharey=True)
    for name in names:
        axs[0].plot(nnx, bc_quarts.xs(name, level=0).loc[idx, "q50"], label=tool_labs[name])
        axs[0].fill_between(nnx, bc_quarts.xs(name, level=0).loc[idx, "q25"],\
                            bc_quarts.xs(name, level=0).loc[idx, "q75"],\
                            alpha=0.5)

        axs[1].plot(nnx, unifrac_quarts.xs(name, level=0).loc[idx, "q50"], label=tool_labs[name])
        axs[1].fill_between(nnx, unifrac_quarts.xs(name, level=0).loc[idx, "q25"],\
                            unifrac_quarts.xs(name, level=0).loc[idx, "q75"],\
                            alpha=0.5)
    for ax in axs:
        ax.set_xticks(nnx)
        ax.set_xticklabels(labs, rotation=45, ha="right")
    axs[0].set_title("Bray-Curtis")
    axs[1].set_title("UniFrac")
    axs[0].set_ylabel("Distance (estimation to reference)")

    handles, labels = ax.get_legend_handles_labels() 
    fig.legend(handles, labels,\
     ncol=nsats, bbox_to_anchor=(0.5, 0.), loc="lower center")

    fig.text(0.5, 0.12, "Simulation + projection space", ha="center")
    fig.text(0.06, .99, "a", va="top")
    fig.text(0.53, .99, "b", va="top")

    plt.subplots_adjust(bottom=bottom, top=top, left=left, wspace=wspace)
    plt.savefig(figdir + "Fig3_quarts_bc_unifrac.png")
    plt.close(fig)

    #--------------------------------------------------------------------------
    #Fig5 : FPRA and FNRA plot
    FPRA = superconf.xs("FPRA", level=1)
    FNRA = superconf.xs("FNRA", level=1)

    ratio=1.; top = 0.93; bottom = 0.13; left=0.08; hspace=0.2
    Ly, lx, ly = plot_height(figwidth, nsats, 2, ratio, bottom=bottom, top=top, left=left, hspace=hspace)
    Ly = np.round(Ly, 2)

    fig, axs = plt.subplots(2, nsats, figsize=(figwidth, Ly), sharey="row", sharex=True)
    for name, ax in zip(names, axs[0]):
        colored_box_swarm_plot(ax, FPRA[name], idx, cols1, ms=1.)
        ax.set_title(tool_labs[name])
        # ax.set_xticks(np.arange(1, len(labs)+1))
        ax.set_xticklabels([])

    for name, ax in zip(names, axs[1]):
        colored_box_swarm_plot(ax, FNRA[name], idx, cols1, ms=1.)
        # ax.set_title(tool_labs[name])
        # ax.set_xticks(np.arange(1, len(labs)+1))
        ax.set_xticklabels([])

    handles = [mpatches.Patch(color=c, label=lab) for c, lab in zip(cols1, labs)]
    fig.legend(handles, labs,\
     ncol=nsats, bbox_to_anchor=(0.5, 0.), loc="lower center")
    axs[0,0].set_ylabel("False positive relative abundance")
    axs[1,0].set_ylabel("False negative relative abundance")
    fig.text(0.5, 0.08, "Simulation + projection space", ha="center", va="bottom")
    fig.text(0.005, .99, "a", va="top", ha="left")
    fig.text(0.005, .55, "b", va="top", ha="left")

    plt.subplots_adjust(bottom=bottom, top=top, left=left, hspace=hspace)
    plt.savefig(figdir + "Fig5_FPRA_boxplot.png")
    plt.close(fig)


    #--------------------------------------------------------------------------
    #plots for confused species

    confused_gtdb = ["d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Bacteroidales;f__Bacteroidaceae;g__Prevotella;s__Prevotella copri",
                    "d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Bacteroidales;f__Bacteroidaceae;g__Prevotella;s__Prevotella sp015074785"]
    confused_names = [spec.split(";s__")[1] for spec in confused_gtdb]

    ratio=1.; top = 0.94; bottom = 0.14; hspace=0.2
    Ly, lx, ly = plot_height(figwidth, nsats, 2, ratio, bottom=bottom, top=top, hspace=hspace)
    Ly = np.round(Ly, 2)

    cmap_loc = plt.rcParams['axes.prop_cycle'].by_key()['color']
    xx = np.linspace(0., 0.8)

    #Prevotella plot
    spec1, spec2 = confused_gtdb
    namespec1, namespec2 = confused_names
    fig, axs = plt.subplots(2, nsats, figsize=(figwidth, Ly), sharey=True, sharex=True)
    for simu, axs_row in zip(['refMet4', 'refKrak'], axs):
        supersat = pd.read_csv("../data/{0}_to_gtdb207/supersat_{0}_gtdb207.sat".format(simu),\
        sep="\t", index_col=[0,1])

        for jname, (name, fname, ax) in enumerate(zip(names, fnames, axs_row)):
            ax.set_title(tool_labs[name])
            ax.plot(xx, xx, ls="--", c="lightgrey")
            x_ref = supersat.xs(spec2, level=0)["reference"]
            y1_est = supersat.xs(spec1, level=0)[fname]
            y2_est = supersat.xs(spec2, level=0)[fname]

            X = pd.concat((x_ref, y1_est, y2_est), axis=1).to_numpy()
            for x, y1, y2 in X:
                ax.plot([x, x], [y1, y2], color="lightgrey", lw=0.25)

            ax.plot(x_ref, y1_est, ls="", marker="o", color=cmap_loc[0], label=namespec1, markersize=2)
            ax.plot(x_ref, y2_est, ls="", marker="s", color=cmap_loc[1], label=namespec2, markersize=2)

            ax.set_yscale("symlog", linthresh=.0001)
            ax.set_xscale("symlog", linthresh=.0001)
        
        axs_row[-1].yaxis.set_label_position("right")
        axs_row[-1].set_ylabel('{}, GTDB'.format(simu_labs[simu]))

    handles, labels = ax.get_legend_handles_labels() 
    fig.legend(handles, labels, bbox_to_anchor=(0.5, 0.),\
     ncol=4, loc="lower center", markerscale=3)

    fig.text(0.5, 0.07, "Relative abundance of {} (reference)".format(namespec2), ha="center", va="bottom")
    fig.text(0.01, 0.5, "Relative abundance (estimated)", va="center", ha="left", rotation="vertical")
    fig.text(0.04, .98, "a", va="top", ha="right")
    fig.text(0.04, .55, "b", va="top", ha="right")

    plt.subplots_adjust(bottom=bottom, top=top, hspace=hspace)
    plt.savefig(figdir + "Fig6A_prevotella.png".format(simu))
    plt.close(fig)


    #Frequently confused species
    dg_wide = pd.read_csv('../analyses/confusion_pairs/frequently_confused_to_plot.tsv',\
     sep='\t', index_col=[0,1], header=[0,1,2])
    FPs = dg_wide.index.get_level_values(0)
    FNs = dg_wide.index.get_level_values(1)
    FP_labels = [FP.split(';s__')[1] for FP in FPs]
    FN_labels = [FN.split(';s__')[1] for FN in FNs]
    # pair_labels = ['{} vs. {}'.format(FP.split(';s__')[1], FN.split(';s__')[1]) for FP, FN in zip(FPs, FNs)]

    ratio=2.; top = 0.89; bottom = 0.25; left = 0.17; right=0.82; wspace=0.1
    Ly, lx, ly = plot_height(figwidth, nsats, 1, ratio, bottom=bottom, top=top, left=left, right=right, wspace=wspace)
    Ly = np.round(Ly, 2)

    nn = np.arange(dg_wide.shape[0])
    nn_sats = np.arange(nsats)
    fig, axs = plt.subplots(1, nsats, figsize=(figwidth, Ly), sharex=True)
    for i, lab, ax in zip(idx, labs, axs):
        dg_loc = dg_wide.xs(i, level=[0,1], axis=1).reindex(columns=fnames)
        im=ax.imshow(dg_loc, cmap='viridis', vmin=.9, vmax=1., aspect='auto')
        ax.grid(False)
        ax.set(
            xticks=nn_sats,
            xticklabels=[tool_labs[name] for name in fnames],
            xlabel=lab
        )
        ax.tick_params(axis='x', labelrotation=90)
        ax.xaxis.set_label_position("top")
    
    axs[0].set(
        yticks=nn,
        yticklabels=FN_labels
    )
    axs[-1].yaxis.set_label_position("right")
    axs[-1].yaxis.tick_right()
    axs[-1].set(
        yticks=nn,
        yticklabels=FP_labels
    )
    [ax.set_yticks([]) for ax in axs[1:-1]]
    cax = fig.add_axes([0.32, 0.05, 0.42, 0.04])
    plt.colorbar(im, cax=cax, orientation='horizontal', label='Pearson r')
    fig.text(0.32, 0.07, 'Pearson R:  ', ha='right', va='center')
    
    fig.text(0.165, 0.9, 'Present in simulation', ha='right', va='bottom')
    fig.text(0.825, 0.9, 'Present in estimate', ha='left', va='bottom')
    fig.text(0.5, 0.99, "Frequent confusion pairs", ha="center", va="top")
    fig.text(0.04, .98, "c", va="top", ha="right")
    plt.subplots_adjust(bottom=bottom, top=top, left=left, right=right, wspace=wspace)
    plt.savefig(figdir + "Fig6B_prevotella.png")
    plt.close(fig)

    #--------------------------------------------------------------------------
    #Beta-diversity
    def pairdist_components_plot0(V1, V2, ax, linewidth=0.5, linecolor="b", markersize=2):
        for (_, v1), (_, v2) in zip(V1.iterrows(), V2.iterrows()):
            ax.plot([v1["PC1"], v2["PC1"]], [v1["PC2"], v2["PC2"]],\
             color=linecolor, lw=linewidth)

        p1 = ax.plot(V1["PC1"], V1["PC2"], ls="", marker="o", markersize=markersize, label="a")
        p2 = ax.plot(V2["PC1"], V2["PC2"], ls="", marker="s", markersize=markersize, label="b")

        return p1, p2

    fnames_loc = ["motus3", "metaphlan4"]
    names_loc = ["motus3", "metaphlan4"]

    ratio=1.; top = 0.93; bottom = 0.18
    Ly, lx, ly = plot_height(figwidth, nsats+1, len(names_loc), ratio, bottom=bottom, top=top)
    Ly = np.round(Ly, 2)

    fig, axs = plt.subplots(len(names_loc), idx.size, figsize=(figwidth, Ly), sharex=True, sharey=True)
    for ax, i, lab, c in zip(axs.T, idx, labs, cols):
        indir_loc = "../analyses/article_results_{}_{}/braycurtis/".format(*i)
        V_ref = pd.read_csv(indir_loc + "V_ref.tsv", sep="\t", index_col=0)
        VV_est = [pd.read_csv(indir_loc + "V_{}_matched.tsv".format(fname), sep="\t", index_col=0) for fname in fnames_loc]

        if i == ("refMet4", "gtdb207"):
            V_ref["PC2"] = - V_ref["PC2"]
            VV_est[0]["PC2"] = - VV_est[0]["PC2"]
            VV_est[1]["PC2"] = - VV_est[1]["PC2"]
        elif i == ("refKrak", "uhgg"):
            # V_ref["PC2"] = - V_ref["PC2"]
            VV_est[0]["PC2"] = - VV_est[0]["PC2"]
        elif i == ("refMet4", "uhgg"):
            V_ref["PC2"] = - V_ref["PC2"]
            VV_est[0]["PC2"] = - VV_est[0]["PC2"]
            VV_est[1]["PC2"] = - VV_est[1]["PC2"]

        for name, V_est, a in zip(names_loc, VV_est, ax):
            p1, p2 = pairdist_components_plot0(V_ref, V_est, a,\
             linecolor="grey", markersize=0.5, linewidth=0.3)
    handles, labels = a.get_legend_handles_labels() 
    fig.legend(handles, ["Reference", "Estimate"], markerscale=3,\
     ncol=2, bbox_to_anchor=(0.5, 0.), loc="lower center")

    fig.text(0.5, 0.09, "Principal coordinate 1", ha="center", va="bottom")
    fig.text(0.01, 0.5, "Principal coordinate 2", va="center", ha="left", rotation="vertical")
    for ax, name in zip(axs[:,-1], names_loc):
        ax.yaxis.set_label_position("right")
        ax.set_ylabel(tool_labs[name])

    [ax.set_title(lab) for ax, lab in zip(axs[0], labs)]
    plt.subplots_adjust(bottom=bottom, top=top)
    plt.savefig(figdir + "FigS2_beta_diversity.png")
    plt.close(fig)

    #--------------------------------------------------------------------------
    #distances between distance matrices
    norms_bc = []
    for i in idx:
        indir_loc = "../analyses/article_results_{}_{}/braycurtis/".format(*i)
        df = pd.read_csv(indir_loc + "norm_ref_matched.tsv", sep="\t", index_col=0)
        df[["simu", "space"]] = list(i)
        norms_bc.append(df)

    norms_bc = pd.concat(norms_bc).set_index(["simu", "space"], append=True)
    norms_bc.to_csv(indir + "norms_bc.tsv", sep="\t")

    norms_uf = []
    for i in idx:
        indir_loc = "../analyses/article_results_{}_{}/unifrac/".format(*i)
        df = pd.read_csv(indir_loc + "norm_ref_matched_wu.tsv", sep="\t", index_col=0)
        df[["simu", "space"]] = list(i)
        norms_uf.append(df)

    norms_uf = pd.concat(norms_uf).set_index(["simu", "space"], append=True)
    norms_uf.to_csv(indir + "norms_uf.tsv", sep="\t")


    norms_ait = []
    for i in idx:
        indir_loc = "../analyses/article_results_{}_{}/aitchison/".format(*i)
        df = pd.read_csv(indir_loc + "norm_ref_matched.tsv", sep="\t", index_col=0)
        df[["simu", "space"]] = list(i)
        norms_ait.append(df)

    norms_ait = pd.concat(norms_ait).set_index(["simu", "space"], append=True)
    norms_ait.to_csv(indir + "norms_ait.tsv", sep="\t")


    # # ratio=.2; top = 0.85; bottom = 0.15
    # ratio=.2; top = 0.94; bottom = 0.16; hspace=0.2
    # Ly, lx, ly = plot_height(figwidth, 1, 2, ratio, bottom=bottom, top=top, hspace=hspace)
    # Ly = np.round(Ly, 2)

    # cmap_default = plt.rcParams['axes.prop_cycle'].by_key()['color']
    # cmap_default = np.array([mpl.colors.to_rgb(x) for x in cmap_default])
    # cmap_default = np.c_[cmap_default, np.ones(cmap_default.shape[0])]
    # cmap_loc = cmap_default[[2, 3, 4, 1, 0]]
    # cmap_loc[:3, 3] = 0.5
    # nn = np.arange(idx.size * (nsats+1))
    # norm = "JSD"
    # norm_name = "Jensen-Shannon distance"
    # fig, axs = plt.subplots(2, 1, figsize=(figwidth, Ly), sharex=True)
    # col = norms_bc[norm]
    # for j, name in enumerate(names):
    #     axs[0].bar(nn[j::(nsats+1)], col.xs(name, level=0).values,\
    #                   align="edge", color=cmap_loc[j])

    # ymax = col.max()
    # for j, i in enumerate(idx):
    #     n = names.index(col.xs(i, level=(1,2)).idxmin())
    #     xmax = j*(nsats+1) + n + 0.5
    #     axs[0].plot(xmax, ymax, marker="*", markersize=4, color=cmap_loc[n])

    # col = norms_uf[norm]
    # for j, name in enumerate(names):
    #     axs[1].bar(nn[j::(nsats+1)], col.xs(name, level=0).values, label=tool_labs[name],\
    #             align="edge", color=cmap_loc[j])

    # ymax = col.max()
    # for j, i in enumerate(idx):
    #     n = names.index(col.xs(i, level=(1,2)).idxmin())
    #     xmax = j*(nsats+1) + n + 0.5
    #     axs[1].plot(xmax, ymax, marker="*", markersize=4, color=cmap_loc[n])

    # # axs[0].set_ylabel(norm_name)
    # axs[1].set_xticks(nn[::(nsats+1)])
    # axs[1].set_xticklabels(labs, ha="left")
    # fig.legend(bbox_to_anchor=(0.5, 0.), ncol=nsats, loc="lower center")
    # fig.text(0.5, 0.08, "Simulation + projection space", ha="center", va="bottom")
    # fig.text(0.01, 0.5, norm_name, va="center", ha="left", rotation="vertical")
    # # fig.text(0.01, .99, "a", va="top")
    # # fig.text(0.5, 0.08, "Relative abundance (reference)", ha="center", va="bottom")
    # # fig.text(0.01, 0.5, "Relative abundance (estimated)", va="center", ha="left", rotation="vertical")
    # fig.text(0.04, .99, "a", va="top", ha="right")
    # fig.text(0.04, .56, "b", va="top", ha="right")

    # plt.subplots_adjust(bottom=bottom, top=top, hspace=hspace)
    # plt.savefig(figdir + "Fig7_JSD.png")
    # plt.close(fig)

    #including aitchison
    ratio=.2; top = 0.94; bottom = 0.11; hspace=0.2
    Ly, lx, ly = plot_height(figwidth, 1, 3, ratio, bottom=bottom, top=top, hspace=hspace)
    Ly = np.round(Ly, 2)

    cmap_default = plt.rcParams['axes.prop_cycle'].by_key()['color']
    cmap_default = np.array([mpl.colors.to_rgb(x) for x in cmap_default])
    cmap_default = np.c_[cmap_default, np.ones(cmap_default.shape[0])]
    cmap_loc = cmap_default[[2, 3, 4, 1, 0]]
    cmap_loc[:3, 3] = 0.5
    nn = np.arange(idx.size * (nsats+1))
    norm = "JSD"
    norm_name = "Jensen-Shannon distance"
    fig, axs = plt.subplots(3, 1, figsize=(figwidth, Ly), sharex=True)
    col = norms_bc[norm]
    for j, name in enumerate(names):
        axs[0].bar(nn[j::(nsats+1)], col.xs(name, level=0).values,\
                      align="edge", color=cmap_loc[j])

    ymax = col.max()
    for j, i in enumerate(idx):
        n = names.index(col.xs(i, level=(1,2)).idxmin())
        xmax = j*(nsats+1) + n + 0.5
        axs[0].plot(xmax, ymax, marker="*", markersize=4, color=cmap_loc[n])

    col = norms_uf[norm]
    for j, name in enumerate(names):
        axs[1].bar(nn[j::(nsats+1)], col.xs(name, level=0).values,\
                align="edge", color=cmap_loc[j])

    ymax = col.max()
    for j, i in enumerate(idx):
        n = names.index(col.xs(i, level=(1,2)).idxmin())
        xmax = j*(nsats+1) + n + 0.5
        axs[1].plot(xmax, ymax, marker="*", markersize=4, color=cmap_loc[n])

    col = norms_ait[norm]
    for j, name in enumerate(names):
        axs[2].bar(nn[j::(nsats+1)], col.xs(name, level=0).values, label=tool_labs[name],\
                align="edge", color=cmap_loc[j])

    ymax = col.max()
    for j, i in enumerate(idx):
        n = names.index(col.xs(i, level=(1,2)).idxmin())
        xmax = j*(nsats+1) + n + 0.5
        axs[2].plot(xmax, ymax, marker="*", markersize=4, color=cmap_loc[n])

    axs[-1].set_xticks(nn[::(nsats+1)])
    axs[-1].set_xticklabels(labs, ha="left")
    fig.legend(bbox_to_anchor=(0.5, 0.), ncol=nsats, loc="lower center")
    fig.text(0.5, 0.06, "Simulation + projection space", ha="center", va="bottom")
    fig.text(0.01, 0.5, norm_name, va="center", ha="left", rotation="vertical")
    fig.text(0.04, .97, "a", va="top", ha="right")
    fig.text(0.04, .67, "b", va="top", ha="right")
    fig.text(0.04, .37, "c", va="top", ha="right")

    plt.subplots_adjust(bottom=bottom, top=top, hspace=hspace)
    plt.savefig(figdir + "Fig7_JSD_aitchison.png")
    plt.close(fig)


    #--------------------------------------------------------------------------
    #Raw richness 

    rich = pd.read_csv('../analyses/raw_richness/raw_richness.tsv', sep='\t', index_col=[0,1]).squeeze()
    shan = pd.read_csv('../analyses/raw_richness/raw_shannon.tsv', sep='\t', index_col=[0,1]).squeeze()

    ratio=1.; top = 0.91; bottom = 0.15; left=0.1; wspace=0.3
    Ly, lx, ly = plot_height(figwidth/2., 2, 1, ratio, bottom=bottom, top=top, left=left, wspace=wspace)
    Ly = np.round(Ly, 2)

    fig, axs = plt.subplots(1, 2, figsize=(figwidth/2., Ly), sharex=True)
    colored_box_swarm_plot(axs[0], rich,\
     pd.Index(['refKrak', 'refMet4']), cols1[1::2], levs=[0], ms=1)
    colored_box_swarm_plot(axs[1], shan,\
     pd.Index(['refKrak', 'refMet4']), cols1[1::2], levs=[0], ms=1)

    axs[0].set(
        ylabel='Richness',
        xticks=[1, 2],
        xticklabels=('refKrak', 'refMet4')
    )

    axs[1].set(
        ylabel='Shannon diversity',
        xticks=[1, 2],
        xticklabels=('refKrak', 'refMet4')
    )

    fig.text(0.5, 0.99, 'Raw species richness and Shannon diversity in simulations', ha='center', va='top')
    fig.text(0.5, 0.02, "Simulation", ha="center", va="bottom")

    plt.subplots_adjust(bottom=bottom, top=top, left=left, wspace=wspace)
    plt.savefig(figdir + "FigS5_raw_richness_shannon.png")
    plt.close(fig)


    #--------------------------------------------------------------------------
    #Phyla abundance and prevalence in simulation
    phylum_abundance = pd.read_csv('../analyses/raw_richness/phylum_abundance.tsv', sep='\t', index_col=[0,1])
    phylum_presence = pd.read_csv('../analyses/raw_richness/phylum_presence.tsv', sep='\t', index_col=[0,1])

    phyla = phylum_abundance.unstack(level=0).mean(axis=1).sort_values(ascending=False).index.tolist()
    ratio=1.; top = 0.94; bottom = 0.13; left=0.08; hspace=0.02; wspace=0.02
    Ly, lx, ly = plot_height(figwidth, 2, 2, ratio, bottom=bottom, top=top, left=left, hspace=hspace, wspace=wspace)
    Ly = np.round(Ly, 2)

    nphyla, nsmpls = phylum_presence.shape
    nn_smpls = np.arange(nsmpls)
    cmap_loc = mpl.cm.get_cmap("tab20")
    cols_loc = [cmap_loc([j])[0] for j in range(20)]

    fig, axs = plt.subplots(2, 2, figsize=(figwidth, Ly), sharey="row", sharex=True)
    for simu, axs_col in zip(['refKrak', 'refMet4'], axs.T):
        df = phylum_abundance.xs(simu, level=0)
        df.sort_values(by=phyla, axis=1, ascending=False, inplace=True)
        base = np.zeros(nsmpls)
        for p, c in zip(phyla, cols_loc):
            x = df.loc[p].to_numpy()
            axs_col[1].bar(nn_smpls, x, width=1., linewidth=0., bottom=base, align='edge', color=c)
            base += x 

        axs_col[1].set(
            xlim=(0., 343),
            ylim=(0., 1.),
            xticks=[]
        )

        df = phylum_presence.xs(simu, level=0)
        df.sort_values(by=phyla, axis=1, ascending=False, inplace=True)
        base = np.zeros(nsmpls)
        for p, c in zip(phyla, cols_loc):
            x = df.loc[p].to_numpy()
            axs_col[0].bar(nn_smpls, x, width=1., linewidth=0., bottom=base, align='edge', color=c)
            base += x

        axs_col[0].set(
            xlim=(0., 343),
            xticks=[],
            title =simu_labs[simu] 
        )

    axs[0,0].set_ylabel('#Species')
    axs[1,0].set_ylabel('Abundance')
    handles = [mpatches.Patch(color=c, label=p) for c, p in zip(cols_loc, phyla)]
    fig.legend(handles, phyla,\
     ncol=4, bbox_to_anchor=(0.5, 0.01), loc="lower center")

    # axs[0,0].set_ylabel("Richness\n(estimate - reference)")
    # axs[1,0].set_ylabel("Shannon diversity\n(estimate - reference)")

    fig.text(0.5, 0.99, "Phyla presence and abundance simulations", ha="center", va="top")
    fig.text(0.5, 0.1, "Samples", ha="center", va="bottom")

    # fig.text(0.01, .99, "a", va="top", ha="right")
    # fig.text(0.005, .57, "b", va="top")
    plt.subplots_adjust(bottom=bottom, top=top, left=left, hspace=hspace, wspace=wspace)
    plt.savefig(figdir + "FigS6_phyla_abundance_presence.png")
    plt.close(fig)
