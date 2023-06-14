#!/usr/bin/python3.6

from __future__ import division
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import os, sys

#--------------------------------------------------------------------------
#figures settings
plt.style.use('ggplot')
plt.style.use('seaborn-whitegrid')

cmap = mpl.cm.get_cmap("Paired")
marks = ["o", "s", "^", "d", "*", "+", "x", "."]
cols = [cmap([j])[0] for j in range(10)]

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
    ny : number fo panels in y direction
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

#--------------------------------------------------------------------------
#names and labels
names = ["kraken2 + bracken", "metaphlan3", "motus3", "metaphlan4", "biomscope"]
fnames = ["kraken2", "metaphlan3", "motus3", "metaphlan4", "biomscope"]
simuls = ["refKrak", "refBioms", "refMet4"]
spaces = ["gtdb207", "igc2"]

space_labs = {"gtdb207": "GTDB", "igc2": "IGC2"}
tool_labs = {"kraken2 + bracken": "Krak/Brac", 
             "kraken2" : "Krak/Brac",
             "metaphlan3": "MPA3", 
             "metaphlan4": "MPA4", 
             "motus3": "mOTUs3", 
             "biomscope": "BiomS", 
             "simulation" : "Ref",
             "reference" : "Ref"}

labs = ["{}, {}".format(simu, space_labs[sp]) for simu in simuls for sp in spaces]
idx =pd.MultiIndex.from_product([simuls, spaces], names=["space", "simu"])

nn_smpls = np.arange(343)
nsats = len(names)

if __name__ == "__main__":
    """
    Comparing simuCRC simulations in spaces gtdb207 and msp
    """
    plt.close("all")
    indir = "../analyses/combined_metrics/"

    figdir = "../figures/final_figures/"
    if not os.path.exists(figdir):
        os.makedirs(figdir)

    #--------------------------------------------------------------------------
    #Fig 2. Lost abundance and clades            

    lost_abundance = pd.read_csv(indir + "lost_abundance.tsv", sep="\t", index_col=[0,1,2]).swaplevel()
    lost_clades = pd.read_csv(indir + "lost_clades.tsv", sep="\t", index_col=[0,1,2]).swaplevel()

    ratio=1.; top = 0.8; bottom = 0.1 
    Ly, lx, ly = plot_height(figwidth, nsats +1, 2, ratio, bottom=bottom, top=top)

    fig, axs = plt.subplots(2, nsats+1, figsize=(figwidth, np.round(Ly, 2)),\
     sharey="row", sharex=True)
    for j, (i, lab) in enumerate(zip(idx, labs)):
        lost_ab = lost_abundance.xs(i, level=[1,2])
        lost_cl = lost_clades.xs(i, level=[1,2])
        # for ax, (name, col) in zip(axs[0], lost_ab.items()):
        for ax, fname in zip(axs.T, ["reference"] + fnames):
            ax[0].plot(nn_smpls, lost_ab[fname].sort_values(ascending=False), color=cmap([j])[0], label=lab)
            ax[1].plot(nn_smpls, lost_cl[fname].sort_values(ascending=False), color=cmap([j])[0], label=lab)

        # for ax, (name, col) in zip(axs[1], lost_cl.items()):
            # ax.plot(nn_smpls, col.sort_values(ascending=False), color=cmap([j])[0], label=lab)

    axs[0,0].set_ylabel("Lost relative abundance")
    axs[1,0].set_ylabel("Lost features")
    [ax.set_yscale("symlog", linthreshy=10.**-5) for ax in axs[0]]
    [ax.set_yscale("symlog", linthreshy=1) for ax in axs[1]]
    [ax.set_title(tool_labs[name]) for ax, name in zip(axs[0], lost_ab.columns)]
    # [ax.legend() for ax in axs[:,-1]]
    handles, labels = axs[0,0].get_legend_handles_labels() 
    fig.legend(handles, labels,\
     ncol=3, bbox_to_anchor=(0.5, 1.), loc="upper center")

    fig.text(0.5, 0.01, "Samples", ha="center")
    # plt.xlabel("Samples", ha="center")
    plt.subplots_adjust(top=top, bottom=bottom)
    plt.savefig(figdir + "Fig2_lost_abundance_clades_log.png")
    plt.close(fig)

    #--------------------------------------------------------------------------
    #Fig 3. Richness and Shannon diversity

    def colored_boxplot(ax, data, idx, cols):
        for j, (i, c) in enumerate(zip(idx, cols)):
            # c = cmap([j])[0]
            ax.boxplot(data.xs(i, level=(1,2)), widths=.5, positions=[j+1.], notch=True,\
                       boxprops=dict(color=c),
                       capprops=dict(color=c),
                       whiskerprops=dict(color=c),
                       flierprops=dict(color=c, markeredgecolor=c),
                       medianprops=dict(color=c))

    def colored_box_swarm_plot(ax, data, idx, cols, ms=.5):
        for j, (i, c) in enumerate(zip(idx, cols)):
            yy = data.xs(i, level=(1,2)).to_numpy()
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


    richness = pd.read_csv(indir + "richness.tsv", sep="\t", index_col=(0,1,2)).swaplevel()
    shannon = pd.read_csv(indir + "shannon.tsv", sep="\t", index_col=(0,1,2)).swaplevel()

    richness_diff = richness.sub(richness["simulation"], axis=0)
    shannon_diff = shannon.sub(shannon["simulation"], axis=0)

    ratio=2.; top = 0.9; bottom = 0.22; left=0.1
    Ly, lx, ly = plot_height(figwidth/2., nsats, 2, ratio, bottom=bottom, top=top, left=left)
    Ly = np.round(Ly, 2)

    fig, axs = plt.subplots(1, nsats, figsize=(figwidth/2., Ly), sharey="row", sharex=True)
    for name, ax in zip(names, axs):
        colored_box_swarm_plot(ax, richness_diff[name], idx, cols)
        ax.set_title(tool_labs[name])
        ax.set_xticks(np.arange(1, len(labs)+1))
        ax.set_xticklabels(labs, rotation=90)

    axs[0].set_ylabel("Richness (estimate-reference)")
    fig.text(0.5, 0.01, "Simulation + projection space", ha="center")
    fig.text(0.01, .99, "a", va="top")
    # plt.suptitle("Species richness")
    plt.subplots_adjust(bottom=bottom, top=top, left=left)
    plt.savefig(figdir + "Fig3A_richness_diff.png")
    plt.close(fig)

    fig, axs = plt.subplots(1, nsats, figsize=(figwidth/2., Ly), sharey="row", sharex=True)
    for name, ax in zip(names, axs):
        colored_box_swarm_plot(ax, shannon_diff[name], idx, cols)
        ax.set_title(tool_labs[name])
        ax.set_xticks(np.arange(1, len(labs)+1))
        ax.set_xticklabels(labs, rotation=90)

    axs[0].set_ylabel("Shannon diversity (estimate-reference)")
    fig.text(0.5, 0.01, "Simulation + projection space", ha="center")
    fig.text(0.01, .99, "b", va="top")
    # plt.suptitle("Shannon diversity")
    plt.subplots_adjust(bottom=bottom, top=top, left=left)
    plt.savefig(figdir + "Fig3B_shannon_diff.png")
    plt.close(fig)

    #--------------------------------------------------------------------------
    #Fig 4. Bray-Curtis and UniFrac distance
    dist_to_ref = pd.read_csv(indir + "dist_to_ref.tsv", sep="\t", index_col=(0,1,2)).swaplevel()
    dist_to_ref_wu = pd.read_csv(indir + "dist_to_ref_wu.tsv", sep="\t", index_col=(0,1,2)).swaplevel()

    ratio=1.; top = 0.92; bottom = 0.1; left=0.1
    Ly, lx, ly = plot_height(figwidth/2., nsats, nsats+1, ratio, bottom=bottom, top=top, left=left)
    Ly = np.round(Ly, 2)

    fig, axs = plt.subplots(idx.size, nsats, figsize=(figwidth/2., Ly), sharex=True, sharey=True)
    for ax, i, lab, c in zip(axs, idx, labs, cols):
        for name, a in zip(names, ax):
            a.scatter(dist_to_ref.xs(i, level=(1,2))[name],\
                      dist_to_ref_wu.xs(i, level=(1,2))[name],\
                      label=lab, color=c)
    fig.text(0.5, 0.01, "Bray-Curtis distance", ha="center")
    fig.text(0.01, 0.5, "UniFrac distance", va="center", rotation="vertical")
    fig.text(0.01, .99, "b", va="top")
    for ax, lab in zip(axs[:,-1], labs):
        ax.yaxis.set_label_position("right")
        ax.set_ylabel(lab, fontsize=7)
    [ax.set_title(tool_labs[name]) for ax, name in zip(axs[0], names)]
    plt.subplots_adjust(bottom=bottom, top=top, left=left)
    # plt.suptitle("UniFrac vs. Bray-Curtis distance")
    plt.savefig(figdir + "Fig4B_distance_scatter.png")
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

    ratio=1.; top = 0.73; bottom = 0.14
    Ly, lx, ly = plot_height(figwidth, nsats, 1, ratio, bottom=bottom, top=top)
    Ly = np.round(Ly, 2)

    # with mpl.rc_context(ctxt2):
    fig, axs = plt.subplots(1, nsats, figsize=(figwidth, Ly), sharex=True, sharey=True)
    for name, ax in zip(names, axs):
        for c, i, lab, m in zip(cols, idx, labs, marks):
            xx, yy = TPR.xs(i, level=(1,2))[name], PPV.xs(i, level=(1,2))[name]
            # starplot(ax, xx, yy, color=c, edgecolor=c, label=lab, marker=m, ms_big=10, small_markers=True)
            ax.scatter(xx, yy, color=c, marker=m, s=2.)
            cface = c.copy()
            cface[-1] = 0.2
            confidence_ellipse(xx, yy, ax, n_std=3., edgecolor=c, facecolor=cface, lw=0.5, label=lab)
        ax.set_title(tool_labs[name])
    axs[0].set_ylabel("Precision")
    fig.text(0.5, 0.01, "Sensitivity", ha="center", va="bottom")

    handles, labels = ax.get_legend_handles_labels() 
    fig.legend(handles, labels,\
     ncol=3, bbox_to_anchor=(0.5, 1.), loc="upper center")
    # axs[-1].legend()
    plt.subplots_adjust(bottom=bottom, top=top)
    plt.savefig(figdir + "Fig5_precision_recall.png")
    plt.close(fig)


    #--------------------------------------------------------------------------
    #Quartile plots

    richness_quarts = pd.read_csv(indir + "quarts/richness_diff.tsv".format(name), sep="\t", index_col=(0,1,2))
    shannon_quarts = pd.read_csv(indir + "quarts/shannon_diff.tsv".format(name), sep="\t", index_col=(0,1,2))
    bc_quarts = pd.read_csv(indir + "quarts/braycurtis.tsv".format(name), sep="\t", index_col=(0,1,2))
    unifrac_quarts = pd.read_csv(indir + "quarts/unifrac.tsv".format(name), sep="\t", index_col=(0,1,2))
    nnx = np.arange(1, len(labs)+1)

    ratio=2.; top = 0.88; bottom = 0.18; left=0.1
    Ly, lx, ly = plot_height(figwidth/2., 2, 1, ratio, bottom=bottom, top=top, left=left)
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
        ax.set_xticklabels(labs, rotation=90)
    axs[0].set_title("Bray-Curtis")
    axs[1].set_title("UniFrac")
    axs[0].set_ylabel("Distance (estimation to reference)")

    handles, labels = ax.get_legend_handles_labels() 
    fig.legend(handles, labels,\
     ncol=3, bbox_to_anchor=(0.5, 1.), loc="upper center")
    # axs[0].legend(fontsize=6)

    fig.text(0.5, 0.01, "Simulation + projection space", ha="center")
    fig.text(0.01, .99, "a", va="top")

    plt.subplots_adjust(bottom=bottom, top=top, left=left)
    plt.savefig(figdir + "Fig4A_quarts_bc_unifrac.png")
    plt.close(fig)

    #--------------------------------------------------------------------------
    #FPRA plot
    FPRA = superconf.xs("FPRA", level=1)
    FPRA_quarts = pd.read_csv(indir + "quarts/FPRA.tsv".format(name), sep="\t", index_col=(0,1,2))

    ratio=2.; top = 0.92; bottom = 0.22
    Ly, lx, ly = plot_height(figwidth, nsats, 1, ratio, bottom=bottom, top=top)
    Ly = np.round(Ly, 2)

    fig, axs = plt.subplots(1, nsats, figsize=(figwidth, Ly), sharey="row", sharex=True)
    for name, ax in zip(names, axs):
        # colored_boxplot(ax, FPRA[name], idx, cols)
        colored_box_swarm_plot(ax, FPRA[name], idx, cols)
        ax.set_title(tool_labs[name])
        ax.set_xticks(np.arange(1, len(labs)+1))
        ax.set_xticklabels(labs, rotation=90)

    # axs[0].set_ylabel("FPRA")
    axs[0].set_ylabel("False positive relative abundance")
    fig.text(0.5, 0.01, "Simulation + projection space", ha="center")
    plt.subplots_adjust(bottom=bottom, top=top)
    plt.savefig(figdir + "Fig6_FPRA_boxplot.png")
    plt.close(fig)

    #--------------------------------------------------------------------------
    #plots for confused species
    confused_gtdb = ["d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Enterobacterales;f__Enterobacteriaceae;g__Escherichia;s__Escherichia coli",
                "d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Enterobacterales;f__Enterobacteriaceae;g__Escherichia;s__Escherichia fergusonii",
                "d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Bacteroidales;f__Bacteroidaceae;g__Prevotella;s__Prevotella copri",
                "d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Bacteroidales;f__Bacteroidaceae;g__Prevotella;s__Prevotella copri_A",
                "d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Bacteroidales;f__Bacteroidaceae;g__Prevotella;s__Prevotella sp015074785"]
    confused_names = [spec.split(";s__")[1] for spec in confused_gtdb]


    ratio=1.; top = 0.79; bottom = 0.16
    Ly, lx, ly = plot_height(figwidth, nsats, 1, ratio, bottom=bottom, top=top)
    Ly = np.round(Ly, 2)

    alpha=0.5
    cmap_loc = plt.rcParams['axes.prop_cycle'].by_key()['color']
    # cmap_loc = cols[::2]
    xx = np.linspace(0., 0.8)
    for jsimu, simu in enumerate(simuls):
        # indir_loc = "../data/{}_to_gtdb/".format(simu)
        supersat = pd.read_csv("../data/{0}_to_gtdb207/supersat_{0}_gtdb207.sat".format(simu),\
         sep="\t", index_col=[0,1])
        fpra = superconf.xs(("FPRA", simu, "gtdb207"), level=(1,2,3))

        #Prevotella plot
        spec1, spec2, spec3 = confused_gtdb[2:] 
        namespec1, namespec2, namespec3 = confused_names[2:]
        leg = [namespec1 + " (reference)", namespec2 + " (reference)", namespec3 + " (reference)",\
         namespec1 + " (estimate)", namespec2 + " (estimate)", namespec3 + " (estimate)" ]
        fig, axs = plt.subplots(1, nsats, figsize=(figwidth, Ly), sharey=True, sharex="row")
        for jname, (name, fname, ax) in enumerate(zip(names, fnames, axs)):
            ax.set_title(tool_labs[name])
            ax.plot(xx, xx, ls="--", c="lightgrey")
            x = fpra[name]
            # y1_ref = supersat.xs(spec1, level=0)["reference"]
            # y2_ref = supersat.xs(spec2, level=0)["reference"]
            y3_ref = supersat.xs(spec3, level=0)["reference"]
            y1_est = supersat.xs(spec1, level=0)[fname]
            y2_est = supersat.xs(spec2, level=0)[fname]
            y3_est = supersat.xs(spec3, level=0)[fname]

            #draw segments
            # y_max = pd.concat((y1_ref, y2_ref, y3_ref, y1_est, y2_est, y3_est), axis=1).max(axis=1)
            # Y = np.array([np.zeros(y_max.shape), y_max.to_numpy()])
            yy = pd.concat((y1_est, y2_est, y3_est, x), axis=1).to_numpy()
            Y = np.sort(yy, axis=1)[:, -2:].T
            X = np.array([y3_ref.to_numpy(), y3_ref.to_numpy()])
            ii = np.nonzero(X[0,:]>0)[0]
            ax.plot(X[:,ii], Y[:,ii], color="lightgrey", lw=0.25)

            ax.plot(y3_ref, y1_est, ls="", marker="o", color=cmap_loc[0], label=namespec1, alpha=alpha)
            ax.plot(y3_ref, y2_est, ls="", marker="s", color=cmap_loc[1], label=namespec2, alpha=alpha)
            ax.plot(y3_ref, y3_est, ls="", marker="^", color=cmap_loc[2], label=namespec3, alpha=alpha)
            ax.plot(y3_ref, x, ls="", marker="d", color=cmap_loc[3], label="FPRA", alpha=alpha)
            # ax.set_xlabel(namespec3)

            ax.set_yscale("symlog", linthreshy=.0001)
            ax.set_xscale("symlog", linthreshx=.0001)

        handles, labels = ax.get_legend_handles_labels() 
        fig.legend(handles, labels, bbox_to_anchor=(0.5, 1.),\
         ncol=4, loc="upper center", markerscale=3)

        # [ax.set_ylabel("Relative species abundance") for ax in axs[:,0]]
        fig.text(0.5, 0.01, "Relative abundance of {} (reference)".format(namespec3), ha="center", va="bottom")
        fig.text(0.01, 0.5, "Relative abundance (estimated)", va="center", ha="left", rotation="vertical")
        fig.text(0.01, .99, "a", va="top")
        plt.subplots_adjust(bottom=bottom, top=top)
        plt.savefig(figdir + "Fig7A_prevotella_copri_1_{}.png".format(simu))
        plt.close(fig)


        fig, axs = plt.subplots(1, nsats, figsize=(figwidth, Ly), sharey=True, sharex="row")
        for jname, (name, fname, ax) in enumerate(zip(names, fnames, axs)):
            ax.set_title(tool_labs[name])
            ax.plot(xx, xx, ls="--", c="lightgrey")
            x = fpra[name]
            # y1_ref = supersat.xs(spec1, level=0)["reference"]
            y2_ref = supersat.xs(spec2, level=0)["reference"]
            # y3_ref = supersat.xs(spec3, level=0)["reference"]
            y1_est = supersat.xs(spec1, level=0)[fname]
            y2_est = supersat.xs(spec2, level=0)[fname]
            y3_est = supersat.xs(spec3, level=0)[fname]

            #draw segments
            # y_max = pd.concat((y1_ref, y2_ref, y3_ref, y1_est, y2_est, y3_est), axis=1).max(axis=1)
            # Y = np.array([np.zeros(y_max.shape), y_max.to_numpy()])
            yy = pd.concat((y1_est, y2_est, y3_est, x), axis=1).to_numpy()
            Y = np.sort(yy, axis=1)[:, -2:].T
            X = np.array([y2_ref.to_numpy(), y2_ref.to_numpy()])
            ii = np.nonzero(X[0,:]>0)[0]
            ax.plot(X[:,ii], Y[:,ii], color="lightgrey", lw=0.25)

            ax.plot(y2_ref, y1_est, ls="", marker="o", color=cmap_loc[0], label=namespec1, alpha=alpha)
            ax.plot(y2_ref, y2_est, ls="", marker="s", color=cmap_loc[1], label=namespec2, alpha=alpha)
            ax.plot(y2_ref, y3_est, ls="", marker="^", color=cmap_loc[2], label=namespec3, alpha=alpha)
            ax.plot(y2_ref, x, ls="", marker="d", color=cmap_loc[3], label="FPRA", alpha=alpha)
            # ax.set_xlabel(namespec3)

            ax.set_yscale("symlog", linthreshy=.0001)
            ax.set_xscale("symlog", linthreshx=.0001)

        handles, labels = ax.get_legend_handles_labels() 
        fig.legend(handles, labels, bbox_to_anchor=(0.5, 1.),\
         ncol=4, loc="upper center", markerscale=3)

        # [ax.set_ylabel("Relative species abundance") for ax in axs[:,0]]
        fig.text(0.5, 0.01, "Relative abundance of {} (reference)".format(namespec2), ha="center", va="bottom")
        fig.text(0.01, 0.5, "Relative abundance (estimated)", va="center", ha="left", rotation="vertical")
        fig.text(0.01, .99, "b", va="top")
        plt.subplots_adjust(bottom=bottom, top=top)
        plt.savefig(figdir + "Fig7B_prevotella_copri_2_{}.png".format(simu))
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

    fnames_loc = ["metaphlan4", "biomscope"]
    names_loc = ["metaphlan4", "biomscope"]

    ratio=1.; top = 0.85; bottom = 0.1
    Ly, lx, ly = plot_height(figwidth, nsats+1, len(names_loc), ratio, bottom=bottom, top=top)
    Ly = np.round(Ly, 2)

    fig, axs = plt.subplots(len(names_loc), idx.size, figsize=(figwidth, Ly), sharex=True, sharey=True)
    for ax, i, lab, c in zip(axs.T, idx, labs, cols):
        indir_loc = "../analyses/article_results_{}_{}/braycurtis/".format(*i)
        V_ref = pd.read_csv(indir_loc + "V_ref.tsv", sep="\t", index_col=0)
        VV_est = [pd.read_csv(indir_loc + "V_{}_matched.tsv".format(fname), sep="\t", index_col=0) for fname in fnames_loc]

        for name, V_est, a in zip(names_loc, VV_est, ax):
            p1, p2 = pairdist_components_plot0(V_ref, V_est, a,\
             linecolor="grey", markersize=0.5, linewidth=0.3)
    handles, labels = a.get_legend_handles_labels() 
    fig.legend(handles, ["Reference", "Estimate"], markerscale=3,\
     ncol=2, bbox_to_anchor=(0.5, 1.), loc="upper center")

    fig.text(0.5, 0.01, "Principal coordinate 1", ha="center", va="bottom")
    fig.text(0.01, 0.5, "Principal coordinate 2", va="center", ha="left", rotation="vertical")
    for ax, name in zip(axs[:,-1], names_loc):
        ax.yaxis.set_label_position("right")
        ax.set_ylabel(tool_labs[name])

    [ax.set_title(lab) for ax, lab in zip(axs[0], labs)]
    plt.subplots_adjust(bottom=bottom, top=top)
    plt.savefig(figdir + "FigS1_beta_diversity.png")
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


    ratio=.2; top = 0.85; bottom = 0.15
    Ly, lx, ly = plot_height(figwidth, 1, 1, ratio, bottom=bottom, top=top)
    Ly = np.round(Ly, 2)

    cmap_default = plt.rcParams['axes.prop_cycle'].by_key()['color']
    cmap_default = np.array([mpl.colors.to_rgb(x) for x in cmap_default])
    cmap_default = np.c_[cmap_default, np.ones(cmap_default.shape[0])]
    cmap_loc = cmap_default[[2, 3, 4, 1, 0]]
    cmap_loc[:3, 3] = 0.5
    nn = np.arange(idx.size * (nsats+1))
    norm = "JSD"
    col = norms_bc[norm]
    fig, ax = plt.subplots(1, 1, figsize=(figwidth, Ly))
    for j, name in enumerate(names):
        ax.bar(nn[j::(nsats+1)], col.xs(name, level=0).values, label=tool_labs[name],\
                      align="edge", color=cmap_loc[j])

    ymax = col.max()
    for j, i in enumerate(idx):
        n = names.index(col.xs(i, level=(1,2)).idxmin())
        xmax = j*(nsats+1) + n + 0.5
        ax.plot(xmax, ymax, marker="*", markersize=4, color=cmap_loc[n])

    ax.set_ylabel(norm)
    ax.set_xticks(nn[::(nsats+1)])
    ax.set_xticklabels(labs, ha="left")
    fig.legend(bbox_to_anchor=(0.5, 1.), ncol=nsats, loc="upper center")
    fig.text(0.5, 0.01, "Simulation + projection space", ha="center", va="bottom")
    fig.text(0.01, .99, "a (Bray-Curtis)", va="top")
    plt.subplots_adjust(bottom=bottom, top=top)
    # plt.suptitle(norm)
    plt.savefig(figdir + "Fig8A_JSD_braycurtis.png")
    plt.close(fig)

    col = norms_uf[norm]
    fig, ax = plt.subplots(1, 1, figsize=(figwidth, Ly))
    for j, name in enumerate(names):
        ax.bar(nn[j::(nsats+1)], col.xs(name, level=0).values, label=tool_labs[name],\
                align="edge", color=cmap_loc[j])

    ymax = col.max()
    for j, i in enumerate(idx):
        n = names.index(col.xs(i, level=(1,2)).idxmin())
        xmax = j*(nsats+1) + n + 0.5
        ax.plot(xmax, ymax, marker="*", markersize=4, color=cmap_loc[n])

    ax.set_ylabel(norm)
    ax.set_xticks(nn[::(nsats+1)])
    ax.set_xticklabels(labs, ha="left")
    fig.legend(bbox_to_anchor=(0.5, 1.), ncol=nsats, loc="upper center")
    fig.text(0.5, 0.01, "Simulation + projection space", ha="center", va="bottom")
    fig.text(0.01, .99, "b (UniFrac)", va="top")
    plt.subplots_adjust(bottom=bottom, top=top)
    # plt.suptitle(norm)
    plt.savefig(figdir + "Fig8B_JSD_unifrac.png")
    plt.close(fig)
