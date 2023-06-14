#!/usr/bin/python3.6

from __future__ import division
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
# import seaborn as sns
from scipy.spatial.distance import squareform, pdist#, braycurtis, cosine, jensenshannon
# from scipy.stats import spearmanr, pearsonr
from skbio.stats.distance import mantel
from skbio.diversity.beta import unweighted_unifrac, weighted_unifrac
from skbio import read
from skbio.tree import TreeNode
from Bio import Phylo
import sys, os

from metrics_aux import confusion_matrix
# from compare_kraken import pcoa_dist

# sns.set_theme()

marks = ["o", "s", "d", "^", "*", "+", "x", "."]

#setting figure parameters
fsize = 2
mpl.rcParams["font.size"]=10

from scipy.linalg import eigh
def pcoa_dist(dat_dist, to_file='', k=2, clusters={}):
    '''
    Perform MDS/PCoA decomposition
    Arguments:
    df : data frame
    to_file : path to file for saving figure
    k : number of MDS components to calculate
    '''
    n = len(dat_dist)
    labs = dat_dist.columns.tolist()

    H = np.eye(n) - np.ones((n,n))/n
    ZZ = -0.5*H.dot(dat_dist**2).dot(H)
    W, V = eigh(ZZ, subset_by_index = (n-k, n-1))
    V = V[:,::-1]
    W = W[::-1]
    if to_file:
        fig, ax = plt.subplots(1, 1, figsize = (10, 10))
        if clusters:
            ax.scatter(V[:,0], V[:,1], marker='.')
            for (name, jj), c in zip(clusters.items(), cols):
                ax.scatter(V[jj,0], V[jj,1], color=c, label=name)
            ax.legend()
        else:
            ax.scatter(V[:,0], V[:,1])
        ax.set_xlabel('component 1')
        ax.set_ylabel('component 2')
        for j, smpl in enumerate(labs):
            ax.annotate(smpl, (V[j,0], V[j,1]), fontsize=6)

        plt.savefig(to_file)
        plt.close(fig)
    return W, V

def pairwise_nmds(sats, smpls, metric):
    DD = []; VV = []
    for sat in sats:
        if metric == "cosine_sqrt":
            D_bc = pd.DataFrame(squareform(pdist(np.sqrt(sat[smpls]).T, "cosine")), index=smpls, columns=smpls)
        else:
            D_bc = pd.DataFrame(squareform(pdist(sat[smpls].T, metric)), index=smpls, columns=smpls)
        _, V_bc = pcoa_dist(D_bc, k=3)

        DD.append(D_bc)
        VV.append(V_bc)
    return DD, VV

def pairdist_components_plot0(V1, V2, ax):
    p1 = ax.scatter(V1[:,0], V1[:,1], marker="s")
    p2 = ax.scatter(V2[:,0], V2[:,1], marker="d")

    for v1, v2 in zip(V1, V2):
        ax.plot([v1[0], v2[0]], [v1[1], v2[1]], color="b", lw=.5)

    return p1, p2

def pairdist_scatterplot0(D1, D2, ax):
    ii = np.triu_indices(D2.shape[0], k=1)
    x = D1.to_numpy()[ii]
    y = D2.to_numpy()[ii]
    # r = pearsonr(x, y)[0]
    # rs = spearmanr(x, y)[0]
    # ax.scatter(x, y, label="Pearson r = {:.4f},\nSpearman r = {:.4f}".format(r, rs))
    r = mantel(D1, D2, method="pearson")[0]
    p = ax.scatter(x, y, label="r = {:.4f}".format(r))

    # xx = np.linspace(x.min(), x.max())
    xx = np.linspace(0.2, 1.)
    ax.plot(xx, xx, ls="--", color="k")

def pairdist_blandaltman0(D1, D2, ax):
    ii = np.triu_indices(D2.shape[0], k=1)
    x = D1.to_numpy()[ii]
    y = D2.to_numpy()[ii]

    mean = (x + y)/2.
    diff = x-y
    ax.scatter(mean, diff, marker="o")
    md = diff.mean()
    sd = diff.std()
    ax.axhline(md, ls="--", color="k")
    ax.axhline(md + 1.96*sd, ls="--", color="k")
    ax.axhline(md - 1.96*sd, ls="--", color="k")

def build_unifrac_matrices(sat, ids, tree, norm=True, validate=True):
    # ids = df.index.tolist()
    # X = sat.loc[df["gtdb_taxonomy"]].to_numpy()
    smpls = sat.columns
    nsmpls = len(smpls)
    X = sat.to_numpy()
    X = np.round(X / X[X>0].min())

    D_uu = np.zeros((nsmpls, nsmpls))
    D_wu = np.zeros((nsmpls, nsmpls))
    for ju, u in enumerate(X.T):
        for jv, v in enumerate(X.T[:ju]):
            D_uu[ju, jv] = unweighted_unifrac(u, v, ids, tree, validate=validate)
            D_wu[ju, jv] = weighted_unifrac(u, v, ids, tree, normalized=norm, validate=validate)
    # D_ref_uu += D_ref_uu.T
    # D_ref_wu += D_ref_wu.T
    D_uu = pd.DataFrame(D_uu + D_uu.T, index=smpls, columns=smpls)
    D_wu = pd.DataFrame(D_wu + D_wu.T, index=smpls, columns=smpls)
    return D_uu, D_wu

def unifrac_dist_to_ref(sat0, sat, ids, tree, norm=True, validate=True):
    smpls = sat.columns
    nsmpls = len(smpls)
    X0 = sat0.to_numpy()
    X = sat.to_numpy()

    fctr = np.min([X0[X0>0].min(), X[X>0].min()])
    X0 = np.round(X0 / fctr)
    X = np.round(X / fctr)

    df = []
    for u, v in zip(X0.T, X.T):
        uu = unweighted_unifrac(u, v, ids, tree, validate=validate)
        wu = weighted_unifrac(u, v, ids, tree, normalized=norm, validate=validate)
        df.append((uu, wu))

    df = pd.DataFrame(df, index=smpls, columns=["unweighted_unifrac", "weighted_unifrac"])
    return df

ctxt1 = {
"figure.figsize": (11, 8),
"axes.titlesize": 12,
"axes.labelsize": 8,
"xtick.labelsize": 6,
"ytick.labelsize": 6,
"legend.fontsize": 6,
# "font.size": 8,
"lines.markersize": 2,
"lines.linewidth": 1,
"figure.subplot.hspace": .3,
"figure.subplot.wspace": 0,
"figure.subplot.right": 0.99,
"figure.subplot.left": 0.05,
"figure.subplot.top": 0.95,
"figure.subplot.bottom": 0.06}
