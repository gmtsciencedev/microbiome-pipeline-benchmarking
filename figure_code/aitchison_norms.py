#!/usr/bin/python3.6

from __future__ import division
import numpy as np
import pandas as pd
from skbio.stats.composition import clr
from scipy.spatial.distance import euclidean
import sys, os

from metrics_aux import *

def replace_zeros(SAT, d=10**(-7), method='simple'):
    if method == 'simple':
        X = SAT.copy()
        X[X == 0] += d
        X = X.divide(X.sum(), axis=1)
    elif method == 'additive':
        Z = (SAT == 0).sum()
        D = SAT.shape[0]
        X1 = SAT[SAT > 0].sub(d * (Z + 1) * Z / D**2, axis=1).filna(0.)
        X2 = SAT[SAT > 0].add(d * (Z + 1) * (D - Z) / D**2).fillna(0.)
        X = X1 + X2
    elif method == 'multiplicative':
        X=SAT[SAT > 0].mul(1. - d*(SAT == 0).sum(), axis=1).fillna(d)
    return X


simu = "refMet4"#"refKrak"#"refMet4"
space = "uhgg"#"uhgg"#"gtdb207"
names = ["kraken2 + bracken", "metaphlan3", "motus3", "metaphlan4"]

indir = "../data/{}_to_{}/".format(simu, space)
fnames = ["kraken2", "metaphlan3", "motus3", "metaphlan4"]

supersat = pd.read_csv(indir + "supersat_{}_{}.sat".format(simu, space), sep="\t", index_col=(0,1))
matched_sats = [supersat[name].unstack() for name in fnames]
sat_ref = supersat["reference"].unstack()
smpls = sat_ref.columns
nsmpls = len(smpls)

nsats = len(matched_sats)
nrows = 1
ncols = nsats

if __name__ == "__main__":
    """
    Article figures
    """
    metric = "aitchison"
    outdir = "../analyses/article_results_{0}_{1}/{2}/".format(simu, space, metric)
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    #--------------------------------------------------------------------
    #non-metric multidimensional scaling (NMDS)
    d = 10.**(-7)#10.**(-7)
    zero_replacement_method = "simple"#"multiplicative"
    nPC = 3 #number of principal components

    X_ref = replace_zeros(sat_ref[smpls], d=d, method=zero_replacement_method)
    logX_ref = clr(X_ref)
    D_ref = pd.DataFrame(squareform(pdist(logX_ref.T, "euclidean")), index=smpls, columns=smpls)
    _, V_ref = pcoa_dist(D_ref, k=nPC)
    V_ref = pd.DataFrame(V_ref, index=smpls, columns=["PC{}".format(j+1) for j in range(nPC)])
    D_ref.to_csv(outdir + "D_ref.tsv", sep="\t")
    V_ref.to_csv(outdir + "V_ref.tsv", sep="\t")


    # DD_matched, VV_matched = pairwise_nmds(matched_sats, smpls, metric)
    DD_matched = []
    VV_matched = []
    XX_matched = []
    for sat in matched_sats:
        X_matched = replace_zeros(sat[smpls], d=d, method=zero_replacement_method)
        XX_matched.append(X_matched.copy())
        logX_matched = clr(X_matched)
        D_matched = pd.DataFrame(squareform(pdist(logX_matched.T, "euclidean")), index=smpls, columns=smpls)
        DD_matched.append(D_matched)
        _, V_matched = pcoa_dist(D_matched, k=nPC)
        V_matched = pd.DataFrame(V_matched, index=smpls, columns=["PC{}".format(j+1) for j in range(nPC)])
        VV_matched.append(V_matched)

    #saving distance matrices and PCoA components
    for fname, D_matched, V_matched in zip(fnames, DD_matched, VV_matched):
        D_matched.to_csv(outdir + "D_{}_matched.tsv".format(fname), sep="\t")
        V_matched.to_csv(outdir + "V_{}_matched.tsv".format(fname), sep="\t")


    #--------------------------------------------------------------------
    #calculating intermatrix norms
    norm_ref_matched = intermatrix_norms([D_ref]*len(names), DD_matched, names, filepath=outdir +"norm_ref_matched.tsv")

    #--------------------------------------------------------------------
    #distance between real sample and simulation
    dist_to_ref = pd.DataFrame(pd.DataFrame({name : [euclidean(x_est, x_ref) for x_est, x_ref in zip(clr(X_matched).T, logX_ref.T)] for name, X_matched in zip(names, XX_matched)}, index=smpls))
    dist_to_ref.to_csv(outdir + "dist_to_ref.tsv", sep="\t")

