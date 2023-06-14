#!/usr/bin/python3.6

from __future__ import division
import numpy as np
import pandas as pd
from scipy.spatial.distance import braycurtis, cosine, pdist, squareform, jensenshannon
# from scipy.stats import spearmanr, pearsonr
# from skbio.stats.distance import mantel
import sys, os

from metrics_aux import *
# from article_figures import *


simu = "refKrak"#"simuCRC2k"#"simuCRC2b"#"simuCRC"
space = "gtdb207"#"msp"#"gtdb207"
names = ["kraken2 + bracken", "metaphlan3", "motus3", "metaphlan4", "biomscope"]

indir = "../data/{}_to_{}/".format(simu, space)
fnames = ["kraken2", "metaphlan3", "motus3", "metaphlan4", "biomscope"]

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
    metric = "braycurtis"#"cosine_sqrt"#"jensenshannon"#"braycurtis"
    outdir = "../analyses/article_results_{0}_{1}/{2}/".format(simu, space, metric)
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    #--------------------------------------------------------------------
    #non-metric multidimensional scaling (NMDS)
    nPC = 3 #number of principal components
    if metric == "cosine_sqrt":
        D_ref = pd.DataFrame(squareform(pdist(np.sqrt(sat_ref[smpls]).T, "cosine")), index=smpls, columns=smpls)
    else:
        D_ref = pd.DataFrame(squareform(pdist(sat_ref[smpls].T, metric)), index=smpls, columns=smpls)
    _, V_ref = pcoa_dist(D_ref, k=nPC)
    V_ref = pd.DataFrame(V_ref, index=smpls, columns=["PC{}".format(j+1) for j in range(nPC)])

    DD_matched, VV_matched = pairwise_nmds(matched_sats, smpls, metric)

    #saving distance matrices and PCoA components
    D_ref.to_csv(outdir + "D_ref.tsv", sep="\t")
    V_ref.to_csv(outdir + "V_ref.tsv", sep="\t")

    for fname, D_matched, V_matched in zip(fnames, DD_matched, VV_matched):
        D_matched.to_csv(outdir + "D_{}_matched.tsv".format(fname), sep="\t")
        V_matched.to_csv(outdir + "V_{}_matched.tsv".format(fname), sep="\t")


    #--------------------------------------------------------------------
    #calculating intermatrix norms
    norm_ref_matched = intermatrix_norms([D_ref]*len(names), DD_matched, names, filepath=outdir +"norm_ref_matched.tsv")

    #--------------------------------------------------------------------
    #distance between real sample and simulation
    if metric == "braycurtis":
        dist_to_ref = pd.DataFrame(pd.DataFrame({name : [braycurtis(sat[smpl], sat_ref[smpl]) for smpl in smpls] for name, sat in zip(names, matched_sats)}, index=smpls))
        # lab="Bray-Curtis distance"
    elif metric == "jensenshannon":
        dist_to_ref = pd.DataFrame(pd.DataFrame({name : [jensenshannon(sat[smpl], sat_ref[smpl]) for smpl in smpls] for name, sat in zip(names, matched_sats)}, index=smpls))
        # lab="Jensen-Shannon distance"
    elif metric == "cosine_sqrt":
        dist_to_ref = pd.DataFrame(pd.DataFrame({name : [cosine_sqrt(sat[smpl], sat_ref[smpl]) for smpl in smpls] for name, sat in zip(names, matched_sats)}, index=smpls))
        # lab="cosine-sqrt distance"

    dist_to_ref.to_csv(outdir + "dist_to_ref.tsv", sep="\t")
    # dist_to_ref.sort_values(by="biomscope", inplace=True)

    #--------------------------------------------------------------------
    #Richness
    richnesses = [(sat[smpls] > 0).sum(axis=0) for sat in matched_sats]
    richness_ref = (sat_ref[smpls] > 0).sum(axis=0)

    df_richness = pd.concat([richness_ref] + richnesses, axis=1)
    df_richness.columns = ["simulation"] + names
    df_richness.to_csv(outdir + "../richness.tsv")


    #Shannon diversity
    h = 10**(-8)
    shannon_ref = - np.sum(sat_ref[smpls] * np.log2(sat_ref[smpls] + h), axis=0)
    shannon_norm_ref = shannon_ref / np.log2(richness_ref)
    shannons = [- np.sum(sat[smpls] * np.log2(sat[smpls] + h), axis=0) for sat in matched_sats]
    shannons_norm = [s / np.log2(r) for s, r in zip(shannons, richnesses)]

    df_shannon = pd.concat([shannon_ref] + shannons, axis=1)
    df_shannon.columns = ["simulation"] + names
    df_shannon.to_csv(outdir + "../shannon.tsv")


    #--------------------------------------------------------------------
    #Sensitivity, specificity, precision and FPRA
    confs = [confusion_matrix(sat[smpls].to_numpy(), sat_ref[smpls].to_numpy(), samples=smpls) for sat in matched_sats]
    superconf = pd.concat([conf.stack().rename(name) for name, conf in zip(names, confs)], axis=1)
    superconf.to_csv(outdir + "../superconf_{}_{}.tsv".format(simu, space), sep="\t")

    if False:
        #saving detailed information on false positives and false negatives
        for fname, conf, sat in zip(fnames, confs, matched_sats):
            print(fname)
            summary = []
            for smpl, FPRA in conf["FPRA"].sort_values(ascending=False).items():
                x, y = sat_ref[smpl], sat[smpl]
                FPs = y[(x==0)&(y>0)].sort_values(ascending=False)
                # summary.extend([[name, smpl, sp, fpra, 0.] for sp, fpra in FPs.items()])
                summary.extend([[fname, smpl, sp, fpra, 0.] for sp, fpra in FPs.items()])
                FNs = x[(x>0)&(y==0)].sort_values(ascending=False)
                # summary.extend([[name, smpl, sp, 0., fnra] for sp, fnra in FNs.items()])
                summary.extend([[fname, smpl, sp, 0., fnra] for sp, fnra in FNs.items()])

            # summary = pd.DataFrame(summary, columns=["tool", "sample", "feature", "FPRA", "FNRA"])
            summary = pd.DataFrame(summary, columns=["tool", "sample", "feature", "FPRA", "FNRA"])
            summary.to_csv(outdir + "../FPRA_summary_{}.tsv".format(fname), sep="\t", index=False)


    # #--------------------------------------------------------------------
    # kept = [pd.read_csv(indir + "{}_kept.tsv".format(name), sep="\t", index_col=0) for name in ["reference"] + fnames]

    # names1 = ["reference"] + names
    # lost_abundance = pd.concat([df.iloc[:,-1] for df in kept], axis=1).clip(lower=10**-15)
    # lost_abundance.columns = names1
    # lost_clades = pd.concat([df.iloc[:,-2] for df in kept], axis=1)
    # lost_clades.columns = names1

    # lost_abundance.to_csv(outdir + "../lost_abundance.tsv", sep="\t")
    # lost_clades.to_csv(outdir + "../lost_clades.tsv", sep="\t")
