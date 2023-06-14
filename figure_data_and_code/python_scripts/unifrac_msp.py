#!/usr/bin/python3.6

from __future__ import division
import numpy as np
import pandas as pd
from skbio.diversity.beta import unweighted_unifrac, weighted_unifrac
from skbio import read
from skbio.tree import TreeNode
from Bio import Phylo
from time import time
import sys, os

from unifrac_aux import build_unifrac_matrices, unifrac_dist_to_ref
from metrics_aux import intermatrix_norms


#------------------------------------------------------------------------------
simu = "simuCRC"#"simuCRC2k"#"simuCRC2b"#"simuCRC"
space = "msp"
names = ["kraken2 + bracken", "metaphlan3", "motus3", "metaphlan4", "biomscope"]
fnames = ["kraken2", "metaphlan3", "motus3", "metaphlan4", "biomscope"]


indir = "../data/{}_to_{}/".format(simu, space)

supersat = pd.read_csv(indir + "supersat_{}_{}.sat".format(simu, space), sep="\t", index_col=(0,1))
matched_sats = [supersat[name].unstack() for name in fnames]
sat_ref = supersat["reference"].unstack()
smpls = sat_ref.columns
nsmpls = len(smpls)

# Full gtdb207 taxonomy
refdir = "../data/taxo_indexes/"
tree_path = refdir + "IGC2.1990MSPs.nwk"

if __name__ == "__main__":
    """
    Calculating UniFrac distance matrices in IGC2/MSP space
    """

    outdir = "../analyses/article_results_{}_msp/unifrac/".format(simu)
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    if True:
        print("Preparing phylogenetic tree")
        t0 = time()
        tree_bac = Phylo.read(tree_path, format="newick")
        acc_bac = pd.Index([t.name for t in tree_bac.get_terminals()])
        tips = sat_ref.index.intersection(pd.Index(acc_bac))
        bad_tips = pd.Index(acc_bac).difference(sat_ref.index)

        for msp in bad_tips:
            tree_bac.prune(msp)

        tree_bac.root.name="root"
        Phylo.write(tree_bac, outdir + "tree_bac.nwk", "newick")
        t1 = time()
        print(t1-t0)

    #--------------------------------------------------------------------------
    #Creating unifrac matrices
    #helpful vignette about rooting the tree https://www.biostars.org/p/467682/#google_vignette

    #load tree in scikit-bio format
    tree_bac = read(outdir + "tree_bac.nwk", format="newick", into=TreeNode, convert_underscores=False)
    tree_bac = tree_bac.root_at("root")
    tips = [t.name for t in tree_bac.tips()]
    sat_bac_ref = sat_ref.loc[tips]
    X = sat_ref.loc[tips].to_numpy()
    X = np.round(X / X[X>0].min())

    #unifrac matrices for reference data
    if True:
        print("reference")
        t0 = time()
        D_ref_uu, D_ref_wu = build_unifrac_matrices(sat_bac_ref.iloc[:,:nsmpls], tips, tree_bac, validate=False)
        D_ref_uu.to_csv(outdir + "Duu_ref.tsv", sep="\t")
        D_ref_wu.to_csv(outdir + "Dwu_ref.tsv", sep="\t")
        t1 = time()
        print(t1-t0)    
    else:
        D_ref_uu = pd.read_csv(outdir + "Duu_ref.tsv", sep="\t", index_col=0)
        D_ref_wu = pd.read_csv(outdir + "Dwu_ref.tsv", sep="\t", index_col=0)


    #calculate unifrac distances for pipelines' outputs
    sats_bac = [sat.loc[tips] for sat in matched_sats]
    DD_uu = []
    DD_wu = []
    if True:
        for fname, sat in zip(fnames, sats_bac):
            print(fname)
            t0 = time()
            D_uu, D_wu = build_unifrac_matrices(sat, tips, tree_bac, validate=False)
            D_uu.to_csv(outdir + "Duu_{}.tsv".format(fname), sep="\t")
            D_wu.to_csv(outdir + "Dwu_{}.tsv".format(fname), sep="\t")
            DD_uu.append(D_uu)
            DD_wu.append(D_wu)
            t1 = time()
            print(t1 - t0)
    else:
        for fname, sat in zip(fnames, matched_sats):
            D_uu = pd.read_csv(outdir + "Duu_{}.tsv".format(fname), sep="\t", index_col=0)
            D_wu = pd.read_csv(outdir + "Dwu_{}.tsv".format(fname), sep="\t", index_col=0)
            DD_uu.append(D_uu)
            DD_wu.append(D_wu)
    

    #--------------------------------------------------------------------
    #calculating intermatrix norms
    norm_ref_matched_uu = intermatrix_norms([D_ref_uu]*len(names), DD_uu, names, filepath=outdir +"norm_ref_matched_uu.tsv")
    norm_ref_matched_wu = intermatrix_norms([D_ref_wu]*len(names), DD_wu, names, filepath=outdir +"norm_ref_matched_wu.tsv")

    #--------------------------------------------------------------------
    #UniFrac distances between real sample and simulation
    if True:
        print("Distances between real sample and simulation")
        # dist_to_ref = [{name: unifrac_dist_to_ref(sat_unifrac_ref, sat, tips, tree)} for name, sat in zip(names, sats_unifrac)]
        dist_to_ref = [unifrac_dist_to_ref(sat_bac_ref, sat, tips, tree_bac, validate=False) for sat in sats_bac]
        dist_to_ref_uu = pd.concat([df["unweighted_unifrac"].rename(name) for name, df in zip(names, dist_to_ref)], axis=1)
        dist_to_ref_wu = pd.concat([df["weighted_unifrac"].rename(name) for name, df in zip(names, dist_to_ref)], axis=1)

        dist_to_ref_uu.to_csv(outdir + "dist_to_ref_uu.tsv", sep="\t")
        dist_to_ref_wu.to_csv(outdir + "dist_to_ref_wu.tsv", sep="\t")
    else:
        dist_to_ref_uu = pd.read_csv(outdir + "dist_to_ref_uu.tsv", sep="\t", index_col=0)
        dist_to_ref_wu = pd.read_csv(outdir + "dist_to_ref_wu.tsv", sep="\t", index_col=0)
