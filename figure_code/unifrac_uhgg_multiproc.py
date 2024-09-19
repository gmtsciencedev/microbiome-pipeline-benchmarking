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
simu = "refMet4"#"refKrak"#"refMet4"
space = "uhgg"
names = ["kraken2 + bracken", "metaphlan3", "motus3", "metaphlan4"]
fnames = ["kraken2", "metaphlan3", "motus3", "metaphlan4"]


indir = "../data/{}_to_{}/".format(simu, space)
supersat = pd.read_csv(indir + "supersat_{}_{}.sat".format(simu, space), sep="\t", index_col=(0,1))
matched_sats = [supersat[name].unstack() for name in fnames]
sat_ref = supersat["reference"].unstack()
smpls = sat_ref.columns#.drop(["Species_rep", "mgnify_taxo"]).tolist()
nsmpls = len(smpls)

refdir = "../data/taxo_indexes/"

if __name__ == "__main__":
    """
    Calculating UniFrac distance matrices in UHGG space
    """
    outdir = "../analyses/article_results_{}_uhgg/unifrac/".format(simu)
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    print("Preparing phylogenetic tree")
    t0 = time()
    uhgg_bac = pd.read_csv(refdir + "genomes-all_metadata.tsv", sep="\t", usecols=['Genome', 'Lineage']).rename(columns={'Lineage': 'UHGG'}).set_index('UHGG')
    tree_bac = Phylo.read(refdir + "bac120_iqtree.nwk", format="newick")
    acc_bac = pd.Index([t.name for t in tree_bac.get_terminals()])

    df_bac = uhgg_bac[uhgg_bac.index.isin(sat_ref.index)&
                uhgg_bac.Genome.isin(acc_bac)]
    idx_MGYG = sat_ref.index[sat_ref.index.str.contains(';s__MGYG')]
    df_bac1 = idx_MGYG.str.extract(';s__(MGYG.*)')
    df_bac1.index=idx_MGYG
    df_bac1.columns = ['Genome']
    df_bac = pd.concat((df_bac, df_bac1))
    df_bac = df_bac[df_bac.Genome.isin(acc_bac)]
    df_bac = df_bac[~df_bac.index.duplicated()]

    ids = df_bac.squeeze().drop_duplicates().tolist()
    for genome in acc_bac:
        if genome not in ids:
            tree_bac.prune(genome)
    
    tree_bac.root.name="root"
    Phylo.write(tree_bac, outdir + "tree_bac.nwk", "newick")
    df_bac.to_csv(outdir + "tree_id_bac.tsv", sep="\t")
    t1 = time()
    print(t1-t0)    

    #--------------------------------------------------------------------------
    #Creating unifrac matrices

    #load tree in scikit-bio format
    tree_bac = read(outdir + "tree_bac.nwk", format="newick", into=TreeNode, convert_underscores=False)
    tree_bac = tree_bac.root_at("root")
    sat_bac_ref = sat_ref.join(df_bac, how='right').groupby('Genome').sum().loc[ids]

    sats_bac = [sat.join(df_bac, how='right').groupby('Genome').sum().loc[ids] for sat in matched_sats]

    def calculate_D_uu_D_wu(fname, sat):
        D_uu, D_wu = build_unifrac_matrices(sat, ids, tree_bac, validate=False)
        D_uu.to_csv(outdir + "Duu_{}.tsv".format(fname), sep="\t")
        D_wu.to_csv(outdir + "Dwu_{}.tsv".format(fname), sep="\t")

    from multiprocessing import Process
    pp = []
    for fname, sat_bac in zip(['ref'] + fnames, [sat_bac_ref] + sats_bac):
        print(fname)
        p = Process(target=calculate_D_uu_D_wu, args=(fname, sat_bac))
        p.start()
        pp.append(p)
    
    [p.join() for p in pp]

    D_ref_uu = pd.read_csv(outdir + "Duu_ref.tsv", sep="\t", index_col=0)
    D_ref_wu = pd.read_csv(outdir + "Dwu_ref.tsv", sep="\t", index_col=0)

    DD_uu = []
    DD_wu = []
    for fname in fnames:
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
    dist_to_ref = [unifrac_dist_to_ref(sat_bac_ref, sat, ids, tree_bac, validate=False) for sat in sats_bac]
    dist_to_ref_uu = pd.concat([df["unweighted_unifrac"].rename(name) for name, df in zip(names, dist_to_ref)], axis=1)
    dist_to_ref_wu = pd.concat([df["weighted_unifrac"].rename(name) for name, df in zip(names, dist_to_ref)], axis=1)

    dist_to_ref_uu.to_csv(outdir + "dist_to_ref_uu.tsv", sep="\t")
    dist_to_ref_wu.to_csv(outdir + "dist_to_ref_wu.tsv", sep="\t")

