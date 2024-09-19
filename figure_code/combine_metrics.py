#!/usr/bin/python3.6

from __future__ import division
import numpy as np
import pandas as pd
import os

names = ["kraken2 + bracken", "metaphlan3", "motus3", "metaphlan4"]
simuls = ["refKrak", "refMet4"]
spaces = ["gtdb207", "uhgg"]

if __name__ == "__main__":
    """
    Combine metrics calculated for different simulations and in different projection space
    """

    outdir = "../analyses/combined_metrics/"
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    #--------------------------------------------------------------------------
    #Distance to reference
    dist_to_ref = []
    for sp in spaces:
        for simu in simuls:
            indir = "../analyses/article_results_{}_{}/braycurtis/".format(simu, sp)
            df = pd.read_csv(indir + "dist_to_ref.tsv", sep="\t", index_col=0)
            df[["space", "simu"]] = [sp, simu]
            dist_to_ref.append(df)

    dist_to_ref = pd.concat(dist_to_ref).set_index(["space", "simu"], append=True)
    dist_to_ref.to_csv(outdir + "dist_to_ref.tsv", sep="\t")

    dist_to_ref_wu = []
    for sp in spaces:
        for simu in simuls:
            indir = "../analyses/article_results_{}_{}/unifrac/".format(simu, sp)
            df = pd.read_csv(indir + "dist_to_ref_wu.tsv", sep="\t", index_col=0)
            df[["space", "simu"]] = [sp, simu]
            dist_to_ref_wu.append(df)

    dist_to_ref_wu = pd.concat(dist_to_ref_wu).set_index(["space", "simu"], append=True)
    dist_to_ref_wu.to_csv(outdir + "dist_to_ref_wu.tsv", sep="\t")

    #--------------------------------------------------------------------------
    #Lost abundance and clades
    # lost_abundance = []
    # lost_clades = []
    # for jsp, sp in enumerate(spaces):
    #     for simu in simuls:
    #         # indir = "../analyses/article_results_{}_{}/".format(simu, sp)
    #         indir = "../data/{}_to_{}/".format(simu, sp)
    #         lost_ab = pd.read_csv(indir + "lost_abundance.tsv", sep="\t", index_col=0)
    #         lost_ab[["space", "simu"]] = [sp, simu]
    #         lost_abundance.append(lost_ab)
            
    #         lost_cl = pd.read_csv(indir + "lost_clades.tsv", sep="\t", index_col=0)
    #         lost_cl[["space", "simu"]] = [sp, simu]
    #         lost_clades.append(lost_cl)

    # lost_abundance = pd.concat(lost_abundance).set_index(["space", "simu"], append=True)
    # lost_clades = pd.concat(lost_clades).set_index(["space", "simu"], append=True)
    # lost_abundance.to_csv(outdir + "lost_abundance.tsv", sep="\t")
    # lost_clades.to_csv(outdir + "lost_clades.tsv", sep="\t")

    lost_summary = []
    for jsp, sp in enumerate(spaces):
        for simu in simuls:
            indir = "../data/{}_to_{}/".format(simu, sp)
            lost_sum = pd.read_csv(indir + "lost_summary.tsv", sep="\t", index_col=(0,1))
            lost_sum[["space", "simu"]] = [sp, simu]
            lost_summary.append(lost_sum)
            
    lost_summary = pd.concat(lost_summary).set_index(["space", "simu"], append=True)
    lost_summary.index.rename("stat", level=1, inplace=True)
    lost_summary.index.rename("sample", level=0, inplace=True)
    lost_summary.to_csv(outdir + "lost_summary.tsv", sep="\t")
    
    #--------------------------------------------------------------------------
    #richness and shannon entropy            
    richness = []
    shannon = []
    for sp in spaces:
        for simu in simuls:
            indir = "../analyses/article_results_{}_{}/".format(simu, sp)
            rich = pd.read_csv(indir + "richness.tsv", sep=",", index_col=0)
            rich[["space", "simu"]] = [sp, simu]
            richness.append(rich)
            shan = pd.read_csv(indir + "shannon.tsv", sep=",", index_col=0)
            shan[["space", "simu"]] = [sp, simu]
            shannon.append(shan)

    richness = pd.concat(richness).set_index(["space", "simu"], append=True)
    shannon = pd.concat(shannon).set_index(["space", "simu"], append=True)
    richness.to_csv(outdir + "richness.tsv", sep="\t")
    shannon.to_csv(outdir + "shannon.tsv", sep="\t")

    #--------------------------------------------------------------------------
    #assembling confusion matrices
    superconf = []
    for sp in spaces:
        for simu in simuls:
            indir = "../analyses/article_results_{}_{}/".format(simu, sp)
            df = pd.read_csv(indir + "superconf_{}_{}.tsv".format(simu, sp), sep="\t", index_col=[0,1])
            df[["space", "simu"]] = [sp, simu]
            superconf.append(df)

    superconf = pd.concat(superconf).set_index(["space", "simu"], append=True)
    superconf.to_csv(outdir + "superconf.tsv", sep="\t")

    #--------------------------------------------------------------------------
    #assembling fpra/fnra summary
    fpra_summary = []
    for sp in spaces:
        for simu in simuls:
            indir = "../analyses/article_results_{}_{}/".format(simu, sp)
            df = pd.read_csv(indir + "FPRA_summary.tsv", sep="\t", index_col=[0,1])
            df[["space", "simu"]] = [sp, simu]
            fpra_summary.append(df)

    fpra_summary = pd.concat(fpra_summary).set_index(["space", "simu"], append=True)
    fpra_summary.to_csv(outdir + "FPRA_summary.tsv", sep="\t")

    #--------------------------------------------------------------------------
    #assembling tables
    richness_diff = richness[names].sub(richness["simulation"], axis=0)
    shannon_diff = shannon[names].sub(shannon["simulation"], axis=0)

    lost_summary.rename(columns={'kraken2': 'kraken2 + bracken', 'reference': 'simulation'}, inplace=True)
    data = [lost_summary.xs('Lost_abundance', level=1),
             lost_summary.xs('Lost_mOTUs', level=1), 
             richness, 
             shannon,
             richness_diff,
             shannon_diff,
             dist_to_ref, 
             dist_to_ref_wu,
             superconf.xs("TPR", level=1),
             superconf.xs("PPV", level=1),
             superconf.xs("FPRA", level=1),
             superconf.xs("FNRA", level=1)]

    data_names = ["lost_abundance",
                   "lost_clades",
                   "richness",
                   "shannon",
                   "richness_diff",
                   "shannon_diff",
                   "braycurtis",
                   "unifrac",
                   "TPR",
                   "PPV",
                   "FPRA",
                   "FNRA"]

    if not os.path.exists(outdir + "quarts/"):
        os.makedirs(outdir + "quarts/")

    def calc_quarts(df):
        df.columns.name="tool"
        out = pd.concat([df.groupby(["space", "simu"]).quantile(q=x).T.stack().stack() for x in [0.25, .5, 0.75]], axis=1)
        out.columns = ["q25", "q50", "q75"]
        out["iqr"] = out["q75"] - out["q25"]
        return out
    
    summary = []
    for df, name in zip(data, data_names):
        dg = calc_quarts(df)
        summary.append(dg)
        dg.to_csv(outdir + "quarts/{}.tsv".format(name), sep="\t")

    out = pd.concat([df["q50"] for df in summary[:4]], axis=1)
    out.columns = data_names[:4]
    out.to_csv(outdir + "median_basic.tsv", sep="\t")

    out1 = pd.concat([df["q50"] for df in summary[4:]], axis=1)
    out1.columns = data_names[4:]
    out1.to_csv(outdir + "median_measures.tsv", sep="\t")
