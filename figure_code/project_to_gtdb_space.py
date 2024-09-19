#!/usr/bin/python3.6

from __future__ import division
import numpy as np
import pandas as pd
import sys, os, re

def lost_stat(sat, idx_kept):
    sat_short = sat.loc[idx_kept]
    df_kept = pd.concat((sat_short.sum(), (sat_short>0).sum(), (sat > 0).sum()), axis=1)
    df_kept.columns = ["Kept_abundance", "Kept_mOTUs", "Raw_mOTUs"]
    df_kept["Lost_mOTUs"] = df_kept["Raw_mOTUs"] - df_kept["Kept_mOTUs"]
    df_kept["Lost_abundance"] = 1. - df_kept["Kept_abundance"]
    df_kept["Lost_mOTUs_frac"] = df_kept["Lost_mOTUs"] / df_kept["Raw_mOTUs"]
    return df_kept

simu = "refMet4"#refKrak#refMet4

if __name__ == "__main__":

    """
    Projecting results of various pipelines to gtdb-r207 space
    """
    sats = []
    kept = []
    names = ["motus3", "metaphlan4", "metaphlan3", "kraken2"]

    refdir = "../data/taxo_indexes/"
    indir = "../data/{}_raw/".format(simu) 
    outdir = "../data/{}_to_gtdb207/".format(simu)
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    # ---------------------------------------------------------------------------
    # Full gtdb207 taxonomy
    gtdb207_bac = pd.read_csv(refdir + "bac120_taxonomy_r207.tsv", sep="\t", header=None, names=["accession", "GTDB-R207"])
    gtdb207_ar = pd.read_csv(refdir + "ar53_taxonomy_r207.tsv", sep="\t", header=None, names=["accession", "GTDB-R207"])
    gtdb207 = pd.concat((gtdb207_bac, gtdb207_ar))#.drop_duplicates()
    idx207 = pd.Index(gtdb207["GTDB-R207"]).drop_duplicates()

    # ---------------------------------------------------------------------------
    #Reference sat
    #Downloading reference data
    SAT0 = pd.read_csv(indir + "reference.tsv", sep="\t", index_col=0).fillna(0.)
    smpls = SAT0.columns.tolist()
    SAT0 = SAT0.divide(SAT0.sum(), axis=1) #normalizing SAT

    corresp_ref = pd.read_csv(refdir + "reference_to_gtdb_final.tsv", sep="\t", index_col="mgnify")
    
    idx = SAT0.index.difference(corresp_ref.index)
    if idx.size > 0:
        pd.DataFrame([], index=idx).to_csv(outdir + "missing_mgnify_to_gtdb.tsv")

    #projecting sat onto gtdb207
    sat0 = SAT0.join(corresp_ref, how="left").groupby("gtdb").sum()
    sat0.index.name = "GTDB-R207"
    sat0 = sat0.divide(sat0.sum(), axis=1)

    df_kept = lost_stat(SAT0, corresp_ref.index.intersection(SAT0.index))
    kept.append(df_kept)

    # ---------------------------------------------------------------------------
    #motus3

    #Loading pipeline data
    SAT = pd.read_csv(indir + "motus3.tsv", sep="\t", index_col="#consensus_taxonomy")
    SAT = SAT[smpls]
    SAT = SAT[SAT.sum(axis=1)>0]
    SAT["mOTU"] = SAT.index.str.extract("\[([ref|meta|ext].*)\]", expand=False)
    SAT.dropna(subset=["mOTU"], inplace=True)
    SAT.set_index("mOTU", inplace=True)
    SAT = SAT.divide(SAT.sum(), axis=1) #normalizing SAT

    corresp_motus = pd.read_csv(refdir + "motus3_to_gtdb_final.tsv", sep="\t", index_col="mOTU")

    sat = SAT.join(corresp_motus, how="left").groupby("GTDB-R207").sum()
    sat = sat.divide(sat.sum(), axis=1)
    sats.append(sat)

    df_kept = lost_stat(SAT, corresp_motus.index.intersection(SAT.index))
    kept.append(df_kept)

    # ---------------------------------------------------------------------------
    #metaphlan4

    #Loading pipeline data
    SAT = pd.read_csv(indir + "metaphlan4.tsv", sep="\t", skiprows=[0], index_col=0).fillna(0.)
    SAT.rename(columns={col : col.split(".")[0] for col in SAT.columns}, inplace=True)
    SAT = SAT[smpls]
    SAT = SAT[SAT.index.str.contains("|t__", regex=False)]
    SAT["metaphlan"] = SAT.index.str.extract("t__(.*)", expand=False)
    SAT.set_index("metaphlan", inplace=True)
    SAT = SAT.divide(SAT.sum(), axis=1) #normalizing SAT

    corresp_metaphlan4 = pd.read_csv(refdir + "mpa42gtdb.tsv", sep="\t", index_col=0)
    corresp_metaphlan4_1 = pd.read_csv(refdir + "mpa42gtdb3.tsv", sep="\t", index_col=0)
    corresp_metaphlan4 = pd.concat((corresp_metaphlan4, corresp_metaphlan4_1))
    
    idx = SAT.index.difference(corresp_metaphlan4.index)
    if idx.size > 0:
        pd.DataFrame([], index=idx).to_csv(outdir + "missing_metaphlan4_to_gtdb.tsv")

    sat = SAT.join(corresp_metaphlan4, how="left").groupby("gtdbTaxonomy").sum()
    sat.index.name = "GTDB-R207"
    sat = sat.divide(sat.sum(), axis=1)
    sats.append(sat)

    df_kept = lost_stat(SAT, corresp_metaphlan4.index.intersection(SAT.index))
    kept.append(df_kept)

    # ---------------------------------------------------------------------------
    #metaphlan3

    #Loading pipeline data
    SAT = pd.read_csv(indir + "metaphlan3.tsv", sep="\t", skiprows=[0], index_col="clade_name")
    SAT.rename(columns={col : col.split(".")[0] for col in SAT.columns}, inplace=True)
    SAT = SAT[smpls]
    SAT = SAT[SAT.index.str.contains("|s__", regex=False)]
    SAT = SAT.divide(SAT.sum(), axis=1) #normalizing SAT

    corresp_metaphlan3 = pd.read_csv(refdir + "metaphlan3_to_gtdb_final.tsv", sep="\t", index_col="clade")
    
    idx = SAT.index.difference(corresp_metaphlan3.index)
    if idx.size > 0:
        pd.DataFrame([], index=idx).to_csv(outdir + "missing_metaphlan3_to_gtdb.tsv")

    sat = SAT.join(corresp_metaphlan3, how="left").groupby("GTDB-R207").sum()
    sat = sat.divide(sat.sum(), axis=1)
    sats.append(sat)

    df_kept = lost_stat(SAT, corresp_metaphlan3.index.intersection(SAT.index))
    kept.append(df_kept)

    # ---------------------------------------------------------------------------
    #kraken2+bracken (already in gtdb207)

    #Loading pipeline data
    SAT = pd.read_csv(indir + "kraken2.tsv", sep="\t", index_col="gtdb_taxonomy")
    SAT = SAT.divide(SAT.sum(axis=0), axis=1)
    SAT = SAT[smpls]
    SAT = SAT[SAT.sum(axis=1)>0]

    SAT.index.name = "GTDB-R207"
    sats.append(SAT)

    df_kept = lost_stat(SAT, SAT.index)
    kept.append(df_kept)

    # ---------------------------------------------------------------------------
    #Combining data for lost abundance and species
    names1 = ["reference"] + names

    lost_summary = pd.concat([df.stack() for df in kept], axis=1)
    lost_summary.columns = names1
    lost_summary.to_csv(outdir + "lost_summary.tsv", sep="\t")

    #Combining abundance tables
    #Matching sats to have the same index
    idx = sat0.index.copy()
    for sat in sats:
        idx = idx.union(sat.index)

    idx = idx.sort_values()
    idx.name = "GTDB-R207"
    sat_ref = sat0.reindex(idx, fill_value=0.)
    matched_sats = [sat.reindex(idx, fill_value=0.) for sat in sats]

    #saving data as a single table
    supersat = pd.concat([sat_ref.stack().rename("reference")] + [sat.stack().rename(name) for name, sat in zip(names, matched_sats)], axis=1)
    supersat.to_csv(outdir + "supersat_{}_gtdb207.sat".format(simu), sep="\t")
