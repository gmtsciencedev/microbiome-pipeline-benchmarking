#!/usr/bin/python3.6

from __future__ import division
import numpy as np
import pandas as pd
import sys, os

def parse_corresp(corresp):
    dM = []
    for motu, row in corresp.squeeze().dropna().items():
        msps = row.split(",")
        if "msp" in msps:
            print(motu, msps)
        if len(msps) == 1:
            df = pd.DataFrame([[motu, msps[0], 1.]], columns=["mOTU", "msp_name", "abundance"])
            dM.append(df)
        else:
            df = pd.DataFrame([[motu] + msp.split(":") for msp in msps], columns=["mOTU", "msp_name", "abundance"])
            dM.append(df)

    dM = pd.concat(dM).astype({"abundance": float})
    return pd.pivot_table(dM, index="msp_name", columns="mOTU", values="abundance", fill_value=0.) 

def lost_stat(sat, idx_kept):
    sat_short = sat.loc[idx_kept]
    df_kept = pd.concat((sat_short.sum(), (sat_short>0).sum(), (sat > 0).sum()), axis=1)
    df_kept.columns = ["Kept_abundance", "Kept_mOTUs", "Raw_mOTUs"]
    df_kept["Lost_mOTUs"] = df_kept["Raw_mOTUs"] - df_kept["Kept_mOTUs"]
    df_kept["Lost_abundance"] = 1. - df_kept["Kept_abundance"]
    return df_kept

simu = "refBioms"#simuCRC2k#simuCRC2b#simuCRC

if __name__ == "__main__":

    """
    Projecting results of various pipelines to gtdb-r207 space
    """
    sats = []
    kept = []
    names = ["reference", "motus3", "metaphlan4", "metaphlan3", "kraken2", "biomscope"]

    refdir = "../data/taxo_indexes/"
    indir = "../data/{}_raw/".format(simu) 
    outdir = "../data/{}_to_igc2/".format(simu)
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    # ---------------------------------------------------------------------------
    #Reference sat
    #Downloading reference data
    SAT0 = pd.read_csv(indir + "reference.tsv", sep="\t", index_col=0).fillna(0.)
    smpls = SAT0.columns.tolist()
    SAT0 = SAT0.divide(SAT0.sum(), axis=1) #normalizing SAT
    # SAT0.to_csv(outdir + "reference_raw_norm.sat", sep="\t") #properly labeled raw sat 

    corresp_ref = pd.read_csv(refdir + "reference_to_igc2_final.tsv", sep="\t", index_col=0)
    dM = parse_corresp(corresp_ref["msp"])

    idx = SAT0.index.difference(dM.columns)
    if idx.size > 0:
        pd.DataFrame([], index=idx).to_csv(outdir + "missing_mgnify_to_igc2.tsv")

    cols = dM.columns.intersection(SAT0.index)
    sat0 = dM[cols].dot(SAT0.loc[cols])
    sat0 = sat0.divide(sat0.sum(), axis=1)
    # sat0.to_csv(outdir + "reference.sat", sep="\t")
    
    #statistics on lost species and clades
    df_kept = lost_stat(SAT0, cols)
    df_kept.to_csv(outdir + "reference_kept.tsv", sep="\t")

    sats.append(sat0)
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
    # SAT.to_csv(outdir + "motus3_raw_norm.sat", sep="\t") #properly labeled raw sat 

    corresp_motus = pd.read_csv(refdir + "motus3_to_igc2_final.tsv", sep="\t", index_col=0)
    dM = parse_corresp(corresp_motus["msp"])

    idx = SAT.index.difference(dM.columns)
    if idx.size > 0:
        pd.DataFrame([], index=idx).to_csv(outdir + "missing_motus3_to_igc2.tsv")

    cols = dM.columns.intersection(SAT.index)
    sat = dM[cols].dot(SAT.loc[cols])
    sat = sat.divide(sat.sum(), axis=1)
    # sat.to_csv(outdir + "motus3.sat", sep="\t")
    
    #statistics on lost species and clades
    df_kept = lost_stat(SAT, cols)
    df_kept.to_csv(outdir + "motus3_kept.tsv", sep="\t")

    sats.append(sat)
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
    # SAT.to_csv(outdir + "metaphlan4_raw_norm.sat", sep="\t") #properly labeled raw sat 

    corresp_metaphlan4 = pd.read_csv(refdir + "metaphlan4_to_igc2_final.tsv", sep="\t", index_col=0)
    dM = parse_corresp(corresp_metaphlan4["msp"])

    idx = SAT.index.difference(dM.columns)
    if idx.size > 0:
        pd.DataFrame([], index=idx).to_csv(outdir + "missing_metaphlan4_to_igc2.tsv")

    cols = dM.columns.intersection(SAT.index)
    sat = dM[cols].dot(SAT.loc[cols])
    sat = sat.divide(sat.sum(), axis=1)
    # sat.to_csv(outdir + "metaphlan4.sat", sep="\t")
    
    #statistics on lost species and clades
    df_kept = lost_stat(SAT, cols)
    df_kept.to_csv(outdir + "metaphlan4_kept.tsv", sep="\t")

    sats.append(sat)
    kept.append(df_kept)

    # ---------------------------------------------------------------------------
    #metaphlan3

    #Loading pipeline data
    SAT = pd.read_csv(indir + "metaphlan3.tsv", sep="\t", skiprows=[0], index_col="clade_name")
    SAT.rename(columns={col : col.split(".")[0] for col in SAT.columns}, inplace=True)
    SAT = SAT[smpls]
    SAT = SAT[SAT.index.str.contains("|s__", regex=False)]
    SAT = SAT.divide(SAT.sum(), axis=1) #normalizing SAT
    # SAT.to_csv(outdir + "metaphlan3_raw_norm.sat", sep="\t")

    corresp_metaphlan3 = pd.read_csv(refdir + "metaphlan3_to_igc2_final.tsv", sep="\t", index_col=0)
    dM = parse_corresp(corresp_metaphlan3["msp"])

    idx = SAT.index.difference(dM.columns)
    if idx.size > 0:
        pd.DataFrame([], index=idx).to_csv(outdir + "missing_metaphlan3_to_igc2.tsv")

    cols = dM.columns.intersection(SAT.index)
    sat = dM[cols].dot(SAT.loc[cols])
    sat = sat.divide(sat.sum(), axis=1)
    # sat.to_csv(outdir + "metaphlan3.sat", sep="\t")
    
    #statistics on lost species and clades
    df_kept = lost_stat(SAT, cols)
    df_kept.to_csv(outdir + "metaphlan3_kept.tsv", sep="\t")

    sats.append(sat)
    kept.append(df_kept)

    # ---------------------------------------------------------------------------
    #kraken2+bracken

    #Loading pipeline data
    SAT = pd.read_csv(indir + "kraken2.tsv", sep="\t", index_col="gtdb_taxonomy")
    SAT = SAT.divide(SAT.sum(axis=0), axis=1)
    SAT = SAT[smpls]
    SAT = SAT[SAT.sum(axis=1)>0]
    # SAT.to_csv(outdir + "kraken2_raw_norm.sat", sep="\t")

    corresp_kraken = pd.read_csv(refdir + "kraken2_to_igc2_final.tsv", sep="\t", index_col=0)
    dM = parse_corresp(corresp_kraken["msp"])

    idx = SAT.index.difference(dM.columns)
    if idx.size > 0:
        pd.DataFrame([], index=idx).to_csv(outdir + "missing_kraken2_to_igc2.tsv")

    cols = dM.columns.intersection(SAT.index)
    sat = dM[cols].dot(SAT.loc[cols])
    sat = sat.divide(sat.sum(), axis=1)
    # sat.to_csv(outdir + "kraken2.sat", sep="\t")

    #statistics on lost species and clades
    df_kept = lost_stat(SAT, cols)
    df_kept.to_csv(outdir + "kraken2_kept.tsv", sep="\t")

    sats.append(sat)
    kept.append(df_kept)

    # ---------------------------------------------------------------------------
    # biomscope (already in msp space)

    #Loading pipeline data
    SAT = pd.read_csv(indir + "biomscope.tsv", sep="\t", index_col=0)
    # SAT = pd.read_csv(indir + "../hyper_to_msp/biomscope_core6_trim0.2_allsample0_msp_matched.sat", sep="\t", index_col=0)
    # SAT = pd.read_csv(indir + "all_samples.sat", sep="\t", index_col=0)
    SAT = SAT[smpls]
    SAT = SAT.divide(SAT.sum(), axis=1) #normalizing SAT
    SAT = SAT[SAT.sum(axis=1)>0]
    # SAT.to_csv(outdir + "biomscope_raw_norm.sat", sep="\t")

    # SAT.to_csv(outdir + "biomscope.sat", sep="\t")
    #statistics on lost species and clades
    df_kept = lost_stat(SAT, SAT.index)
    df_kept.to_csv(outdir + "biomscope_kept.tsv", sep="\t")

    sats.append(SAT)
    kept.append(df_kept)


    msp_taxo = pd.read_csv(refdir + "IGC2.1989MSPs.taxo.tsv", sep="\t", index_col="msp_name")
    idx_msp = msp_taxo.index
    if not all([sat.index.difference(idx_msp).empty for sat in sats]):
        print("WARNING: some index values are not in msp_catalogue")

    # ---------------------------------------------------------------------------
    #Combining data for lost abundance and species
    lost_abundance = pd.concat([df.iloc[:,-1] for df in kept], axis=1).clip(lower=10**-15)
    lost_abundance.columns = names
    lost_clades = pd.concat([df.iloc[:,-2] for df in kept], axis=1)
    lost_clades.columns = names

    lost_abundance.to_csv(outdir + "lost_abundance.tsv", sep="\t")
    lost_clades.to_csv(outdir + "lost_clades.tsv", sep="\t")


    #Combining abundance tables
    #Matching sats to have the same index
    idx = pd.Index([])
    for sat in sats:
        idx = idx.union(sat.index)

    idx = idx.sort_values()
    idx.name = "msp_name"
    # sat_ref = sats_short[0].reindex(idx, fill_value=0.)
    matched_sats = [sat.reindex(idx, fill_value=0.) for sat in sats]

    # sat_ref.to_csv(outdir + "reference_gtdb207_matched.tsv", sep="\t")
    for name, sat in zip(names, matched_sats):
        sat.to_csv(outdir + "{}_matched.sat".format(name), sep="\t") 

    #saving data as a single table
    supersat = pd.concat([sat.stack().rename(name) for name, sat in zip(names, matched_sats)], axis=1)
    supersat.to_csv(outdir + "supersat_{}_igc2.sat".format(simu), sep="\t")
