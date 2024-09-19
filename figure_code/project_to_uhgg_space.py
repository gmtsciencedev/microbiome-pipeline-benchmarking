#!/usr/bin/python3.6

from __future__ import division
import numpy as np
import pandas as pd
import sys, os, re

from project_to_gtdb_space import lost_stat

def parse_corresp(corresp):
    dM = []
    for motu, row in corresp.squeeze().dropna().items():
        msps = row.split(",")
        # if "msp" in msps:
        #     print(motu, msps)
        n = len(msps)
        if n == 1:
            df = pd.DataFrame([[motu, msps[0], 1.]], columns=["mOTU", "msp_name", "abundance"])
            dM.append(df)
        else:
            # df = pd.DataFrame([[motu] + msp.split(":") for msp in msps], columns=["mOTU", "msp_name", "abundance"])
            df = pd.DataFrame([[motu, msp.strip(), 1./n] for msp in msps], columns=["mOTU", "msp_name", "abundance"])
            dM.append(df)

    dM = pd.concat(dM).astype({"abundance": float})
    return pd.pivot_table(dM, index="msp_name", columns="mOTU", values="abundance", fill_value=0.) 

simu = "refKrak"#refKrak#refMet4

if __name__ == "__main__":

    """
    Projecting results of various pipelines to UHGG space
    """
    sats = []
    kept = []
    names = ["motus3", "metaphlan4", "metaphlan3", "kraken2"]

    refdir = "../data/taxo_indexes/"
    indir = "../data/{}_raw/".format(simu) 
    outdir = "../data/{}_to_uhgg/".format(simu)
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    # ---------------------------------------------------------------------------
    #Reference sat
    #Downloading reference data
    SAT0 = pd.read_csv(indir + "reference.tsv", sep="\t", index_col=0).fillna(0.)
    smpls = SAT0.columns.tolist()
    SAT0 = SAT0.divide(SAT0.sum(), axis=1) #normalizing SAT

    corresp_ref = pd.read_csv(refdir + "reference_to_uhgg.tsv", sep="\t")
    idx = corresp_ref[corresp_ref.Lineage.str.contains(';g__;s__')].index
    corresp_ref['Lineage'] = corresp_ref.apply(lambda x: x['Lineage'] + x['reference'] if x['Lineage'].endswith(';s__') else x['Lineage'], axis=1)
    corresp_ref['Lineage'] = corresp_ref.apply(lambda x: x['Lineage'].replace('__;', '__' + x['reference'] + ';'), axis=1)
    corresp_ref.set_index('reference', inplace=True)

    idx = SAT0.index.difference(corresp_ref.index)
    if idx.size > 0:
        pd.DataFrame([], index=idx).to_csv(outdir + "missing_mgnify_to_uhgg.tsv")

    #projecting sat onto uhgg
    sat0 = SAT0.join(corresp_ref, how="left").groupby("Lineage").sum()
    sat0.index.name = "UHGG"
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
    
    corresp_motus = pd.read_csv(refdir + "motus2uhgg.tsv", sep="\t", index_col=0)
    corresp_motus1 = pd.read_csv(refdir + "motus2uhgg-missing.tsv", sep="\t", index_col=0)
    corresp_motus = pd.concat((corresp_motus, corresp_motus1))
    corresp_motus['UHGG_taxonomy'] = corresp_motus['UHGG_taxonomy'].str.replace('r1__', 'd__')
    corresp_motus = corresp_motus[~corresp_motus['UHGG_taxonomy'].str.contains('unclassified')]

    idx = SAT.index.difference(corresp_motus.index)
    if idx.size > 0:
        pd.DataFrame([], index=idx).to_csv(outdir + "missing_mOTU_to_uhgg.tsv")

    sat = SAT.join(corresp_motus, how="left").groupby("UHGG_taxonomy").sum()
    sat.index.name = "UHGG"
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

    corresp_metaphlan4 = pd.read_csv(refdir + "mpa42uhgg.tsv", sep="\t", index_col=0)
    corresp_metaphlan4_1 = pd.read_csv(refdir + "mpa42uhgg2.tsv", sep="\t", index_col=0, header=0, names=corresp_metaphlan4.columns)
    corresp_metaphlan4 = pd.concat((corresp_metaphlan4, corresp_metaphlan4_1))
    corresp_metaphlan4['UHGG_taxonomy'] = corresp_metaphlan4['UHGG_taxonomy'].str.replace('r1__', 'd__')
    corresp_metaphlan4 = corresp_metaphlan4[~corresp_metaphlan4['UHGG_taxonomy'].str.contains('unclassified')]

    idx = SAT.index.difference(corresp_metaphlan4.index)
    if idx.size > 0:
        pd.DataFrame([], index=idx).to_csv(outdir + "missing_metaphlan4_to_uhgg.tsv")

    sat = SAT.join(corresp_metaphlan4, how="left").groupby("UHGG_taxonomy").sum()
    sat.index.name = "UHGG"
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

    corresp_metaphlan3 = pd.read_csv(refdir + "mpa3_to_uhgg_final.tsv", sep="\t", index_col=0)
    corresp_metaphlan3['UHGG_taxonomy'] = corresp_metaphlan3['UHGG_taxonomy'].str.replace('r1__', 'd__')
    corresp_metaphlan3 = corresp_metaphlan3[~corresp_metaphlan3['UHGG_taxonomy'].str.contains('unclassified')]
    dM = parse_corresp(corresp_metaphlan3)

    idx = SAT.index.difference(dM.columns)
    if idx.size > 0:
        pd.DataFrame([], index=idx).to_csv(outdir + "missing_metaphlan3_to_uhgg.tsv")

    cols = dM.columns.intersection(SAT.index)
    sat = dM[cols].dot(SAT.loc[cols]).groupby(dM.index.name).sum()
    sat.index.name = "UHGG"
    sat = sat.divide(sat.sum(), axis=1)
    sats.append(sat)

    df_kept = lost_stat(SAT, cols)
    kept.append(df_kept)

    # ---------------------------------------------------------------------------
    #kraken2+bracken

    #Loading pipeline data
    SAT = pd.read_csv(indir + "kraken2.tsv", sep="\t", index_col="gtdb_taxonomy")
    SAT = SAT.divide(SAT.sum(axis=0), axis=1)
    SAT = SAT[smpls]
    SAT = SAT[SAT.sum(axis=1)>0]

    corresp_kraken = pd.read_csv(refdir + "kraken2uhgg.tsv", sep="\t")
    corresp_kraken1 = pd.read_csv(refdir + "kraken_gtdb2uhgg_new.tsv", sep="\t")
    corresp_kraken1.rename(columns={'uhgg': 'UHGG_taxonomy'}, inplace=True)
    corresp_kraken = pd.concat((corresp_kraken, corresp_kraken1))
    corresp_kraken['UHGG_taxonomy'] = corresp_kraken['UHGG_taxonomy'].str.replace('r1__', 'd__')
    corresp_kraken['taxo'] = corresp_kraken['taxo'].str.replace('r1__', 'd__')
    corresp_kraken = corresp_kraken[~corresp_kraken['UHGG_taxonomy'].str.contains('unclassified')]
    corresp_kraken.drop_duplicates(inplace=True)

    dM = []
    for taxo, group in corresp_kraken.groupby('taxo'):
        group['abundance'] = 1. / group.shape[0]
        dM.append(group)

    dM = pd.concat(dM).astype({"abundance": float})
    dM = pd.pivot_table(dM, index="UHGG_taxonomy", columns="taxo", values="abundance", fill_value=0.) 

    idx = SAT.index.difference(dM.columns)
    if idx.size > 0:
        pd.DataFrame([], index=idx).to_csv(outdir + "missing_kraken_to_uhgg.tsv")

    cols = dM.columns.intersection(SAT.index)
    sat = dM[cols].dot(SAT.loc[cols]).groupby(dM.index.name).sum()
    sat.index.name = "UHGG"
    sat = sat.divide(sat.sum(), axis=1)
    sats.append(sat)

    df_kept = lost_stat(SAT, cols)
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
    idx.name = "UHGG"
    sat_ref = sat0.reindex(idx, fill_value=0.)
    matched_sats = [sat.reindex(idx, fill_value=0.) for sat in sats]

    #saving data as a single table
    supersat = pd.concat([sat_ref.stack().rename("reference")] + [sat.stack().rename(name) for name, sat in zip(names, matched_sats)], axis=1)
    supersat.to_csv(outdir + "supersat_{}_uhgg.sat".format(simu), sep="\t")
