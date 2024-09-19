#!/usr/bin/python3.6

from __future__ import division
import numpy as np
import pandas as pd
from scipy.stats import entropy
import os

simuls = ["refMet4", "refKrak"]
simu_labs = {"refMet4": "refMet4", "refKrak": "refKrak"}

if __name__ == "__main__":

    """
    Raw richness of simulation data
    """

    outdir = '../analyses/raw_richness/'
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    
    sats = [pd.read_csv(f"../data/{simu}_raw/reference.tsv", sep="\t", index_col=0).fillna(0.) for simu in simuls]
    richnesses = [(sat > 0).sum(axis=0) for sat in sats]
    shannons = [sat.apply(lambda x: entropy(x, base=2.), axis=0) for sat in sats]


    rich = pd.concat(richnesses, keys=[simu_labs[simu] for simu in simuls]).reset_index()
    rich.columns=['Simulation', 'Sample', 'Richness']
    shan = pd.concat(shannons, keys=[simu_labs[simu] for simu in simuls]).reset_index()
    shan.columns=['Simulation', 'Sample', 'Shannon']

    rich.to_csv(outdir + 'raw_richness.tsv', sep='\t', index=False)
    shan.to_csv(outdir + 'raw_shannon.tsv', sep='\t', index=False)


    corresp_ref = pd.read_csv("../data/taxo_indexes/reference_to_uhgg.tsv", sep="\t")
    idx = corresp_ref[corresp_ref.Lineage.str.contains(';g__;s__')].index
    corresp_ref['Lineage'] = corresp_ref.apply(lambda x: x['Lineage'] + x['reference'] if x['Lineage'].endswith(';s__') else x['Lineage'], axis=1)
    corresp_ref['Lineage'] = corresp_ref.apply(lambda x: x['Lineage'].replace('__;', '__' + x['reference'] + ';'), axis=1)
    corresp_ref.set_index('reference', inplace=True)


    corresp_ref['Phylum'] = corresp_ref.Lineage.str.extract(';p__([A-Za-z]*)[_;]')

    smpls = sats[0].columns.intersection(sats[1].columns)
    phylum_abundance = [sat[smpls].join(corresp_ref['Phylum'], how='left').groupby('Phylum').sum() for sat in sats]
    phylum_presence = [(sat[smpls] > 0).astype(np.float64).join(corresp_ref['Phylum'], how='left').groupby('Phylum').sum() for sat in sats]
    pd.concat(phylum_abundance, axis=0, keys=simuls).to_csv(outdir + 'phylum_abundance.tsv', sep='\t')
    pd.concat(phylum_presence, axis=0, keys=simuls).to_csv(outdir + 'phylum_presence.tsv', sep='\t')
