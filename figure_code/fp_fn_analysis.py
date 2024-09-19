#!/usr/bin/python3.6

from __future__ import division
import numpy as np
import pandas as pd
import os


names = ["kraken2 + bracken", "metaphlan3", "motus3", "metaphlan4"]
fnames = ["kraken2", "metaphlan3", "motus3", "metaphlan4"]
simuls = ["refKrak", "refMet4"]
spaces = ["gtdb207", "uhgg"]

simu_labs = {"refMet4": "refMet4", "refKrak": "refKrak"}
space_labs = {"gtdb207": "GTDB", "uhgg": "UHGG"}
tool_labs = {"kraken2 + bracken": "Krak/Brac", 
             "kraken2" : "Krak/Brac",
             "metaphlan3": "MPA3", 
             "metaphlan4": "MPA4", 
             "motus3": "mOTUs3", 
             "biomscope": "BiomS", 
             "simulation" : "Ref",
             "reference" : "Ref"}

labs = ["{}, {}".format(simu_labs[simu], space_labs[sp]) for simu in simuls for sp in spaces]
idx_labs =pd.MultiIndex.from_product([simuls, spaces], names=["space", "simu"])

nn_smpls = np.arange(343)
nsats = len(names)
nlabs = len(labs)

if __name__ == "__main__":
    """
    Combine metrics calculated for different simulations and in different projection space
    """

    nmin = 25
    outdir = '../analyses/confusion_pairs/'
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    if not os.path.exists(outdir + 'FPRA_summary_all.tsv'):
        fpra_summaries = []
        for i in idx_labs:
            fpra_summary = pd.read_csv("../analyses/article_results_{}_{}/FPRA_summary.tsv".format(*i), sep="\t", index_col=[0,1])
            fpra_summaries.append(fpra_summary)

        fpra_summaries = pd.concat(fpra_summaries, keys=idx_labs)
        fpra_summaries.to_csv(outdir + 'FPRA_summary_all.tsv', sep='\t')
    else:
        fpra_summaries = pd.read_csv(outdir + 'FPRA_summary_all.tsv', sep='\t', index_col=[0, 1, 2, 3])


    #compiling coincidence data
    if os.path.exists(outdir + f'fp_fn_corr_nmin{nmin}.tsv'):
        res = pd.read_csv(outdir + f'fp_fn_corr_nmin{nmin}.tsv', sep='\t')
    else:
        res = []
        for simu in simuls:
            for space in spaces:
                fpra_summary = pd.read_csv("../analyses/article_results_{}_{}/FPRA_summary.tsv".format(simu, space), sep="\t", index_col=[0,1])
                for tool in fnames:
                    df = fpra_summary.xs(tool, level=0)
                    fpra = pd.pivot_table(df, index="feature", columns="sample", values="FPRA", fill_value=0.)
                    fpra = fpra[(fpra > 0).sum(axis=1) > nmin]
                    fnra = pd.pivot_table(df, index="feature", columns="sample", values="FNRA", fill_value=0.)
                    fnra = fnra[(fnra > 0).sum(axis=1) >= nmin]
                    coincidence = (fpra > 0).astype(np.float64).dot((fnra > 0).astype(np.float64).T)
                    coinc = coincidence.stack()
                    coinc = coinc[coinc > nmin].sort_values(ascending=False)

                    corr = []
                    for (taxo_fp, taxo_fn), n in coinc.items():
                        y = fpra.loc[taxo_fp]
                        x = fnra.loc[taxo_fn]
                        idx = x[x > 0].index.intersection(y[y > 0].index)
                        corr.append((taxo_fp, taxo_fn, n, np.corrcoef(x[idx], y[idx])[0, 1])) 

                    corr = pd.DataFrame(corr, columns=['FP', 'FN', 'n', 'r'])#.set_index(['FP', 'FN'])
                    corr[['simu', 'space', 'tool']] = [simu, space, tool]
                    res.append(corr)

        res = pd.concat(res)
        res = res[res['r'] > 0.5]
        res.to_csv(outdir + f'fp_fn_corr_nmin{nmin}.tsv', sep='\t', index=False)


    #most frequently confused pairs
    rmin = 0.9
    res = res[res['r'] > rmin]

    df = res.groupby(by=['FP', 'FN']).size()
    df.to_csv(outdir + 'frequently_confused_pairs_nmin{}_rmin{:.2f}.tsv'.format(nmin, rmin), sep='\t')
    df = df[df > 2].rename('counts').sort_values(ascending=False).reset_index()
    df['FP'] = df['FP'].str.replace('unclassified', ';s__unclassified')
    df['label'] = df.apply(lambda x: x['FP'].split(';s__')[1] + ' vs. ' + x['FN'].split(';s__')[1], axis=1)

    nn = np.arange(df.shape[0])

    dg = []
    for idx, row in df.iterrows():
        dg.append(res[(res['FP'] == row['FP']) & (res['FN'] == row['FN'])])
    dg = pd.concat(dg)

    dg_wide = dg.pivot_table(index=['FP', 'FN'], columns=['simu', 'space', 'tool'], values='r')
    dg_wide = dg_wide.loc[df.set_index(['FP', 'FN']).index]
    dg_wide.to_csv(outdir + 'frequently_confused_to_plot.tsv', sep='\t')
