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
space = "gtdb207"
names = ["kraken2 + bracken", "metaphlan3", "motus3", "metaphlan4", "biomscope"]
fnames = ["kraken2", "metaphlan3", "motus3", "metaphlan4", "biomscope"]


indir = "../data/{}_to_{}/".format(simu, space)
supersat = pd.read_csv(indir + "supersat_{}_{}.sat".format(simu, space), sep="\t", index_col=(0,1))
matched_sats = [supersat[name].unstack() for name in fnames]
sat_ref = supersat["reference"].unstack()
smpls = sat_ref.columns#.drop(["Species_rep", "mgnify_taxo"]).tolist()
nsmpls = len(smpls)

# Full gtdb207 taxonomy
refdir = "../data/taxo_indexes/"

if __name__ == "__main__":
    """
    Calculating UniFrac distance matrices in GTDB207 space
    """
    outdir = "../analyses/article_results_{}_gtdb207/unifrac/".format(simu)
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    # if True:
    print("Preparing phylogenetic tree")
    t0 = time()
    gtdb207_bac = pd.read_csv(refdir + "bac120_taxonomy_r207.tsv", sep="\t", header=None, names=["gtdb_accession", "gtdb_taxonomy"], index_col="gtdb_accession")
    # treepath_bac = "/home/puller/Documents/GMT_science/Kraken2/gtdb/bac120_r207.tree"
    tree_bac = Phylo.read(refdir + "bac120_r207.tree", format="newick")
    acc_bac = pd.Index([t.name for t in tree_bac.get_terminals()])

    df_bac = gtdb207_bac.loc[acc_bac]
    df_bac = df_bac[df_bac["gtdb_taxonomy"].isin(sat_ref.index)]
    df_bac.index.name="gtdb_accession"
    taxo_bac = df_bac["gtdb_taxonomy"].tolist()

    for genome, taxo in gtdb207_bac.loc[acc_bac].squeeze().items():
        if taxo not in taxo_bac:
            tree_bac.prune(genome)

    Phylo.write(tree_bac, outdir + "tree_bac.nwk", "newick")
    df_bac.to_csv(outdir + "tree_id_bac.tsv", sep="\t")
    t1 = time()
    print(t1-t0)    
    # else:
    #     df_bac = pd.read_csv(outdir + "tree_id_bac.tsv", sep="\t", index_col="gtdb_accession")

    #--------------------------------------------------------------------------
    #Creating unifrac matrices

    #load tree in scikit-bio format
    tree_bac = read(outdir + "tree_bac.nwk", format="newick", into=TreeNode, convert_underscores=False)
    ids = df_bac.index.tolist()
    sat_bac_ref = sat_ref.loc[df_bac["gtdb_taxonomy"]]
    X = sat_ref.loc[df_bac["gtdb_taxonomy"]].to_numpy()
    X = np.round(X / X[X>0].min())

    #unifrac matrices for reference data
    # if True:
    print("reference")
    t0 = time()    
    D_ref_uu, D_ref_wu = build_unifrac_matrices(sat_bac_ref.iloc[:,:nsmpls], ids, tree_bac)
    D_ref_uu.to_csv(outdir + "Duu_ref.tsv", sep="\t")
    D_ref_wu.to_csv(outdir + "Dwu_ref.tsv", sep="\t")
    t1 = time()
    print(t1-t0)    
    # else:
    #     D_ref_uu = pd.read_csv(outdir + "Duu_ref.tsv", sep="\t", index_col=0)
    #     D_ref_wu = pd.read_csv(outdir + "Dwu_ref.tsv", sep="\t", index_col=0)


    #calculate unifrac distances for pipelines' outputs
    sats_bac = [sat.loc[df_bac["gtdb_taxonomy"]] for sat in matched_sats]
    DD_uu = []
    DD_wu = []
    # if True:
    for fname, sat in zip(fnames, sats_bac):
        print(fname)
        t0 = time()
        D_uu, D_wu = build_unifrac_matrices(sat, ids, tree_bac)
        D_uu.to_csv(outdir + "Duu_{}.tsv".format(fname), sep="\t")
        D_wu.to_csv(outdir + "Dwu_{}.tsv".format(fname), sep="\t")
        DD_uu.append(D_uu)
        DD_wu.append(D_wu)
        t1 = time()
        print(t1 - t0)
    # else:
    #     for fname, sat in zip(fnames, matched_sats):
    #         D_uu = pd.read_csv(outdir + "Duu_{}.tsv".format(fname), sep="\t", index_col=0)
    #         D_wu = pd.read_csv(outdir + "Dwu_{}.tsv".format(fname), sep="\t", index_col=0)
    #         DD_uu.append(D_uu)
    #         DD_wu.append(D_wu)
    

    # #PCoA on unifrac distance matrices
    # nPC = 3
    # _, V_ref_uu = pcoa_dist(D_ref_uu, k=nPC)
    # _, V_ref_wu = pcoa_dist(D_ref_wu, k=nPC)

    # VV_uu = [pcoa_dist(D_uu, k=nPC)[1] for D_uu in DD_uu]
    # VV_wu = [pcoa_dist(D_wu, k=nPC)[1] for D_wu in DD_wu]

    # j = names.index("kraken2 + bracken")
    # # VV_uu[j][:,0] = -VV_uu[j][:,0]
    # VV_wu[j][:,1] = -VV_wu[j][:,1]

    # j = names.index("metaphlan3")
    # VV_uu[j][:,1] = -VV_uu[j][:,1]
    # # VV_wu[j][:,0] = -VV_wu[j][:,0]
    # # VV_wu[j][:,1] = -VV_wu[j][:,1]

    # j = names.index("motus3")
    # # VV_wu[j][:,0] = -VV_wu[j][:,0]

    # # j = names.index("biomscope")
    # # VV_uu[j][:,1] = -VV_uu[j][:,1]
    # # VV_wu[j][:,1] = -VV_wu[j][:,1]

    # with mpl.rc_context(ctxt1):
    #     fig, axs = plt.subplots(3, nsats, sharey="row", sharex="row")
    #     for ax, name, V_uu in zip(axs[0], names, VV_uu):
    #         p1, p2 = pairdist_components_plot0(V_ref_uu, V_uu, ax)
    #         ax.set_title(name)
    #         ax.set_xlabel('component 1')
    #     axs[0, 0].legend([p1, p2], ["Reference", "Result"], loc="lower left")#, bbox_to_anchor=(1.5, 0.5))
    #     axs[0, 0].set_ylabel('component 2')

    #     for ax, name, D_uu in zip(axs[1], names, DD_uu):
    #         pairdist_scatterplot0(D_ref_uu, D_uu, ax)
    #         ax.set_xlabel("distance ({})".format("Native space"))
    #     axs[1,0].set_ylabel("distance ({})".format("MSP space"))

    #     for ax, name, D_uu in zip(axs[2], names, DD_uu):
    #         pairdist_blandaltman0(D_ref_uu, D_uu, ax)
    #         ax.set_xlabel("distance (mean)")
    #     axs[2,0].set_ylabel("distance (difference)")

    #     plt.savefig(figdir + "figure3_uu.png")
    #     plt.close(fig)

    #     fig, axs = plt.subplots(3, nsats, sharey="row", sharex="row")
    #     for ax, name, V_wu in zip(axs[0], names, VV_wu):
    #         p1, p2 = pairdist_components_plot0(V_ref_wu, V_wu, ax)
    #         ax.set_title(name)
    #         ax.set_xlabel('component 1')
    #     axs[0, 0].legend([p1, p2], ["Reference", "Result"], loc="lower left")#, bbox_to_anchor=(1.5, 0.5))
    #     axs[0, 0].set_ylabel('component 2')

    #     for ax, name, D_wu in zip(axs[1], names, DD_wu):
    #         pairdist_scatterplot0(D_ref_wu, D_wu, ax)
    #         ax.set_xlabel("distance ({})".format("Native space"))
    #     axs[1,0].set_ylabel("distance ({})".format("MSP space"))

    #     for ax, name, D_wu in zip(axs[2], names, DD_wu):
    #         pairdist_blandaltman0(D_ref_wu, D_wu, ax)
    #         ax.set_xlabel("distance (mean)")
    #     axs[2,0].set_ylabel("distance (difference)")

    #     plt.savefig(figdir + "figure3_wu.png")
    #     plt.close(fig)

    #--------------------------------------------------------------------
    #calculating intermatrix norms
    norm_ref_matched_uu = intermatrix_norms([D_ref_uu]*len(names), DD_uu, names, filepath=outdir +"norm_ref_matched_uu.tsv")
    norm_ref_matched_wu = intermatrix_norms([D_ref_wu]*len(names), DD_wu, names, filepath=outdir +"norm_ref_matched_wu.tsv")

    #--------------------------------------------------------------------
    #UniFrac distances between real sample and simulation
    # dist_to_ref = [{name: unifrac_dist_to_ref(sat_unifrac_ref, sat, ids, tree)} for name, sat in zip(names, sats_unifrac)]
    dist_to_ref = [unifrac_dist_to_ref(sat_bac_ref, sat, ids, tree_bac) for sat in sats_bac]
    dist_to_ref_uu = pd.concat([df["unweighted_unifrac"].rename(name) for name, df in zip(names, dist_to_ref)], axis=1)
    dist_to_ref_wu = pd.concat([df["weighted_unifrac"].rename(name) for name, df in zip(names, dist_to_ref)], axis=1)

    dist_to_ref_uu.to_csv(outdir + "dist_to_ref_uu.tsv", sep="\t")
    dist_to_ref_wu.to_csv(outdir + "dist_to_ref_wu.tsv", sep="\t")

    sys.exit()
    # fig, ax = plt.subplots(1, 1, figsize=(6, 6))
    # ax.boxplot(dist_to_ref_uu.T)
    # ax.set_ylabel("unweighted unifrac")
    # ax.set_xticks(np.arange(1, nsats+1))
    # ax.set_xticklabels(names, rotation=90)
    # plt.tight_layout()
    # plt.savefig(figdir + "dist_to_ref_box_uu.png")
    # plt.close(fig)

    # fig, ax = plt.subplots(1, 1, figsize=(6, 6))
    # ax.boxplot(dist_to_ref_wu.T)
    # ax.set_ylabel("weighted unifrac")
    # ax.set_xticks(np.arange(1, nsats+1))
    # ax.set_xticklabels(names, rotation=90)
    # plt.tight_layout()
    # plt.savefig(figdir + "dist_to_ref_box_wu.png")
    # plt.close(fig)

    # dist_to_ref_uu.sort_values(by="biomscope", inplace=True)
    # nn_smpls = np.arange(len(smpls))
    # fig, axs = plt.subplots(nsats, 1, figsize=(30, 3*nsats), sharey=True, sharex=True)
    # ax_iter = iter(axs)
    # for name, col in dist_to_ref_uu.items():
    #     ax = next(ax_iter)
    #     ax.bar(nn_smpls, col)
    #     ax.set_xticks(nn_smpls)
    #     ax.set_xticklabels(col.index, rotation=90, fontsize=10)
    #     ax.set_title(name)
    #     ax.set_ylabel("unweighted unifrac")
    # axs[-1].set_xlabel("samples")
    # plt.tight_layout()
    # plt.savefig(figdir + "dist_to_ref_uu.png")
    # plt.close(fig)

    # dist_to_ref_wu.sort_values(by="biomscope", inplace=True)
    # fig, axs = plt.subplots(nsats, 1, figsize=(30, 3*nsats), sharey=True, sharex=True)
    # ax_iter = iter(axs)
    # for name, col in dist_to_ref_wu.items():
    #     ax = next(ax_iter)
    #     ax.bar(nn_smpls, col)
    #     ax.set_xticks(nn_smpls)
    #     ax.set_xticklabels(col.index, rotation=90, fontsize=10)
    #     ax.set_title(name)
    #     ax.set_ylabel("weighted unifrac")
    # axs[-1].set_xlabel("samples")
    # plt.tight_layout()
    # plt.savefig(figdir + "dist_to_ref_wu.png")
    # plt.close(fig)
    # #--------------------------------------------------------------------
    # #non-metric multidimensional scaling (NMDS) using alternative metrics
    # metric = "braycurtis"#"cosine_sqrt"#"jensenshannon"#"braycurtis"

    # if metric == "cosine_sqrt":
    #     D_ref = pd.DataFrame(squareform(pdist(np.sqrt(sat_ref[smpls]).T, "cosine")), index=smpls, columns=smpls)
    # else:
    #     D_ref = pd.DataFrame(squareform(pdist(sat_ref[smpls].T, metric)), index=smpls, columns=smpls)
    # _, V_ref = pcoa_dist(D_ref, k=3)

    # DD_bac, VV_bac = pairwise_nmds(sats_bac, smpls, metric)

    # if metric == "braycurtis":
    #     j = names.index("kraken2 + bracken")
    #     VV_bac[j][:,1] = -VV_bac[j][:,1]

    #     j = names.index("metaphlan3")
    #     VV_bac[j][:,1] = -VV_bac[j][:,1]
    #     VV_bac[j][:,0] = -VV_bac[j][:,0]

    #     # j = names.index("biomscope")
    #     # VV_bac[j][:,1] = -VV_bac[j][:,1]
        
    # with mpl.rc_context(ctxt1):
    #     fig, axs = plt.subplots(3, nsats, sharey="row", sharex="row")
    #     for ax, name, V_bac in zip(axs[0], names, VV_bac):
    #         p1, p2 = pairdist_components_plot0(V_ref, V_bac, ax)
    #         ax.set_title(name)
    #         ax.set_xlabel('component 1')
    #     axs[0, 0].legend([p1, p2], ["Native space", "MSP space"], loc="lower left")#, bbox_to_anchor=(1.5, 0.5))
    #     axs[0, 0].set_ylabel('component 2')

    #     for ax, name, D_bac in zip(axs[1], names, DD_bac):
    #         pairdist_scatterplot0(D_ref, D_bac, ax)
    #         ax.set_xlabel("distance ({})".format("Native space"))
    #     axs[1,0].set_ylabel("distance ({})".format("MSP space"))

    #     for ax, name, D_bac in zip(axs[2], names, DD_bac):
    #         pairdist_blandaltman0(D_ref, D_bac, ax)
    #         ax.set_xlabel("distance (mean)")
    #     axs[2,0].set_ylabel("distance (difference)")

    #     plt.savefig(figdir + "figure3.png")
    #     plt.close(fig)

    # norm_ref_matched = intermatrix_norms([D_ref]*len(names), DD_bac, names, filepath=outdir +"norm_ref_matched_bac.tsv")
    # #--------------------------------------------------------------------
    # #distance between real sample and simulation
    # if metric == "braycurtis":
    #     dist_to_ref = pd.DataFrame(pd.DataFrame({name : [braycurtis(sat[smpl], sat_bac_ref[smpl]) for smpl in smpls] for name, sat in zip(names, sats_bac)}, index=smpls))
    #     lab="Bray-Curtis distance"
    # elif metric == "jensenshannon":
    #     dist_to_ref = pd.DataFrame(pd.DataFrame({name : [jensenshannon(sat[smpl], sat_bac_ref[smpl]) for smpl in smpls] for name, sat in zip(names, sats_bac)}, index=smpls))
    #     lab="Jensen-Shannon distance"
    # elif metric == "cosine_sqrt":
    #     dist_to_ref = pd.DataFrame(pd.DataFrame({name : [cosine_sqrt(sat[smpl], sat_bac_ref[smpl]) for smpl in smpls] for name, sat in zip(names, sats_bac)}, index=smpls))
    #     lab="cosine-sqrt distance"

    # dist_to_ref.to_csv(outdir + "dist_to_ref_bac.tsv", sep="\t")

    # dist_to_ref.sort_values(by="biomscope", inplace=True)
    # nn_smpls = np.arange(len(smpls))
    # fig, axs = plt.subplots(nsats, 1, figsize=(30, 3*nsats), sharey=True, sharex=True)
    # ax_iter = iter(axs)
    # for name, col in dist_to_ref.items():
    #     ax = next(ax_iter)
    #     ax.bar(nn_smpls, col)
    #     ax.set_xticks(nn_smpls)
    #     ax.set_xticklabels(col.index, rotation=90, fontsize=10)
    #     ax.set_title(name)
    #     ax.set_ylabel(lab)
    # axs[-1].set_xlabel("samples")
    # plt.tight_layout()
    # plt.savefig(figdir + "dist_to_ref_bac.png")
    # plt.close(fig)

    # fig, ax = plt.subplots(1, 1, figsize=(6, 6))
    # ax.boxplot(dist_to_ref.T)
    # ax.set_ylabel(lab)
    # # ax.set_xlabel("#samples")
    # ax.set_xticks(np.arange(1, dist_to_ref.columns.size+1))
    # ax.set_xticklabels(dist_to_ref.columns, rotation=90)
    #     # ax.yaxis.set_tick_params(which='both', labelleft=True)
    #     # ax.xaxis.set_tick_params(which='both', labelbottom=True)
    # plt.tight_layout()
    # plt.savefig(figdir + "dist_to_ref_box_bac.png")
    # plt.close(fig)

