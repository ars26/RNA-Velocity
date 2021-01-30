#This file allows the usage of RNA Velocity to visualise the expression of the tumor micro-envoirnment in Hodgkin’s Lymphoma.
#RNA velocity is predicated upon a simple fact: unspliced and spliced Messenger RNA (mRNA) can be used to find the time derivative of the gene expression state.

#Import datasets and packages.

import numpy as np
import pandas as pd
import scanpy as sc
import scvelo as scv
sc.settings.verbosity = 3 
sc.logging.print_header()
sc.settings.set_figure_params(dpi=80, facecolor='white')
scv.logging.print_version()
scv.settings.verbosity = 3  # show errors(0), warnings(1), info(2), hints(3)
scv.settings.presenter_view = True  # set max width size for presenter view
scv.settings.set_figure_params('scvelo')  # for beautified visualization
results_file = 'results_HL.h5ad'
adata_r = sc.read_10x_mtx('/Users/bahawarsdhillon/Desktop/BIO 257 - Applied Genomics/Scanpy-Project/filtered_feature_bc_matrix/', var_names='gene_symbols',cache=False)
adata_r.var_names_make_unique()
adata_r

velocity = scv.read("/Users/bahawarsdhillon/Desktop/BIO 257 - Applied Genomics/Scanpy-Project/Parent_NGSC3_DI_HodgkinsLymphoma_possorted_genome_bam_JLA4X.loom")

adata = scv.utils.merge(adata_r, velocity)

adata_r = scv.utils.merge(adata_r, velocity)

sc.pp.filter_cells(adata_r, min_genes=200)
sc.pp.filter_genes(adata_r, min_cells=3)
sc.pp.filter_cells(adata_r, max_counts=39766)
sc.pp.filter_cells(adata_r, max_genes=5942)
adata_r.var['mt'] = adata_r.var_names.str.startswith('MT-')  # annotate the group of mitochondrial genes as 'mt'
sc.pp.calculate_qc_metrics(adata_r, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
adata = adata_r[adata_r.obs.pct_counts_mt < 10, :]

sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
adata.raw = adata
sc.pp.regress_out(adata, 'pct_counts_mt')
sc.pp.scale(adata, max_value=10)
cell_cycle_genes = [x.strip() for x in open('cell_cycle.txt')]

# Split into 2 lists
s_genes = cell_cycle_genes[:43]
g2m_genes = cell_cycle_genes[43:]

cell_cycle_genes = [x for x in cell_cycle_genes if x in adata.var_names]
sc.tl.score_genes_cell_cycle(adata, s_genes=s_genes, g2m_genes=g2m_genes)
sc.tl.pca(adata, svd_solver='arpack')
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
sc.tl.umap(adata)
sc.tl.leiden(adata, resolution=0.5)
sc.pl.umap(adata, color=['leiden'])
sc.tl.score_genes(adata, ["CD4"], score_name="CD4X")
cd4s = adata[adata.obs.CD4X > 0, :]
cd4s
sc.pl.umap(cd4s, color="leiden")
sc.pl.umap(adata, color="leiden")
sc.tl.score_genes(adata, ["CD8A"], score_name="CD8X")
tcells = adata[adata.obs.CD8X > 0, :]
tcells
sc.pl.umap(tcells, color="leiden")
sc.pl.umap(adata, color="leiden")
sc.pl.umap(cd4s, color="CD8A", color_map="magma")
sc.pl.umap(tcells, color="CD4", color_map="magma")
sc.pl.umap(tcells, color="CD4", color_map="magma")
sc.tl.score_genes(adata, ["CD3E"], score_name="tvel")
tvel = adata[adata.obs.tvel > 0, :]
tvel
sc.pl.umap(tvel, color="leiden")
sc.pl.umap(adata, color="leiden")
tvel
scv.pl.proportions(tvel, groupby="leiden")
scv.pp.moments(tvel, n_pcs=30, n_neighbors=30)
scv.tl.recover_dynamics(tvel)
scv.tl.velocity(tvel, mode='dynamical')
scv.tl.velocity_graph(tvel)
scv.pl.velocity_embedding_stream(tvel, basis='umap', color="leiden")
scv.pl.velocity_embedding(tvel, arrow_length=3, arrow_size=2, dpi=120, color="leiden")
scv.tl.velocity_confidence(tvel)
keys = 'velocity_length', 'velocity_confidence'
scv.pl.scatter(tvel, c=keys, cmap='coolwarm', perc=[5, 95])
scv.tl.latent_time(tvel)
top_genes = tvel.var['fit_likelihood'].sort_values(ascending=False).index
scv.pl.scatter(tvel, basis=top_genes[:15], ncols=5, frameon=False, color="leiden")
scv.pl.velocity(tvel, var_names= ["IL2RA"], color="leiden")
scv.logging.print_version()
scv.settings.verbosity = 3  # show errors(0), warnings(1), info(2), hints(3)
scv.settings.presenter_view = True  # set max width size for presenter view
scv.settings.set_figure_params('scvelo')  # for beautified visualization
scv.pl.velocity(tvel, var_names= ["IL2RA"], color="leiden")
scv.pl.velocity(tvel, var_names= ["IFNG", "TBX21"], color="leiden")
#naive cells
scv.pl.velocity(tvel, var_names= ["CCR7", "IL7R", "LEF1"], color="leiden")
#Treg markers
scv.pl.velocity(tvel, var_names= ["IL2RA", "FOXP3", "CTLA4", "LAG3", "TNFRSF18", "IKZF2", "IKZF4", "LGMN"], color="leiden")

sc.pl.umap(tvel, color="CD8B", color_map="magma")
#CD8 cytotoxic markers
scv.pl.velocity(tvel, var_names= ["GZMK", "GNLY", "GZMA"], color="leiden")
sc.tl.score_genes(cd4s, ["CD8A"], score_name="CD8X2")
dp = cd4s[cd4s.obs.CD8X2 > 0, :]
dp
sc.pl.umap(dp, color="leiden")
sc.pl.umap(cd4s, color="leiden")
#batch 0 is CD4 and CD8 double positive cells
sc.pl.violin(merge, keys=["IL21", "IL4", "IL2RA"], groupby="batch")
sc.tl.rank_genes_groups(merge, "batch", method="wilcoxon")

#variably expressed genes for the double positive batch are mostly associated with various cancers???
pd.DataFrame(merge.uns['rank_genes_groups']['names']).head(20)
sc.pl.umap(adata, color="leiden")
bcells = adata[adata.obs["leiden"].isin(["2"])]

sc.tl.score_genes(bcells, ["TNFRSF8", "CD274", "IRF4", "PAX5", "FUT4"], score_name="HRS")

sc.pl.umap(bcells, color="HRS", color_map="magma")
allt = tvel

sc.tl.pca(allt, svd_solver='arpack')
sc.pp.neighbors(allt, n_neighbors=10, n_pcs=40)
sc.tl.umap(allt)
sc.tl.leiden(allt, resolution=0.5)
sc.pl.umap(allt, color=['leiden'])
cd8_markers= ['GZMA', 'CCL5', 'GZMB', 'KLRD1', 'NKG7', 'CD8A', 'KLRK1', 'LGALS1', 'KLRC1', 'KLRG1', 'LGALS3', 'GZMK', 'CX3CR1', 'PLAC8', 'ZEB2', 'CTSD', 'H2AFZ', 'EPSTI1', 'DNAJC15', 'CD48', 'CTSW', 'SELPLG', 'TMSB4X', 'PLEK', 'TM6SF1', 'S100A10', 'VIM', 'ZYX', 'ABRACL', 'CRIP1', 'RACGAP1', 'CYBA', 'CCR2', 'PYCARD', 'SUB1', 'RAP1B', 'GLIPR2', 'EMP3', 'KRTCAP2', 'ID2', 'RPA2', 'STMN1', 'AHNAK', 'DAPK2', 'OSTF1', 'ITGAX', 'ANXA6', 'CD47', 'ARL4C', 'RUNX3', 'IFNGR1', 'LFNG', 'RASGRP2', 'REEP5', 'CHSY1', 'IL18RAP', 'JAK1', 'SGK1', 'KLRC2', 'TBX21', 'SEMA4A', 'SLAMF7', 'TXNDC5', 'ST3GAL6', 'CCNB2', 'CTSC', 'BIRC5', 'ITGB2', 'SP100', 'LAMB3', 'UBE2C', 'ANXA2', 'CENPW', 'CALM1', 'COX5B', 'COX5A', 'SNX10', 'EFHD2', 'LSP1', 'SERPINB9', 'HMGB2', 'THY1', 'RBM3', 'SNX5', 'PPP3CA', 'ARL6IP5', 'CKS1B', 'SEC61B', 'GNA15', 'NPM3', 'UBE2G2', 'ACTB', 'UGCG', 'S1PR4', 'CLIC1', 'LMNB1', 'CKS2', 'TPM4', 'TAGLN2', 'NDUFB7', 'PTP4A3', 'H2AFY', 'HMGN2', 'YWHAQ', 'GGH', 'DNAJC9', 'CMPK1', 'GIMAP7', 'CCND3', 'HCST', 'SMC2', 'ANAPC5', 'DEK', 'LSM5', 'CMTM7', 'PCNA', 'H2AFV', 'RAN', 'ANP32E', 'CENPA', 'SLBP', 'DUT', 'TMPO', 'H2AFX', 'TUBB4B', 'UBE2S', 'TUBA1B']

sc.tl.score_genes(allt, cd8_markers, score_name="cd8_ciucci")
sc.pl.violin(allt, keys="cd8_ciucci", groupby="leiden")
treg_markers = ['FOXP3', 'IKZF2', 'TNFRSF4', 'CAPG', 'TNFRSF18', 'CD74', 'CTLA4', 'RGS1', 'CCND2', 'SELL', 'SERINC3', 'SAMSN1', 'IFNGR1', 'GIMAP7', 'LTB', 'BTG1', 'IL7R', 'SDF4', 'CD2', 'SHISA5', 'GPX4', 'MBNL1', 'PELI1']

sc.tl.score_genes(allt, treg_markers, score_name="treg_ciucci")
sc.pl.violin(allt, keys="treg_ciucci", groupby="leiden")
tfh_markers = ['PDCD1', 'RGS10', 'CXCR5', 'TOX', 'TOX2', 'PTRH1', 'ANGPTL2', 'MARCKSL1', 'SMCO4', 'TNFSF8', 'PPP1R14B', 'TNFAIP8', 'HIF1A', 'LPP', 'MAF', 'CD160', 'SH2D1A', 'CD200', 'TBC1D4', 'TPI1', 'RPSA', 'HMGB1', 'ICOS', 'GAPDH', 'PRKCA', 'PTPRCAP', 'ZAP70', 'PTPN11', 'DENND2D', 'ZFP36L1', 'FYN', 'GDI2', 'LIMD2', 'PKM', 'NT5E', 'CTSB', 'PFKL', 'PGAM1', 'MATK', 'TRIM8', 'COX14', 'ASAP1', 'GNG2', 'P2RX7', 'LRMP', 'CD3G', 'ALDOA', 'GNA13', 'ISG15', 'RAB37', 'MMD', 'FAM162A', 'BORCS8', 'DDIT4', 'PTP4A2', 'CD82', 'MIF', 'PFKP', 'GIMAP5', 'EEA1', 'BATF']

sc.tl.score_genes(allt, tfh_markers, score_name="tfh_ciucci")
sc.pl.violin(allt, keys="tfh_ciucci", groupby="leiden")
list_TCells=["CD8", "CD3E", "CD4"]
list_NKCells=["NCAM1","GZMK", "GNLY", "GZMA"]
list_cytoTOXIC=["GZMK", "GNLY", "GZMA", "CD8A", "CD8B"]

list_Tregs=["IL2RA", "FOXP3", "CTLA4", "LAG3", "TNFRSF18", "IKZF2", "IKZF4", "LGMN"]
list_Th1=["IFNG", "TBX21", "CD4"]
list_Th2=["GATA3", "CRTH2", "IL4", "IL13", "CD4"]
list_Th17=["CD161", "CCR4", "CD4"]
list_TFH=["PDCD1", "CXCR5", "BCL6", "CD4"]

list_NaiveT=["CCR7", "IL7R", "LEF1"]
list_ActivatedT=["IL2RA"] #CD25


sc.tl.score_genes(allt, list_Tregs, score_name="treg_shah")
sc.pl.violin(allt, keys="treg_shah", groupby="leiden")

sc.tl.score_genes(allt, list_Th1, score_name="th1_shah")
sc.pl.violin(allt, keys="th1_shah", groupby="leiden")

sc.tl.score_genes(allt, list_Th2, score_name="th2_shah")
sc.pl.violin(allt, keys="th2_shah", groupby="leiden")

sc.tl.score_genes(allt, list_Th17, score_name="th17_shah")
sc.pl.violin(allt, keys="th17_shah", groupby="leiden")

sc.tl.score_genes(allt, list_TFH, score_name="tfh_shah")
sc.pl.violin(allt, keys="tfh_shah", groupby="leiden")

sc.tl.score_genes(allt, list_NaiveT, score_name="naive_shah")
sc.pl.violin(allt, keys="naive_shah", groupby="leiden")

sc.tl.score_genes(allt, list_ActivatedT, score_name="active_shah")
sc.pl.violin(allt, keys="active_shah", groupby="leiden")

list_cytoTOXIC=["GZMK", "GNLY", "GZMA", "CD8A", "CD8B"]
sc.tl.score_genes(allt, list_cytoTOXIC, score_name="cd8_shah")
sc.pl.violin(allt, keys="cd8_shah", groupby="leiden")

