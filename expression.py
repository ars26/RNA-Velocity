#This file allows the usage of ScanPy to visualise the expression of the tumor micro-envoirnment in Hodgkin’s Lymphoma.
#Please refer to attached HTML file for outputs.

import numpy as np
import pandas as pd
import scanpy as sc
sc.settings.verbosity = 3 
sc.logging.print_header()
sc.settings.set_figure_params(dpi=80, facecolor='white')

results_file = 'results_HL.h5ad'
adata_r = sc.read_10x_mtx('/Users/bahawarsdhillon/Desktop/BIO 257 - Applied Genomics/Scanpy-Project/filtered_feature_bc_matrix/', var_names='gene_symbols',cache=False)
adata_r.var_names_make_unique()
adata_r
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
sc.pl.highly_variable_genes(adata)
adata.raw = adata
sc.pp.regress_out(adata, 'pct_counts_mt')
sc.pp.scale(adata, max_value=10)
cell_cycle_genes = [x.strip() for x in open('cell_cycle.txt')]
print(len(cell_cycle_genes))

# Split into 2 lists
s_genes = cell_cycle_genes[:43]
g2m_genes = cell_cycle_genes[43:]

cell_cycle_genes = [x for x in cell_cycle_genes if x in adata.var_names]
print(len(cell_cycle_genes))
sc.tl.score_genes_cell_cycle(adata, s_genes=s_genes, g2m_genes=g2m_genes)
sc.pl.violin(adata, ['S_score', 'G2M_score'],
             jitter=0.4)
sc.tl.pca(adata, svd_solver='arpack')
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
sc.tl.umap(adata)
sc.tl.leiden(adata, resolution=0.5)
sc.pl.umap(adata, color=['leiden'])
tfh = []
th1 = []
treg = []
memory = []
cd4 = []
cd8= []

with open("/Users/bahawarsdhillon/Desktop/BIO 257 - Applied Genomics/Scanpy-Project/CD8_ARS.txt",'r') as myfile:
    for line in myfile:
        line=line.strip('\n')
        cd8.append(line)

with open("/Users/bahawarsdhillon/Desktop/BIO 257 - Applied Genomics/Scanpy-Project/CD4_Ars.txt",'r') as myfile:
    for line in myfile:
        line=line.strip('\n')
        cd4.append(line)

with open("/Users/bahawarsdhillon/Desktop/BIO 257 - Applied Genomics/Scanpy-Project/TH1_ARS.txt",'r') as myfile:
    for line in myfile:
        line=line.strip('\n')
        th1.append(line)

with open("/Users/bahawarsdhillon/Desktop/BIO 257 - Applied Genomics/Scanpy-Project/tfh_Ars.txt",'r') as myfile:
    for line in myfile:
        line=line.strip('\n')
        tfh.append(line)

with open("/Users/bahawarsdhillon/Desktop/BIO 257 - Applied Genomics/Scanpy-Project/Treg_ARS.txt",'r') as myfile:
    for line in myfile:
        line=line.strip('\n')
        treg.append(line)

with open("/Users/bahawarsdhillon/Desktop/BIO 257 - Applied Genomics/Scanpy-Project/Memory_ARS.txt",'r') as myfile:
    for line in myfile:
        line=line.strip('\n')
        memory.append(line)
print("TH1: \n", len(th1))
th1 = [x for x in th1 if x in adata.var_names]
print(len(th1))

print("TFH: \n", len(tfh))
tfh = [x for x in tfh if x in adata.var_names]
print(len(tfh))

print("CD4: \n", len(cd4))
cd4 = [x for x in cd4 if x in adata.var_names]
print(len(cd4))

print("CD8: \n", len(cd8))
cd8 = [x for x in cd8 if x in adata.var_names]
print(len(cd8))

print("TREG: \n", len(treg))
treg = [x for x in treg if x in adata.var_names]
print(len(treg))

print("Memory: \n", len(memory))
memory = [x for x in memory if x in adata.var_names]
print(len(memory))

allgenes = adata.var_names
iglist = []
for i in range(1, len(allgenes)):
    if allgenes[i].startswith("IG"):
        print(allgenes[i])
        iglist.append(allgenes[i])

sc.tl.score_genes(adata, iglist, score_name="IG")
list_BCells=['CD19', 'SYK', 'PIK3AP1', 'AKT1', 'AKT2', 'AKT3', 'CHUK', 'IKBKB', 'IKBKG']
sc.tl.score_genes(adata, iglist, score_name="bcell")
sc.pl.umap(adata, color=["IG", "bcell"])

sc.tl.score_genes(adata, th1, score_name="th1")
sc.tl.score_genes(adata, tfh, score_name="tfh")
sc.tl.score_genes(adata, treg, score_name="treg")
sc.tl.score_genes(adata, cd8, score_name="cd8")
sc.tl.score_genes(adata, cd4, score_name="cd4")
sc.tl.score_genes(adata, memory, score_name="memory")


sc.pl.umap(adata, color=["cd4", "cd8"])
sc.pl.umap(adata, color=["tfh", "th1", "treg", "memory"])
naive = ["CD8A", "CD4"]
sc.tl.score_genes(adata, naive, score_name="naive")
sc.pl.umap(adata, color="naive")
list_NKCells=["NCAM1","GZMK", "GNLY", "GZMA"]
sc.tl.score_genes(adata, list_NKCells, score_name="nk")

sc.pl.umap(adata, color=["leiden", "nk"], ncols=1)

sc.pl.umap(adata, color="CD8A", size=25, color_map="magma")
sc.pl.umap(adata, color="NCAM1", size=25, color_map="magma" )
sc.pl.umap(adata, color="CD68", size=25, color_map="magma")
list_NaiveT=["CCR7", "IL7R", "LEF1"]
sc.tl.score_genes(adata, list_NaiveT, score_name="naiveT")
sc.pl.umap(adata, color="naiveT", size=25, color_map="magma")
list_TCells=["CD8", "CD3E", "CD4"]
sc.tl.score_genes(adata, list_TCells, score_name="allT")
sc.pl.umap(adata, color="allT", size=25, color_map="magma")
list_BCells=['CD19', 'SYK', 'PIK3AP1', 'AKT1', 'AKT2', 'AKT3', 'CHUK', 'IKBKB', 'IKBKG', 'IGSF9', 'IGSF8', 'IGKC', 'IGKV3-15', 'IGHGP', 'IGFLR1', 'IGLV3-21', 'IGBP1-AS1']
print("IG: \n", len(iglist))
treg = [x for x in iglist if x in adata.uns.hvg]
print(len(iglist))
hvgs = adata.var_vector(k="highly_variable")
allgenes = adata.var_names
thelists = []
for d in range (0, len(allgenes)):
    if (hvgs[d]==True):
        thelists.append(allgenes[d])
        
print(len(thelists))


cov=0
iglistv2=[]
for d in range (0,len(iglist)):
    if iglist[d] in thelists:
        iglistv2.append(iglist[d])

print(iglistv2)
cd8v2 = []
for i in range (0, len(cd8)):
    if (cd8[i] not in cd4):
        cd8v2.append(cd8[i])
print(cd8v2)
print(len(cd8))
print(len(cd8v2))
cd4v2 = []
for i in range (0, len(cd4)):
    if (cd4[i] not in cd8):
        cd4v2.append(cd4[i])

print(len(cd4))
print(len(cd4v2))
list_BIG=['CD19', 'SYK', 'PIK3AP1', 'AKT1', 'AKT2', 'AKT3', 'CHUK', 'IKBKB', 'IKBKG', 'IGSF9', 'IGSF8', 'IGKC', 'IGKV3-15', 'IGHGP', 'IGFLR1', 'IGLV3-21', 'IGBP1-AS1']


sc.tl.leiden(adata, resolution=0.3)
sc.tl.rank_genes_groups(adata, "leiden", method="wilcoxon")
sc.tl.marker_gene_overlap(adata, allmarkers)
allmarkers = {"Tfh": tfh, "Th1": th1, "NK": list_NKCells, "B": list_BIG, "CD4": cd4v2, "CD8": cd8v2, "4": ["CD4"], "8": ["CD8A"]}
sc.tl.marker_gene_overlap(adata, allmarkers, top_n_markers=50, normalize="reference")
sc.pl.umap(adata, color="leiden")
sc.tl.score_genes(adata, list_BIG, score_name="Bv2")
sc.pl.umap(adata, color="Bv2")

pd.DataFrame(adata.uns['rank_genes_groups']['names']).head(10)
sc.pl.umap(adata, color="leiden", size=25)
sc.pl.umap(adata, color=["n macs", "allT"])
sc.pl.umap(adata, color="phase")
sc.pl.umap(adata, color="total_counts")
sc.pl.umap(adata, color=["CD4","CD163","n_genes"], color_map="magma")
list_HRS= ["TNFRSF8", "CD274", "IRF4", "PAX5", "FUT4"]
sc.tl.score_genes(adata, list_HRS, score_name="HRS")
sc.pl.umap(adata, color="HRS", color_map="magma")
sc.pl.umap(adata, color=["CD3E", "allT"], color_map="magma")
sc.pl.umap(adata, color=["IL2RA", "CD69"], color_map="magma")
sc.pl.umap(adata, color=["CCR7", "IL7R", "LEF1"], color_map="magma")
sc.pl.umap(adata, color=["IL2RA", "FOXP3", "CTLA4", "LAG3", "TNFRSF18", "IKZF2", "IKZF4", "LGMN"], color_map="magma")
list_Macrophages=["CD68", "ITGAM", "MPEG1", "CD163", "IRF7", "IFI44"]
macrophagedict = {"Macrophage": list_Macrophages}
sc.pl.matrixplot(adata, macrophagedict, groupby='leiden')
sc.pl.umap(adata, color=["IL9R", "CD33", "IL5RA"], color_map="magma")
sc.pl.violin(adata, keys=["IL9R", "CD33", "IL5RA"], groupby="leiden")
sc.pl.violin(adata, keys=["CD4", "CD8A"], groupby="leiden")
sc.pl.umap(adata, color=["leiden", "CD8A", "CD4"], color_map="magma")
sc.pl.umap(adata, color=["leiden", "IL2RA", "CD4", "CD8A"], color_map="magma")
sc.tl.score_genes(adata, ["IL2RA", "CD4", "CD8A"], score_name="myscore")
sc.pl.umap(adata, color="myscore", color_map="magma")
sc.pl.violin(adata, keys=["th1", "tfh", "treg", "memory", "cd4", "cd8"], groupby="leiden")
sc.pl.violin(adata, keys=["CD19", "CR2", "CD83", "NT5E"] , groupby="leiden")
sc.pl.umap(adata, color=["leiden", "CD19", "CR2", "CD83", "NT5E"], color_map="magma")
sc.pl.umap(adata, color="leiden", color_map="magma")
t = adata[adata.obs['leiden'].isin(["0", "1", "3", "5", "6", "8"])]
sc.pl.umap(t, color="leiden")
t2 = t
sc.pp.highly_variable_genes(t2, min_mean=0.0125, max_mean=3, min_disp=0.5)
sc.pl.highly_variable_genes(t2)
t2.raw = t2
sc.tl.pca(t, svd_solver='arpack')
sc.tl.leiden(t2, resolution=0.4)
sc.tl.rank_genes_groups(t2, 'leiden', method='wilcoxon')
sc.pl.rank_genes_groups(t2, n_genes=25, sharey=False)
pd.DataFrame(t2.uns['rank_genes_groups']['names']).head(10)
sc.pp.neighbors(t2, n_neighbors=10, n_pcs=40)
sc.tl.umap(t2)
sc.pl.umap(t2, color="CD8A", color_map="magma")
sc.pl.umap(t2, color="CD4", color_map="magma")
sc.pl.umap(t2, color="leiden", color_map="magma")
sc.tl.score_genes(t2, cd4v2, score_name="cd4T")
sc.tl.score_genes(t2, cd8v2, score_name="cd8T")
sc.tl.score_genes(t2, th1, score_name="th1T")
sc.tl.score_genes(t2, treg, score_name="tregT")
sc.tl.score_genes(t2, memory, score_name="memoryT")

sc.pl.umap(t2, color=["cd4T", "cd8T"], color_map="magma")
sc.pl.umap(t2, color=["th1T", "tregT", "memoryT"], color_map="magma")
list_Th1=["IFNG", "TBX21", "CD4"]
sc.tl.score_genes(t2, list_Th1, score_name="Th1v2")
sc.pl.umap(t2, color="Th1v2", color_map="magma")
sc.pl.umap(t2, color=["leiden", "CD4", "CD8A"], color_map="magma")
sc.pl.violin(t2, keys=["IFNG", "TBX21"], groupby="leiden")
sc.pl.violin(t2, keys=["S_score","G2M_score", "n_genes"], groupby="leiden")
t3 = t2[t2.obs['leiden'].isin(["0", "1", "2", "4", "5"])]

sc.pp.highly_variable_genes(t3, min_mean=0.0125, max_mean=3, min_disp=0.5)
t3.raw = t3
sc.tl.pca(t, svd_solver='arpack')

sc.tl.leiden(t3, resolution=0.3)
sc.tl.rank_genes_groups(t3, 'leiden', method='wilcoxon')
pd.DataFrame(t3.uns['rank_genes_groups']['names']).head(10)
sc.pl.violin(t3, keys=["CD4", "CD8A"], groupby="leiden")
sc.pl.violin(t3, keys=["GZMK", "GZMB", "GNLY"], groupby="leiden")
sc.pl.umap(t3, color="leiden")
sc.pl.violin(t3, keys=["IL2RA", "CD69"], groupby="leiden")
sc.pl.umap(t, color="leiden")
t4 = adata[adata.obs['leiden'].isin(["0", "1", "3", "5", "6", "8"])]
sc.pl.umap(t4, color=["leiden", "CD4", "CD8A"], color_map="magma")
sc.tl.score_genes(t4, ["CD8A"], score_name="CD8X")
t4
sc8 = t4.obs_vector(k='CD8X')
print(np.median(sc8))
print(np.std(sc8))
print(np.median(sc8)+ np.std(sc8))
t4.obs.columns
cd8 = t4[t4.obs.CD8X > 0.431045, :]
cd8
sc.pl.umap(t4, color="leiden")
sc.pl.umap(cd8, color="leiden")
sc.tl.leiden(cd8, resolution=0.6)
sc.pl.umap(cd8, color=['leiden'])
sc.tl.rank_genes_groups(cd8, 'leiden', method='wilcoxon')
pd.DataFrame(cd8.uns['rank_genes_groups']['names']).head(5)
sc.tl.score_genes(t4, ["CD4"], score_name="CD4X")
sc4 = t4.obs_vector(k='CD4X')
print(np.median(sc4))
print(np.std(sc4))
print(np.median(sc4)+ np.std(sc4))
cd4 = t4[t4.obs.CD4X > 0.17377234, :]
cd4
sc.pl.umap(t4, color="leiden")
sc.pl.umap(cd4, color="leiden")
t_cytotoxic=["GZMK", "GNLY", "GZMA"]
t_Tregs=["IL2RA", "FOXP3", "CTLA4", "LAG3", "TNFRSF18", "IKZF2", "IKZF4", "LGMN"]
t_Th1=["IFNG", "TBX21", "CD3E", "CD4"]
t_Th2=["GATA3", "CRTH2", "IL4", "IL13", "CD3E", "CD4"]
t_Th17=["CD161", "CCR4", "CD3E", "CD4"]
t_TFH=["PDCD1", "CXCR5", "BCL6", "CD3E", "CD4"]


sc.tl.score_genes(t4, t_Tregs, score_name="Tregs")
sc.tl.score_genes(t4, t_Th1, score_name="TH1")
sc.tl.score_genes(t4, t_Th17, score_name="TH17")
sc.tl.score_genes(t4, t_Th2, score_name="TH2")
sc.tl.score_genes(t4, t_TFH, score_name="TFH")
sc.tl.score_genes(t4, t_cytotoxic, score_name="cytotoxic")

sc.pl.violin(t4, keys=["Tregs","TH1", "TH2", "TH17", "cytotoxic", "TFH"], groupby="leiden")
sc.pl.umap(adata, color="leiden")
