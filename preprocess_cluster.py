# Core scverse libraries
from __future__ import annotations

import anndata as ad

# Data retrieval
import pooch
import scanpy as sc

EXAMPLE_DATA = pooch.create( 
    path=pooch.os_cache("scverse_tutorials"), ##获取操作系统标准的缓存目录路径，并在其中创建一个名为 "scverse_tutorials" 的子目录
    base_url="doi:10.6084/m9.figshare.22716739.v1/",##指定数据集的 DOI 地址
) ## 创建一个新的 Pooch 数据管理对象

EXAMPLE_DATA.load_registry_from_doi()##从指定的 DOI 加载数据注册表，注册表包含了数据文件的名称、URL 和校验信息等元数据

samples = {
    "s1d1": "s1d1_filtered_feature_bc_matrix.h5",
    "s1d3": "s1d3_filtered_feature_bc_matrix.h5",
}
adatas = {}

for sample_id, filename in samples.items():
    path = EXAMPLE_DATA.fetch(filename) ##使用 Pooch 对象的 fetch 方法下载指定的文件，并返回文件的本地路径
    sample_adata = sc.read_10x_h5(path) ##使用 Scanpy 库的 read_10x_h5 函数读取下载的 HDF5 文件，并将其转换为 AnnData 对象
    sample_adata.var_names_make_unique() ##自动为所有重复的基因名添加后缀，以确保每个基因名在 AnnData 对象中是唯一的
    adatas[sample_id] = sample_adata ##将处理后的 AnnData 对象存储在一个字典中，键为样本 ID，值为对应的 AnnData 对象

adata = ad.concat(adatas, label="sample") ##使用 AnnData 库的 concat 函数将多个 AnnData 对象合并成一个新的 AnnData 对象，并添加一个新的 obs 列 "sample" 来标识每个细胞所属的样本
adata.obs_names_make_unique() ##自动为所有重复的细胞 ID 添加后缀，以确保每个细胞 ID 在 AnnData 对象中是唯一的
print(adata.obs["sample"].value_counts()) 
adata
adata.var["mt"] = adata.var_names.str.startswith("MT-") ##在 AnnData 对象的 var 数据框中创建一个新的布尔列 "mt"，如果基因名以 "MT-" 开头，则该列的值为 True，否则为 False。
adata.var["ribo"] = adata.var_names.str.startswith(("RPS", "RPL"))
adata.var["hb"] = adata.var_names.str.contains("^HB[^(P)]")
sc.pp.calculate_qc_metrics(adata, qc_vars=["mt", "ribo", "hb"], inplace=True, log1p=True) 
sc.pl.violin(
    adata,
    ["n_genes_by_counts", "total_counts", "pct_counts_mt"],
    jitter=0.4,
    multi_panel=True,
) ##使用 Scanpy 库的 violin 函数绘制小提琴图，展示每个细胞的基因数量、总计数和线粒体基因占比等质量控制指标的分布情况
sc.pl.scatter(adata, "total_counts", "n_genes_by_counts", color="pct_counts_mt") 
sc.pp.filter_cells(adata, min_genes=100)
sc.pp.filter_genes(adata, min_cells=3) 
sc.pp.scrublet(adata, batch_key="sample") 
##归一化
adata.layers["counts"] = adata.X.copy()
sc.pp.normalize_total(adata)
sc.pp.log1p(adata)
##特征选择
sc.pp.highly_variable_genes(adata, n_top_genes=2000, batch_key="sample")
sc.pl.highly_variable_genes(adata)
##降维
sc.tl.pca(adata) ##默认使用高变基因
sc.pl.pca_variance_ratio(adata, n_pcs=50, log=True)
sc.pl.pca(
    adata,
    color=["sample", "sample", "pct_counts_mt", "pct_counts_mt"],
    dimensions=[(0, 1), (2, 3), (0, 1), (2, 3)],
    ncols=2,
    size=2,
) ##主成分图
##最近邻图的构造与可视化
sc.pp.neighbors(adata)## 默认使用 PCA 结果
sc.tl.umap(adata)
sc.pl.umap(
    adata,
    color="sample",
    # Setting a smaller point size to get prevent overlap
    size=2,
)
##聚类
sc.tl.leiden(adata, flavor="igraph", n_iterations=2)
sc.pl.umap(adata, color=["leiden"])
##重新评估质量控制和细胞过滤
sc.pl.umap(
    adata,
    color=["leiden", "predicted_doublet", "doublet_score"],
    # increase horizontal space between panels
    wspace=0.5,
    size=3,
)
## 手动单元类型注释
for res in [0.02, 0.5, 2.0]:
    sc.tl.leiden(adata, key_added=f"leiden_res_{res:4.2f}", resolution=res, flavor="igraph")
sc.pl.umap(
    adata,
    color=["leiden_res_0.02", "leiden_res_0.50", "leiden_res_2.00"],
    legend_loc="on data",
)
#标记基因集
marker_genes = {
    "CD14+ Mono": ["FCN1", "CD14"],
    "CD16+ Mono": ["TCF7L2", "FCGR3A", "LYN"],
    # Note: DMXL2 should be negative
    "cDC2": ["CST3", "COTL1", "LYZ", "DMXL2", "CLEC10A", "FCER1A"],
    "Erythroblast": ["MKI67", "HBA1", "HBB"],
    # Note HBM and GYPA are negative markers
    "Proerythroblast": ["CDK6", "SYNGR1", "HBM", "GYPA"],
    "NK": ["GNLY", "NKG7", "CD247", "FCER1G", "TYROBP", "KLRG1", "FCGR3A"],
    "ILC": ["ID2", "PLCG2", "GNLY", "SYNE1"],
    "Naive CD20+ B": ["MS4A1", "IL4R", "IGHD", "FCRL1", "IGHM"],
    # Note IGHD and IGHM are negative markers
    "B cells": [
        "MS4A1",
        "ITGB1",
        "COL4A4",
        "PRDM1",
        "IRF4",
        "PAX5",
        "BCL11A",
        "BLK",
        "IGHD",
        "IGHM",
    ],
    "Plasma cells": ["MZB1", "HSP90B1", "FNDC3B", "PRDM1", "IGKC", "JCHAIN"],
    # Note PAX5 is a negative marker
    "Plasmablast": ["XBP1", "PRDM1", "PAX5"],
    "CD4+ T": ["CD4", "IL7R", "TRBC2"],
    "CD8+ T": ["CD8A", "CD8B", "GZMK", "GZMA", "CCL5", "GZMB", "GZMH", "GZMA"],
    "T naive": ["LEF1", "CCR7", "TCF7"],
    "pDC": ["GZMB", "IL3RA", "COBLL1", "TCF4"],
}
sc.pl.dotplot(adata, marker_genes, groupby="leiden_res_0.02", standard_scale="var")
adata.obs["cell_type_lvl1"] = adata.obs["leiden_res_0.02"].map(
    {
        "0": "Lymphocytes",
        "1": "Monocytes",
        "2": "Erythroid",
        "3": "B Cells",
    }
)
sc.pl.dotplot(adata, marker_genes, groupby="leiden_res_0.50", standard_scale="var")
# Obtain cluster-specific differentially expressed genes
sc.tl.rank_genes_groups(adata, groupby="leiden_res_0.50", method="wilcoxon")
sc.pl.rank_genes_groups_dotplot(adata, groupby="leiden_res_0.50", standard_scale="var", n_genes=5)