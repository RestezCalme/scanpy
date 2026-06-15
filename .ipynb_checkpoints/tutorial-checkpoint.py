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
help(adata)
