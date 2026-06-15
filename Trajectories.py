from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
import scanpy as sc
sc.settings.verbosity = 3  # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.logging.print_versions()
results_file = "./write/paul15.h5ad"
# low dpi (dots per inch) yields small inline figures
sc.settings.set_figure_params(dpi=80, frameon=False, figsize=(3, 3), facecolor="white")
adata = sc.datasets.paul15()
adata
# this is not required and results will be comparable without it
adata.X = adata.X.astype("float64")
sc.pp.recipe_zheng17(adata)
sc.tl.pca(adata, svd_solver="arpack")
sc.pp.neighbors(adata, n_neighbors=4, n_pcs=20)
sc.tl.draw_graph(adata)
sc.pl.draw_graph(adata, color="paul15_clusters", legend_loc="on data")
sc.tl.diffmap(adata)
sc.pp.neighbors(adata, n_neighbors=10, use_rep="X_diffmap")
sc.tl.draw_graph(adata)
sc.pl.draw_graph(adata, color="paul15_clusters", legend_loc="on data")