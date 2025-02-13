# 🧬 clustermole_py: Cluster Annotation for Single-Cell Data in Python

[![License](https://img.shields.io/badge/License-MIT-blue.svg)](https://opensource.org/licenses/MIT)
[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/nmlakra/clustermole-py/HEAD?labpath=examples%2Fclustermole_py_usage.ipynb)

`clustermole_py` is a Python package inspired by the original [clustermole](https://github.com/igordot/clustermole) R package. It's designed to help you **annotate cell clusters from single-cell RNA-seq data** 📊 using powerful gene set enrichment analysis.

✨ Give your clusters biological meaning with **Enrichr** and **GSVA**! ✨


## 🚀 Key Features

* **Enrichr Integration :**
    * Direct query of the [Enrichr API](https://maayanlab.cloud/Enrichr/) for gene set enrichment analysis.
    * Specialized `Enrichr` module with methods like `get_enrichment()` and `get_cell_type_enrichment()`.
    *  ⚡️ **Multi-threaded `get_cell_type_enrichment()`** for fast cell type annotation using curated libraries.
* **GSVA Module :**
    *  Perform **Gene Set Variation Analysis (GSVA)** directly on your cluster data.
    *  `aggregate_expression()` function for easy normalization of cluster-level expression data from Scanpy `AnnData` objects.
    *  `gsva.get_enrichment()` for DE-free cluster annotation.
* **Scanpy Integration :**
    *  Designed to work seamlessly with [Scanpy](https://scanpy.readthedocs.io/en/stable/) `AnnData` objects.
    *  Example workflow uses Scanpy for data loading, clustering, and marker gene identification.
* **Clear and Interpretable Results :**
    *  Returns enrichment results in easy-to-read Pandas DataFrame.
    *  Output tables include p-values, adjusted p-values, combined scores, and overlapping genes for enriched terms.

## 📦 Installation

For this version, you can install `clustermole_py` directly from GitHub using pip:

```bash
pip install git+https://github.com/nmlakra/clustermole-py.git
```

Once installed, you can import and use `clustermole_py` in your Python environment.


## 🕹️ Quick Usage Example

Here's a simplified example of how to use `clustermole_py` to annotate cell clusters. For a more detailed walkthrough, check out the Jupyter Notebook in the `examples` directory\!

**(1) Prepare your single-cell data and perform clustering using Scanpy:**

```python
import scanpy as sc

# Load data and perform Leiden clustering (example using pbmc3k dataset)
adata = sc.datasets.pbmc3k_processed()
sc.tl.leiden(adata, flavor='igraph', n_iterations=2, resolution=0.5)

# Identify marker genes for each cluster (example using Wilcoxon rank-sum test)
sc.tl.rank_genes_groups(adata, 'leiden', method='wilcoxon')
top_n_markers = 25
cluster_marker_genes = {}
for cluster_id in adata.obs['leiden'].cat.categories:
    markers = sc.get.rank_genes_groups_df(adata, group=cluster_id, key='rank_genes_groups', return_names=True)
    top_genes = markers['names'][:top_n_markers].tolist()
    cluster_marker_genes[f"Leiden_Cluster_{cluster_id}"] = top_genes
```

**(2) Annotate clusters using the `Enrichr` module:**

```python
from clustermole_py.enrichr import Enrichr

# Get marker genes for a cluster (e.g., "Leiden_Cluster_1")
b_cell_markers = cluster_marker_genes["Leiden_Cluster_1"]

# Initialize Enrichr and get cell type enrichment
enrichr_cell_type = Enrichr(gene_list=b_cell_markers, pval_cutoff=0.05)
cell_type_results = enrichr_cell_type.get_cell_type_enrichment()

print(cell_type_results.head()) # Display top results
```

**(Example Output - first few rows of the table)**

```
|    | term name                             |     p-value |   odds ratio |   combined score | overlapping genes                                   |   adjusted p-value | gene_set                    |
|---:|:--------------------------------------|------------:|-------------:|-----------------:|:----------------------------------------------------|-------------------:|:----------------------------|
|  0 | B cell:Kidney                         | 5.08603e-18 |     103.437  |         4118.88  | ['SMIM14', 'EAF2', 'CD79B', 'CD79A', ...]           | 8.54452e-16        | CellMarker_Augmented_2021   |
|  1 | B Cell:Kidney Human                   | 4.98039e-18 |     103.601  |         4127.57  | ['SMIM14', 'EAF2', 'CD79B', 'CD79A', ...]           | 9.21372e-16        | CellMarker_2024             |
|  2 | B Cell:Lung Human                     | 1.92402e-16 |    1664.677  |        60239.22  | ['CD79B', 'CD79A', 'TCL1A', 'MZB1', ...]            | 1.77972e-14        | CellMarker_2024             |
... (and so on)
```

**(3) Annotate clusters using the `gsva` module (DE-free approach):**

```python
from clustermole_py import gsva

# Aggregate expression data
normalized_data_leiden = gsva.aggregate_expression(
    adata, groupby="leiden", scale_factor=1e4, layer="counts", apply_log1p=True
)

# Run GSVA enrichment
gsva_results_panglaodb = gsva.get_enrichment(
    normalized_data_leiden, gene_sets=['PanglaoDB_Augmented_2021']
)

print(gsva_results_panglaodb.head()) # Display top results
```

**(Example Output - first few rows of the table)**

```
|    | group        | term                     |   enrichment_score | gene_set                 |
|---:|:-------------|:-------------------------|-------------------:|:-------------------------|
|  0 | B cells      | B cell                   |          0.198753  | PanglaoDB_Augmented_2021 |
|  1 | NK cells     | NK_cell                  |          0.174589  | PanglaoDB_Augmented_2021 |
|  2 | T cells      | T_cell                   |          0.167892  | PanglaoDB_Augmented_2021 |
... (and so on)
```

## 📚 Documentation

Check out the [Example Notebook](examples/clustermole_py_usage.ipynb) for more information!
