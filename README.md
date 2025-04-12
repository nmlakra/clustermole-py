# 🧬 clustermolepy: Fast Enrichment for Annotating Single-Cell Clusters


[![License](https://img.shields.io/badge/License-MIT-blue.svg)](https://opensource.org/licenses/MIT)
[![Python 3.10](https://img.shields.io/badge/python-3.10-blue.svg)](https://www.python.org/downloads/release/python-3100/)
[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/nmlakra/clustermole-py/HEAD?labpath=examples%2Fclustermolepy_usage.ipynb)

`clustermolepy` is a light weight Python package inspired by the original [clustermole](https://github.com/igordot/clustermole) R package. It's designed to help you **annotate cell clusters from single-cell RNA-seq data** using powerful gene set enrichment analysis.

## 🚀 Key Features

* **Enrichr Integration :**
    * Direct query of the [Enrichr API](https://maayanlab.cloud/Enrichr/) for gene set enrichment analysis.
    * Multi-threaded `get_cell_type_enrichment()` for fast cell type enrichment using curated gene set libraries.
* **Scanpy Integration :**
    * Designed to work seamlessly with [Scanpy](https://scanpy.readthedocs.io/en/stable/) `AnnData` objects.
    * Example workflow uses Scanpy for data loading, clustering, and marker gene identification.
* **Biomart Integration :**
    * Easily convert gene symbols across species using the Ensembl Biomart API.


## 📦 Installation

You can install the stable release from PyPI using pip:

```bash
pip install clustermolepy
```

Once installed, you can import and use `clustermolepy` in your Python environment.


## 🕹️ Quick Usage Example

Here's a simplified example of how to use `clustermolepy` to annotate cell clusters. For a more detailed walkthrough, check out the Jupyter Notebook in the `examples` directory\!

**Example: Interpreting a PBMC Cluster**

```python
import scanpy as sc
from clustermolepy.enrichr import Enrichr

# Load PBMC data and cluster
adata = sc.datasets.pbmc3k_processed()
sc.tl.leiden(adata, resolution=0.5, n_iterations=2, flavor='igraph')

# Get top marker genes for a cluster
sc.tl.rank_genes_groups(adata, 'leiden', method='wilcoxon')
markers = sc.get.rank_genes_groups_df(adata, group='1')
top_genes = markers['names'][:25].tolist()

# Run enrichment
enr = Enrichr(gene_list=top_genes)
df = enr.get_cell_type_enrichment()
print(df.head(3))
```


This will return a table with top matching cell types across multiple reference libraries like CellMarker, PanglaoDB, Azimuth, and more.

**(Example Output - first 3 rows of the table):**
```markdown
|    | term name               |     p-value |   odds ratio |   combined score | overlapping genes                                                                                            |   adjusted p-value |   old p-value |   old adjusted p-value | gene_set                 |
|---:|:------------------------|------------:|-------------:|-----------------:|:-------------------------------------------------------------------------------------------------------------|-------------------:|--------------:|-----------------------:|:-------------------------|
|  0 | B Cells Naive           | 9.79493e-24 |      571.108 |          30257.4 | ['CD79B', 'VPREB3', 'CD74', 'CD79A', 'FCER2', 'TCL1A', 'BANK1', 'LINC00926', 'LY86', 'CD37', 'LTB', 'MS4A1'] |        3.82002e-22 |             0 |                      0 | PanglaoDB_Augmented_2021 |
|  1 | B Cells                 | 9.91341e-23 |      466.235 |          23622.1 | ['CD79B', 'VPREB3', 'CD74', 'CD79A', 'FCER2', 'BANK1', 'HVCN1', 'LY86', 'CXCR4', 'CD37', 'LTB', 'MS4A1']     |        1.93312e-21 |             0 |                      0 | PanglaoDB_Augmented_2021 |
|  2 | B Cell Liver CL:0000236 | 7.54304e-22 |      622.531 |          30277.6 | ['CD79B', 'VPREB3', 'CD74', 'CD79A', 'BANK1', 'HVCN1', 'CXCR4', 'CD37', 'LTB', 'MS4A1']                      |        8.80014e-21 |             0 |                      0 | Tabula_Muris             |
```
**Example: Cross-Species Gene Mapping**

```python
from clustermolepy.utils import Biomart

bm = Biomart()
result = bm.convert_gene_names(
    genes=["TP53", "CD4", "FOXP3"],
    from_organism="hsapiens",
    to_organism="mmusculus"
)
print(result)
```

**(Example Output)**

```
{
  'TP53': ['Trp53'],
  'CD4': ['Cd4'],
  'FOXP3': ['Foxp3']
}
```
## 📚 Documentation

Check out the [Example Notebook](examples/clustermolepy_usage.ipynb) for more information!
