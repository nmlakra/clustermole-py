# cluster-mole-py
# 🧬 clustermole_py: Cluster Annotation for Single-Cell Data in Python (Pre-release)

[![License](https://img.shields.io/badge/License-MIT-blue.svg)](https://opensource.org/licenses/MIT)

`clustermole_py` is a Python package (currently in **pre-release**) inspired by the original `clustermole` R package. It's designed to help you **annotate cell clusters from single-cell RNA-seq data** 📊 using powerful gene set enrichment analysis.

✨ Give your clusters biological meaning with **Enrichr** and **GSVA**! ✨

**⚠️  Important: `clustermole_py` is currently under development and not yet available on PyPI. Installation instructions below are for local development or installation from a source repository.  Stay tuned for the PyPI release! ⚠️**

## 🚀 Key Features

* **Enrichr Integration :**
    * Direct query of the [Enrichr API](https://maayanlab.cloud/Enrichr/) for gene set enrichment analysis.
    * Specialized `Enrichr` module with methods like `get_enrichment()` and `get_cell_type_enrichment()`.
    *  ⚡️ **Multi-threaded `get_cell_type_enrichment()`** for fast cell type annotation using curated libraries.
    *  Supports a wide range of Enrichr gene set libraries.
* **GSVA Module :**
    *  Perform **Gene Set Variation Analysis (GSVA)** directly on your cluster data.
    *  `aggregate_expression()` function for easy normalization of cluster-level expression data from Scanpy `AnnData` objects.
    *  `gsva.get_enrichment()` for DE-free cluster annotation.
* **Scanpy Integration :**
    *  Designed to work seamlessly with [Scanpy](https://scanpy.readthedocs.io/en/stable/) `AnnData` objects.
    *  Example workflow uses Scanpy for data loading, clustering, and marker gene identification.
* **Clear and Interpretable Results :**
    *  Returns enrichment results in easy-to-read tables (Pandas DataFrames - *if applicable to your package*).
    *  Output tables include p-values, adjusted p-values, combined scores, and overlapping genes for enriched terms.

## 📦 Installation (Development/Local Installation - Not on PyPI yet)

For this pre-release version, you can install `clustermole_py` directly from GitHub using pip:

```bash
!pip install scanpy igraph leidenalg git+[https://github.com/nmlakra/clustermole-py.git](https://github.com/nmlakra/clustermole-py.git)
```
