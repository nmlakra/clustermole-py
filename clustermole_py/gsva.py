import gseapy as gp
from anndata import AnnData
import pandas as pd
import numpy as np
from typing import List


def get_enrichment(data: pd.DataFrame, gene_sets: str | List[str]) -> pd.DataFrame:
    enrichment_df = gp.gsva(data, gene_sets, outdir=None).res2d
    if enrichment_df is None:
        return pd.DataFrame()
    if isinstance(gene_sets, str):
        enrichment_df["gene_set"] = gene_sets
    else:
        temp_split = enrichment_df["Term"].str.split("__", expand=True)
        enrichment_df["Term"] = temp_split[1]
        enrichment_df["gene_set"] = temp_split[0]
    enrichment_df.rename(
        columns={"Name": "group", "Term": "term", "ES": "enrichment_score"},
        inplace=True,
    )

    return enrichment_df


def aggregate_expression(
    adata: AnnData,
    groupby: str,
    scale_factor: int,
    layer: str | None = None,
    use_raw: bool = False,
    apply_log1p: bool = True,
) -> pd.DataFrame | None:

    if layer and use_raw:
        raise Exception("Cannot have `layer` and have `use_raw=True`")

    count_matrix = None
    if use_raw:
        count_matrix = adata.raw.X
    elif layer is not None:
        assert adata.layers[layer] is not None, "Invailid `layer` key"
        count_matrix = adata.layers[layer]
    else:
        count_matrix = adata.X

    groups = adata.obs[groupby].unique()

    if not groups:
        return None
    if count_matrix is None:
        return None

    normalized_df = pd.DataFrame(columns=groups, index=adata.var.index)
    for group in groups:
        group_idx = adata.obs[adata.obs[groupby] == group].reset_index().index.tolist()
        total_counts = count_matrix[group_idx].sum()  # type: ignore
        normalized_vector = (count_matrix[group_idx].sum(axis=0) / total_counts) * scale_factor  # type: ignore
        normalized_df.loc[:, group] = normalized_vector

    if apply_log1p:
        return normalized_df.map(np.log1p)

    return normalized_df
