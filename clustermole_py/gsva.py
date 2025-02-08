from gseapy import gp
from anndata import AnnData
import pandas as pd
import numpy as np
from scipy import sparse


class GSVA:

    def __init__(self):
        pass

    def aggregate_expression(
        self, adata: AnnData, group_by: str, scale_factor: int
    ) -> pd.DataFrame | None:

        count_matrix = adata.X
        groups = adata.obs[group_by].unique()

        if not groups:
            return None
        if not count_matrix:
            return None

        normalized_df = pd.DataFrame(columns=groups, index=adata.var.index)
        for group in groups:
            group_idx = adata.obs.reset_index()[adata.obs[group_by] == group][
                group_by
            ].tolist()
            total_counts = count_matrix[group_idx].sum()  # type: ignore
            normalized_vector = (count_matrix[group_idx].sum(axis=1) / total_counts) * scale_factor  # type: ignore

            normalized_df.loc[:, group] = normalized_vector

        return normalized_df
