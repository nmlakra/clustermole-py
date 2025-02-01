import requests
import pandas as pd
from typing import List, Dict
import json


class Enrichr:

    ENRICHR_URL = "https://maayanlab.cloud/Enrichr/"
    SPEEDRICHR_URL = "https://maayanlab.cloud/speedrichr/"

    def __init__(
        self,
        gene_list: List[str],
        pval_cutoff: float | None = None,
        adj_pval_cutoff: float | None = None,
        description: str = "",
    ):
        self.gene_list = gene_list
        self.description = description
        self.user_list_id = None
        self.pval_cutoff = pval_cutoff
        self.adj_pval_cutoff = adj_pval_cutoff

        if self.user_list_id is None:
            self.send_gene_list()

    def get_gene_list(self) -> List[str]:
        return [gene.upper() for gene in self.gene_list if gene.isalnum()]

    def send_gene_list(self) -> int:

        url = Enrichr.ENRICHR_URL + "addList"
        if self.gene_list:
            payload = {
                "list": (None, "\n".join(self.get_gene_list())),
                "description": (None, self.description),
            }
        else:
            raise ValueError("Empty gene list")

        response = requests.post(url, files=payload, verify=True)

        if not response.ok:
            raise Exception(
                f"Error sending gene list, status code: {response.status_code}"
            )

        data = json.loads(response.text)

        if not data["userListId"]:
            raise ValueError("Could not submit gene list")

        self.user_list_id = data["userListId"]
        return self.user_list_id

    def get_enrichment_results(self, gene_set: str) -> pd.DataFrame:

        url = Enrichr.ENRICHR_URL + "enrich"

        params = {"userListId": self.user_list_id, "backgroundType": gene_set}

        response = requests.get(url, params)

        if not response.ok:
            raise Exception(
                f"Error fetching enrichment results, status code: {response.status_code}"
            )
        data = json.loads(response.text)

        return self.format_enrichment_result(data, gene_set)

    def format_enrichment_result(self,
        enrichment_result: Dict, gene_set: str
    ) -> pd.DataFrame:
        cols = [
            "rank",
            "term name",
            "p-value",
            "odds ratio",
            "combined score",
            "overlapping genes",
            "adjusted p-value",
            "old p-value",
            "old adjusted p-value",
        ]
        df = pd.DataFrame(enrichment_result[gene_set], columns=cols)
        df["gene_set"] = gene_set

        if self.adj_pval_cutoff:
            return df[df["adjusted p-value"] <= self.adj_pval_cutoff]
        elif self.pval_cutoff:
            return df[df["p-value"] <= self.pval_cutoff]
        return df
