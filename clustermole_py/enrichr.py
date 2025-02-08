import requests
import pandas as pd
from typing import Iterable, List, Dict, Literal
import json
from concurrent.futures import ThreadPoolExecutor


GeneSetLibrary = Literal[
    "CellMarker_2024",
    "CellMarker_Augmented_2021",
    "Descartes_Cell_Types_and_Tissue_2021",
    "PanglaoDB_Augmented_2021",
    "Azimuth_Cell_Types_2021",
    "Azimuth_2023",
    "Tabula_Sapiens",
    "Human_Gene_Atlas",
    "Tabula_Muris",
    "Mouse_Gene_Atlas",
]


class Enrichr:

    enrichr_url = "https://maayanlab.cloud/Enrichr/"
    speedrichr_url = "https://maayanlab.cloud/speedrichr/"

    default_cell_type_libraries: Iterable[GeneSetLibrary] = [
        "CellMarker_2024",
        "CellMarker_Augmented_2021",
        "Descartes_Cell_Types_and_Tissue_2021",
        "PanglaoDB_Augmented_2021",
        "Azimuth_Cell_Types_2021",
        "Azimuth_2023",
        "Tabula_Sapiens",
        "Human_Gene_Atlas",
        "Tabula_Muris",
        "Mouse_Gene_Atlas",
    ]

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

        url = Enrichr.enrichr_url + "addList"

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

    def get_enrichment(
        self, gene_set: GeneSetLibrary, sort_order: str = "adjusted p-value"
    ) -> pd.DataFrame | pd.Series:

        url = Enrichr.enrichr_url + "enrich"

        params = {"userListId": self.user_list_id, "backgroundType": gene_set}
        response = requests.get(url, params)
        if not response.ok:
            raise Exception(
                f"Error fetching enrichment results, status code: {response.status_code}"
            )
        data = json.loads(response.text)
        result = self.format_enrichment_result(data, gene_set, sort_order)
        if result is None:
            raise Exception("Could not get enrichment results")
        return result

    def format_enrichment_result(
        self,
        enrichment_result: Dict,
        gene_set: GeneSetLibrary,
        sort_order: str | List[str],
    ) -> pd.DataFrame | pd.Series | None:
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
        df = pd.DataFrame(enrichment_result[gene_set], columns=pd.Index(cols))
        df["gene_set"] = gene_set

        if self.adj_pval_cutoff:
            df = df[df["adjusted p-value"] <= self.adj_pval_cutoff]
        elif self.pval_cutoff:
            df = df[df["p-value"] <= self.pval_cutoff]
        if sort_order in df.columns:
            return df.sort_values(by=sort_order).drop("rank", axis=1)  # type: ignore
        else:
            raise KeyError(f"{sort_order} not found")

    def get_cell_type_enrichment(
        self,
        gene_sets: Iterable[GeneSetLibrary] | None = None,
        max_workers: int | None = None,
        sort_order: str = "adjusted p-value",
    ):

        if gene_sets is None:
            gene_sets = Enrichr.default_cell_type_libraries
        with ThreadPoolExecutor(max_workers=max_workers):
            results = (
                pd.concat(  # type: ignore
                    map(
                        lambda gene_set: self.get_enrichment(gene_set, sort_order),
                        gene_sets,
                    )
                )
                .sort_values(by=sort_order)
                .reset_index(drop=True)
            )

        return results

    @staticmethod
    def load_gene_set(gene_set: GeneSetLibrary) -> Dict[str, Dict[str, List]]:
        url = Enrichr.enrichr_url + "geneSetLibrary"
        params = {"mode": "json", "libraryName": gene_set}

        response = requests.get(url, params)
        if not response.ok:
            raise Exception(
                f"Error fetching enrichment results, status code: {response.status_code}"
            )
        data = json.loads(response.text)
        if not data[gene_set]:
            raise Exception(f"Failed fetching the correct {gene_set}")

        return {gene_set: data["terms"]}
