import json
from concurrent.futures import ThreadPoolExecutor
from functools import lru_cache
from typing import Dict, Iterable, List, Literal

import pandas as pd
import requests

GeneSetCategory = Literal[
    "cell_types",
    "crowd",
    "diseases_drugs",
    "legacy",
    "misc",
    "ontologies",
    "pathways",
    "transcription",
]

CATEGORY_NAME_MAP = {
    "cell_types": "Cell Types",
    "crowd": "Crowd",
    "diseases_drugs": "Diseases/Drugs",
    "legacy": "Legacy",
    "misc": "Misc",
    "ontologies": "Ontologies",
    "pathways": "Pathways",
    "transcription": "Transcription",
}


class Enrichr:

    enrichr_url = "https://maayanlab.cloud/Enrichr/"
    speedrichr_url = "https://maayanlab.cloud/speedrichr/"

    _valid_libraries: set[str] | None = None

    default_cell_type_libraries = [
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
        """
        Initializes an Enrichr object and submits the gene list to the Enrichr API.

        Args:
            gene_list: List of gene symbols to analyze.
            pval_cutoff: Optional unadjusted p-value cutoff to filter enrichment results.
            adj_pval_cutoff: Optional adjusted p-value cutoff to filter enrichment results.
            description: Optional description of the gene list.
        Raises:
            ValueError: If the gene list is empty or submission fails.
        """

        self.gene_list = gene_list
        self.description = description
        self.user_list_id = None
        self.pval_cutoff = pval_cutoff
        self.adj_pval_cutoff = adj_pval_cutoff

        if self.user_list_id is None:
            self.send_gene_list()

    def get_gene_list(self) -> List[str]:
        """
        Returns a cleaned version of the gene list in uppercase.

        Returns:
            A list of gene symbols in uppercase, excluding non-alphanumeric entries.
        """
        return [gene.upper() for gene in self.gene_list if gene.isalnum()]

    def send_gene_list(self) -> int:
        """
        Submits the gene list to Enrichr and retrieves a user list ID.

        Args:
            gene_list: List of gene symbols to analyze.
            description: Optional description of the gene list.
        Raises:
            ValueError: If the gene list is empty or submission fails.
        Returns:
            user_list_id: Unique identifier for the submitted gene list.
        """
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
        self, gene_set: str, sort_order: str = "adjusted p-value"
    ) -> pd.DataFrame:
        """
        Fetches enrichment results for a specific gene set.

        Args:
            gene_set: The gene set library to use for enrichment analysis.
            sort_order: Column name to sort the results by.
        Returns:
            A DataFrame containing the enrichment results.
        Raises:
            Exception: If the request fails or if the results cannot be formatted.
        """
        if gene_set not in self.get_valid_libraries():
            raise ValueError(f"Invalid gene set library name: {gene_set}")

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
        gene_set: str,
        sort_order: str | List[str],
    ) -> pd.DataFrame | None:
        """
        Formats and filters enrichment results based on p-value cutoffs and sorting preferences.

        Args:
            enrichment_result: The raw enrichment results from the API.
            gene_set: The gene set library used for enrichment analysis.
            sort_order: Column name or list of column names to sort the results by.
        Returns:
            A DataFrame containing the formatted enrichment results.
        Raises:
            KeyError: If the specified sort order is not found in the results.
        """
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
        gene_sets: Iterable[str] | None = None,
        max_workers: int | None = None,
        sort_order: str = "adjusted p-value",
    ):
        """
        Runs enrichment analysis for multiple gene sets using multithreading.

        Args:
            gene_sets: List of gene set libraries to analyze.
            max_workers: Maximum number of threads to use for parallel processing.
            sort_order: Column name to sort the results by.
        Returns:
            A DataFrame containing the combined enrichment results for all gene sets.
        Raises:
            ValueError: If the gene sets are empty or invalid.
        """
        if gene_sets is None:
            gene_sets = Enrichr.default_cell_type_libraries

        invalid_gene_sets = [
            gs for gs in gene_sets if gs not in self.get_valid_libraries()
        ]
        if invalid_gene_sets:
            raise ValueError(f"Invalid gene set libraries: {invalid_gene_sets}")

        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            results = (
                pd.concat(  # type: ignore
                    executor.map(
                        lambda gene_set: self.get_enrichment(gene_set, sort_order),
                        gene_sets,
                    )
                )
                .sort_values(by=sort_order)
                .reset_index(drop=True)
            )

        return results

    @staticmethod
    def load_gene_set(gene_set: str) -> Dict[str, Dict[str, List]]:
        """
        Loads the gene set library from Enrichr.

        Args:
            gene_set: The gene set library to load.
        Returns:
            A dictionary containing the gene set library.
        Raises:
            Exception: If the request fails or if the results cannot be formatted.
        """
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

    @classmethod
    def _set_valid_libraries(cls, library_metadata: Dict[str, List]) -> None:
        """
        Sets the valid libraries based on the fetched metadata.
        This method is called internally to initialize the valid libraries.

        Args:
            library_metadata: The metadata containing library statistics and categories.
        """
        library_details = library_metadata["statistics"]
        cls._valid_libraries = {lib["libraryName"] for lib in library_details}

    @classmethod
    def get_valid_libraries(cls) -> set[str]:
        """
        This method fetches the library metadata if it has not been previously set.
        It then sets the valid libraries based on the fetched metadata.

        Returns:
            A set of valid library names.
        Raises:
            Exception: If the request fails or if the results cannot be formatted.
        """
        if cls._valid_libraries is None:
            library_metadata = cls.fetch_libraries()
            cls._set_valid_libraries(library_metadata)
        assert cls._valid_libraries is not None
        return cls._valid_libraries

    @staticmethod
    @lru_cache(maxsize=1)
    def fetch_libraries() -> Dict[str, List]:
        """
        Fetches the list of available libraries from Enrichr.

        Returns:
            A dictionary containing the library statistics and categories.
        Raises:
            Exception: If the request fails or if the results cannot be formatted.
        """
        url = Enrichr.enrichr_url + "datasetStatistics"
        response = requests.get(url)

        if response.status_code == 200:
            library_metadata = response.json()
            return library_metadata
        else:
            raise Exception(f"Failed to fetch libraries: {response.status_code}")

    @staticmethod
    def get_libraries(
        name: str | None = None, category: GeneSetCategory | None = None
    ) -> pd.DataFrame:
        """
        Fetches the list of available libraries from Enrichr and filters them based on the provided name or category.

        Args:
            name: Optional name or partial name of the library to filter.
            category: Optional category to filter the libraries.
        Returns:
            A DataFrame containing the filtered library information.
        Raises:
            ValueError: If the provided category is invalid or if no libraries match the name.
        """

        library_metadata = Enrichr.fetch_libraries()
        library_info = library_metadata["statistics"]
        category_info = library_metadata["categories"]

        library_df = pd.DataFrame(library_info).drop(["appyter"], axis=1)

        # Fuzzy matching the get set library names
        if name is not None:
            library_df = library_df[
                library_df["libraryName"].str.contains(name, case=False)
            ].reset_index(drop=True)

            if library_df.empty:
                raise ValueError(f"{name} could not be matched with any library")

        # Filtering results by category
        if category is not None:
            if category not in CATEGORY_NAME_MAP:
                raise ValueError(f"Invalid category: {category}")

            category_name = CATEGORY_NAME_MAP[category]
            category_id = next(
                (c["categoryId"] for c in category_info if c["name"] == category_name),
                None,
            )

            if category_id is None:
                raise ValueError(f"{category} could not be mapped to a valid Id")

            library_df = library_df.query("categoryId == @category_id").reset_index(
                drop=True
            )

        assert isinstance(library_df, pd.DataFrame)
        return library_df
