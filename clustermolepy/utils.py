import json
import os
from io import StringIO
from json.decoder import JSONDecodeError
from typing import Dict, List

import pandas as pd
from biomart import BiomartServer
from platformdirs import user_cache_dir

APP_NAME = "clustermolepy"
CACHE_FILE = os.path.join(user_cache_dir(APP_NAME), "valid_organisms.json")

BIOMART_SERVER_URL = "http://www.ensembl.org/biomart"


class Biomart:

    def __init__(
        self,
        url: str = BIOMART_SERVER_URL,
        verbose: bool = True,
        cache_file: str = CACHE_FILE,
    ):
        """
        Initializes the Biomart class.

        Args:
            url: URL of the Biomart server.
            verbose: If True, prints additional information.
            cache_file: Path to the cache file for valid organisms.
        """
        self.url = url
        self.server = BiomartServer(url)
        self.server.verbose = verbose
        self.cache_file = cache_file
        self.valid_organisms = self._load_valid_organisms()

    def fetch_biomart_database(self, ensembl_dataset_id: str) -> BiomartServer:
        """
        Fetches the Biomart database for a given Ensembl dataset ID.

        Args:
            ensembl_dataset_id: Ensembl dataset ID.
        Returns:
            BiomartServer object for the specified dataset.
        Raises:
            KeyError: If the dataset ID is not valid.
        """
        return self.server.datasets[ensembl_dataset_id]

    def _load_valid_organisms(self) -> set:
        """
        Load the set of valid organisms from a cache file.

        This method checks if the specified cache file exists. If it does, it attempts
        to load the contents of the file as a JSON object and convert it into a set.
        If the file cannot be decoded as JSON, an empty set is returned. If the cache
        file does not exist, an empty set is also returned.

        Returns:
            set: A set of valid organisms loaded from the cache file, or an empty set
            if the file does not exist or cannot be decoded.
        """
        if os.path.exists(self.cache_file):
            try:
                with open(self.cache_file) as f:
                    return set(json.load(f))
            except JSONDecodeError:
                return set()
        return set()

    def _save_valid_organisms(self):
        """
        Saves the list of valid organisms to a cache file.

        This method ensures that the directory for the cache file exists,
        and then writes the sorted list of valid organisms to the file in JSON format.

        Raises:
            OSError: If there is an issue creating the directory or writing to the file.
        """
        os.makedirs(os.path.dirname(self.cache_file), exist_ok=True)
        with open(self.cache_file, "w") as f:
            json.dump(sorted(self.valid_organisms), f)

    def _validate_organsim(self, organism: str):
        """
        Validates whether the given organism is supported.

        This method checks if the provided organism is in the list of valid organisms.
        If not, it attempts to query the Biomart server to verify its availability.
        If the organism is found in Biomart, it is added to the list of valid organisms
        and the updated list is saved. Otherwise, a ValueError is raised.

        Args:
            organism (str): The name of the organism to validate.

        Raises:
            ValueError: If the organism is not valid or not available in Biomart.
        """

        if organism in self.valid_organisms:
            return

        # Try querying in biomart
        try:
            _ = self.server.datasets[f"{organism}_gene_ensembl"]
            self.valid_organisms.add(organism)
            self._save_valid_organisms()
        except KeyError:
            raise ValueError(
                f"Organism `{organism}` is not a valid or available in Biomart"
            )

    def _convert(
        self,
        query_values: List[str],
        from_organism: str,
        to_organism: str,
        filter_name: str,
        query_column: str,
        return_ensembl_ids: bool = False,
    ) -> Dict[str, List[str | None]]:
        """
        Internal method to convert gene symbols or Ensembl IDs using Biomart.

        Args:
            query_values: List of gene symbols or Ensembl IDs.
            from_organism: Source organism Ensembl ID.
            to_organism: Target organism Ensembl ID.
            filter_name: Name of the filter to apply in Biomart.
            query_column: Column to group results on.
            return_ensembl_ids: If True, returns homolog Ensembl IDs instead of gene symbols.

        Returns:
            Dict mapping input IDs to homolog gene symbols or Ensembl IDs in the target organism.
        """
        if from_organism == to_organism:
            raise ValueError("from_organism and to_organism cannot be the same")
        if not query_values:
            raise ValueError("Input query list is empty")
        self._validate_organsim(from_organism)
        self._validate_organsim(to_organism)

        result = {val: [] for val in query_values}
        dataset_id = f"{from_organism}_gene_ensembl"
        database = self.fetch_biomart_database(dataset_id)

        # Choose target column based on whether to return ensembl IDs or gene names
        target_attr = (
            f"{to_organism}_homolog_ensembl_gene"
            if return_ensembl_ids
            else f"{to_organism}_homolog_associated_gene_name"
        )

        attributes = ["ensembl_gene_id", "external_gene_name", target_attr]
        filters = {filter_name: query_values}

        try:
            response = database.search({"attributes": attributes, "filters": filters})
            df = pd.read_csv(StringIO(response.text), sep="\t", names=attributes)
        except Exception as e:
            raise RuntimeError(f"Biomart query failed: {e}") from e

        if df.empty:
            return result

        grouped = (
            df.groupby(query_column)[target_attr]
            .apply(lambda x: [v for v in x if not pd.isna(v)])
            .to_dict()
        )

        for val, homologs in grouped.items():
            result[val] = homologs if homologs else []

        return result

    def convert_ensembl_ids(
        self,
        ensembl_ids: List[str],
        from_organism: str,
        to_organism: str,
    ) -> Dict[str, List[str | None]]:
        """
        Converts Ensembl gene IDs from one organism to Ensembl IDs of homologs in another organism.
        """
        return self._convert(
            query_values=ensembl_ids,
            from_organism=from_organism,
            to_organism=to_organism,
            filter_name="ensembl_gene_id",
            query_column="ensembl_gene_id",
            return_ensembl_ids=True,
        )

    def convert_gene_names(
        self,
        genes: List[str],
        from_organism: str,
        to_organism: str,
    ) -> Dict[str, List[str | None]]:
        """
        Converts gene names from one organism to gene names of homologs in another organism.
        """
        return self._convert(
            query_values=genes,
            from_organism=from_organism,
            to_organism=to_organism,
            filter_name="external_gene_name",
            query_column="external_gene_name",
        )
