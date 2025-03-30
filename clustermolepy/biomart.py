import typing
from io import StringIO
from typing import Dict, List, Literal

import pandas as pd
from biomart import BiomartServer

EnsemblOrganismId = Literal[
    "hsapiens",  # Human
    "mmusculus",  # Mouse
    "rnorvegicus",  # Rat
    "dmelanogaster",  # Fruit fly
    "celegans",  # Nematode
    "scerevisiae",  # Yeast
    "athaliana",  # Thale cress
    "drerio",  # Zebrafish
    "ggallus",  # Chicken
    "btaurus",  # Cow
    "sscrofa",  # Pig
    "cfamiliaris",  # Dog
    "ecaballus",  # Horse
    "oaries",  # Sheep
    "mfascicularis",  # Crab-eating macaque
    "ptroglodytes",  # Chimpanzee
    "mdomestica",  # Opossum
    "xtropicalis",  # Frog
    "olatipes",  # Medaka fish
    "trubripes",  # Pufferfish
    "cintestinalis",  # Sea squirt
    "spombe",  # Fission yeast
    "mabscessus",  # Bacteria
    "pfalciparum",  # Malaria parasite
    "lmajor",  # Parasite
    "tbrucei",  # Parasite
    "ecoli",  # Bacteria
    "pberghei",  # Malaria parasite
    "senterica",  # Bacteria
    "paeruginosa",  # Bacteria
]
ENSEMBLE_ORGANISM_IDS = typing.get_args(EnsemblOrganismId)
BIOMART_SERVER_URL = "http://www.ensembl.org/biomart"


class Biomart:

    def __init__(self, url: str = BIOMART_SERVER_URL, verbose: bool = True):
        self.server = BiomartServer(url)
        self.server.verbose = verbose

    def fetch_biomart_database(self, ensembl_dataset_id):
        return self.server.datasets[ensembl_dataset_id]

    def convert_ensembl_genes(self):
        # TODO: make convert_gene_names and convert_ensembl_genes seperate wrappers around a generic function
        pass

    def convert_gene_names(
        self,
        genes: List[str],
        from_organism: EnsemblOrganismId,
        to_organism: EnsemblOrganismId,
    ) -> Dict[str, List[str | None]]:
        """
        Converts gene symbols from one organism to another using Biomart.

        Args:
            genes: A list of gene symbols (e.g., ["TP53", "BRCA1"]) to convert.
            from_organism: The Ensembl organism ID of the source organism (e.g., "hsapiens" for human).
            to_organism: The Ensembl organism ID of the target organism (e.g., "mmusculus" for mouse).

        Returns:
            A dictionary where:
            - Keys are the original gene symbols from the `genes` input.
            - Values are lists of corresponding homolog gene symbols in the target organism.
              If no homolog is found, the value is `[None]`.

        Raises:
            ValueError: If `from_organism` and `to_organism` are the same, or if `genes` is empty.
            RuntimeError: If the Biomart query fails or returns invalid data.
        """
        assert (
            from_organism != to_organism
        ), "from_organism and to_organism cannot be same"
        if not genes:
            raise ValueError("Invalid gene list provided")
        if from_organism not in ENSEMBLE_ORGANISM_IDS:
            raise ValueError("Invalid from_organism")
        if to_organism not in ENSEMBLE_ORGANISM_IDS:
            raise ValueError("Invalid to_organism")

        # Initialize result with all genes -> []
        result = {gene: [] for gene in genes}

        ensembl_dataset_id = f"{from_organism}_gene_ensembl"
        database = self.fetch_biomart_database(ensembl_dataset_id)
        attributes = [
            "ensembl_gene_id",
            "external_gene_name",
            f"{to_organism}_homolog_associated_gene_name",
        ]
        filters = {"external_gene_name": genes}

        try:
            response = database.search({"attributes": attributes, "filters": filters})
            df = pd.read_csv(StringIO(response.text), sep="\t", names=attributes)
        except Exception as e:
            raise RuntimeError(f"Biomart query failed: {e}") from e

        if df.empty:
            return result  # No homologs found

        # Group by gene and collect homologs, filtering NaNs
        grouped = (
            df.groupby("external_gene_name")[
                f"{to_organism}_homolog_associated_gene_name"
            ]
            .apply(lambda x: [v for v in x if not pd.isna(v)])
            .to_dict()
        )

        # Update result with found homologs
        for gene, homologs in grouped.items():
            result[gene] = homologs if homologs else []

        return result
