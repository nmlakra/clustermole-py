import typing
from functools import cache
from io import StringIO
from typing import Dict, List, Literal

import numpy as np
import pandas as pd
from biomart import BiomartServer

EnsemblOrgansimId = Literal["hsapiens", "mmusculus", "mfasicularis"]
ENSEMBLE_ORGANSIM_IDS = typing.get_args(EnsemblOrgansimId)
BIOMART_SERVER_URL = "http://www.ensembl.org/biomart"


@cache
def fetch_biomart_database(ensembl_dataset_id, url):
    server = BiomartServer(url)
    server.verbose = False
    return server.datasets[ensembl_dataset_id]


def convert_ensembl_genes():
    # TODO: make convert_gene_names and convert_ensembl_genes seperate wrappers around a generic function
    pass


def convert_gene_names(
    genes: List[str],
    current_organsim: EnsemblOrgansimId,
    target_organism: EnsemblOrgansimId,
    url: str = BIOMART_SERVER_URL,
) -> Dict[List[str], None]:
    """
    Converts gene names from one organism to another using Biomart.

    Args:
        genes (List[str]): List of gene names to convert.
        current_organsim (EnsemblOrgansimId): The Ensembl organism ID of the current gene names.
        target_organism (EnsemblOrgansimId): The Ensembl organism ID to convert the gene names to.
        url (str, optional): The URL of the Biomart server. Defaults to BIOMART_SERVER_URL.

    Returns:
        Dict: A dictionary where keys are the original gene names and values are the corresponding converted gene names.
    """
    assert (
        current_organsim != target_organism
    ), "current_organism and target_organism cannot be same"
    if not genes:
        raise ValueError("Invalid gene list provided")
    if current_organsim not in ENSEMBLE_ORGANSIM_IDS:
        raise ValueError("Invalid current_organsim")
    if target_organism not in ENSEMBLE_ORGANSIM_IDS:
        raise ValueError("Invalid target_organism")

    ensembl_dataset_id = f"{current_organsim}_gene_ensembl"
    database = fetch_biomart_database(ensembl_dataset_id, url)

    attributes = [
        "ensembl_gene_id",
        "external_gene_name",
        f"{target_organism}_homolog_associated_gene_name",
    ]

    search_response = database.search({"attributes": attributes})
    gene_encoding_df = pd.read_csv(
        StringIO(search_response.text), names=attributes, sep="\t"
    ).set_index("ensembl_gene_id")

    gene_encodings = (
        gene_encoding_df[gene_encoding_df["external_gene_name"].isin(genes)]
        .groupby("external_gene_name")
        .agg(list)
        .to_dict()[f"{target_organism}_homolog_associated_gene_name"]
    )
    gene_encodings = {
        key: value if value != [np.nan] else None
        for key, value in gene_encodings.items()
    }

    return gene_encoding
