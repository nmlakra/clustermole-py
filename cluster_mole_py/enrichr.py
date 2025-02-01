import requests
import pandas as pd
from typing import List, Dict
import json
from .utils import (
    DEFAULT_ALL_CELL_TYPE_LIBRARIES,
    DEFAULT_HS_CELL_TYPE_LIBRARIES,
    DEFAULT_MM_CELL_TYPE_LIBRARIES,
)


ENRICHR_URL = "https://maayanlab.cloud/Enrichr/"
SPEEDRICHR_URL = "https://maayanlab.cloud/speedrichr/"

def send_gene_list(gene_list: List[str], description: str) -> int:

    url = ENRICHR_URL + "addList"
    payload = {
            "list": (None, "\n".join(gene_list)),
            "description": (None, description)
    }

    response = requests.post(url, files=payload)

    if not response.ok:
        raise Exception(f'Error sending gene list, status code: {response.status_code}')

    data = json.loads(response.text)

    if not data['userListId']:
        raise ValueError('Could not submit gene list')

    return data['userListId']

def get_enrichment_results(gene_list: List[str], gene_set_library: str, description: str = "") -> Dict[str, List]:

    gene_list = [gene.upper() for gene in gene_list]
    url = ENRICHR_URL + "enrich"
    params = {
            'userListId': send_gene_list(gene_list, description),
            'backgroundType': gene_set_library
    }

    response = requests.get(url, params)

    if not response.ok:
        raise Exception(f'Error fetching enrichment results, status code: {response.status_code}')
    data = json.loads(response.text)
    return data
