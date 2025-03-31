import json
import os
import tempfile
import unittest
from unittest.mock import MagicMock

from clustermolepy.utils import Biomart


class TestBiomart(unittest.TestCase):
    def setUp(self):
        # Use a temporary cache file so we don't touch ~/.cache
        self.temp_cache = tempfile.NamedTemporaryFile(delete=False)
        with open(self.temp_cache.name, "w") as f:
            json.dump([], f)

        self.mock_server = MagicMock()
        self.biomart = Biomart(cache_file=self.temp_cache.name)
        self.biomart.server = self.mock_server

        self.hsapiens_dataset = MagicMock()
        self.mock_server.datasets = {
            "hsapiens_gene_ensembl": self.hsapiens_dataset,
            "mmusculus_gene_ensembl": self.hsapiens_dataset,
        }

        self.mock_response_text = (
            "ENSG00000141510\tTP53\tTrp53\n" "ENSG00000012048\tBRCA1\tBrca1\n"
        )

        self.hsapiens_dataset.search.return_value.text = self.mock_response_text

    def tearDown(self):
        os.unlink(self.temp_cache.name)

    def test_convert_gene_names_success(self):
        genes = ["TP53", "BRCA1"]
        result = self.biomart.convert_gene_names(genes, "hsapiens", "mmusculus")
        self.assertIn("TP53", result)
        self.assertIsInstance(result["TP53"], list)
        self.assertTrue("Trp53" in result["TP53"])

    def test_convert_ensembl_ids_success(self):
        ensembl_ids = ["ENSG00000141510", "ENSG00000012048"]
        result = self.biomart.convert_ensembl_ids(ensembl_ids, "hsapiens", "mmusculus")
        self.assertIn("ENSG00000141510", result)
        self.assertTrue("Trp53" in result["ENSG00000141510"])

    def test_invalid_organism_raises(self):
        with self.assertRaises(ValueError):
            self.biomart.convert_gene_names(["TP53"], "hsapiens", "fakeorganism")

    def test_empty_gene_list_raises(self):
        with self.assertRaises(ValueError):
            self.biomart.convert_gene_names([], "hsapiens", "mmusculus")

    def test_same_organism_raises(self):
        with self.assertRaises(ValueError):
            self.biomart.convert_ensembl_ids(
                ["ENSG00000141510"], "hsapiens", "hsapiens"
            )

    def test_cache_file_writes(self):
        # Make sure it saves validated organisms
        self.biomart._validate_organsim("hsapiens")
        with open(self.temp_cache.name) as f:
            data = json.load(f)
            self.assertIn("hsapiens", data)


if __name__ == "__main__":
    unittest.main()
