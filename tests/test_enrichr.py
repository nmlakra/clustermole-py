import unittest
from unittest.mock import patch, MagicMock
from clustermolepy.enrichr import Enrichr


class TestEnrichr(unittest.TestCase):

    @patch("clustermolepy.enrichr.requests.post")
    def test_send_gene_list_success(self, mock_post):
        mock_post.return_value.ok = True
        mock_post.return_value.text = '{"userListId": 12345}'

        enrichr = Enrichr(gene_list=["CD3E", "CD4"])
        self.assertEqual(enrichr.user_list_id, 12345)

    @patch("clustermolepy.enrichr.requests.post")
    def test_send_gene_list_failure(self, mock_post):
        mock_post.return_value.ok = False

        with self.assertRaises(Exception):
            Enrichr(gene_list=["CD3E", "CD4"])

    def test_get_gene_list_filters_invalid(self):
        enrichr = Enrichr.__new__(Enrichr)
        enrichr.gene_list = ["CD3E", "CD4", "!!!", " "]
        result = enrichr.get_gene_list()
        self.assertEqual(result, ["CD3E", "CD4"])

    def test_empty_gene_list_raises(self):
        with self.assertRaises(ValueError):
            Enrichr(gene_list=[])

    @patch("clustermolepy.enrichr.Enrichr.fetch_libraries")
    @patch("clustermolepy.enrichr.requests.get")
    @patch("clustermolepy.enrichr.requests.post")
    def test_get_enrichment_and_format(self, mock_post, mock_get, mock_fetch_libs):
        mock_post.return_value.ok = True
        mock_post.return_value.text = '{"userListId": 1}'

        mock_get.return_value.ok = True
        mock_get.return_value.text = """
        {
            "Test_Library": [
                [1, "termA", 0.01, 1.5, 2.0, ["GENE1"], 0.05, 0.01, 0.05]
            ]
        }
        """
        mock_fetch_libs.return_value = {
            "statistics": [
                {"libraryName": "Test_Library", "categoryId": 1, "appyter": "dummy"}
            ],
            "categories": [{"name": "Cell Types", "categoryId": 1}],
        }

        enrichr = Enrichr(gene_list=["GENE1"])
        df = enrichr.get_enrichment("Test_Library")

        self.assertEqual(df.shape[0], 1)
        self.assertIn("term name", df.columns)

    @patch("clustermolepy.enrichr.requests.get")
    def test_fetch_libraries(self, mock_get):
        mock_get.return_value.ok = True
        mock_get.return_value.status_code = 200
        mock_get.return_value.json.return_value = {
            "statistics": [
                {"libraryName": "TestLib", "categoryId": 1, "appyter": "something"}
            ],
            "categories": [{"name": "Cell Types", "categoryId": 1}],
        }

        result = Enrichr.fetch_libraries()
        self.assertIn("statistics", result)

    @patch("clustermolepy.enrichr.Enrichr.fetch_libraries")
    def test_get_libraries_by_name(self, mock_fetch):
        mock_fetch.return_value = {
            "statistics": [
                {"libraryName": "TestLib", "categoryId": 1, "appyter": "some_url"}
            ],
            "categories": [{"name": "Cell Types", "categoryId": 1}],
        }

        df = Enrichr.get_libraries(name="TestLib")
        self.assertEqual(df.iloc[0]["libraryName"], "TestLib")

    @patch("clustermolepy.enrichr.Enrichr.fetch_libraries")
    def test_get_libraries_invalid_category(self, mock_fetch):
        mock_fetch.return_value = {
            "statistics": [
                {"libraryName": "TestLib", "categoryId": 1, "appyter": "dummy"}
            ],
            "categories": [{"name": "Cell Types", "categoryId": 1}],
        }

        with self.assertRaises(ValueError) as cm:
            Enrichr.get_libraries(category="not_a_category")

        self.assertIn("Invalid category", str(cm.exception))


if __name__ == "__main__":
    unittest.main()
