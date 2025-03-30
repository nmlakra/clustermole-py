import unittest
from pprint import pprint

from clustermolepy.enrichr import Enrichr


class TestEnrichr(unittest.TestCase):

    def setUp(self):
        self.genes = ["PXK", "MS4A1", "CD19", "CD74", "CD79A"]
        self.enrichr = Enrichr(self.genes)

    def test_get_libraries_name(self):
        libraries = self.enrichr.get_libraries(name="KEGG")
        pprint(libraries)

    def test_get_libraries_category(self):
        libraries = self.enrichr.get_libraries(category="cell_types")
        pprint(libraries)


if __name__ == "__main__":
    unittest.main()
