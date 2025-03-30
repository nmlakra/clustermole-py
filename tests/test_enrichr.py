import unittest
from clustermolepy.enrichr import Enrichr
from pprint import pprint


class TestEnrichr(unittest.TestCase):

    def setUp(self):
        self.genes = ['PXK', 'MS4A1', 'CD19', 'CD74', 'CD79A']
        self.enrichr = Enrichr(self.genes)


    def test_get_libraries(self):
        libraries = self.enrichr.get_libraries("cell_types")
        pprint(libraries)


if __name__ == '__main__':
    unittest.main()

