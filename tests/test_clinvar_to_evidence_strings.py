import unittest

from eva_cttv_pipeline import clinvar_to_evidence_strings
from eva_cttv_pipeline import utilities


class GetTermsFromFile(unittest.TestCase):
    def setUp(self):
        ignore_file = utilities.get_resource_file("eva_cttv_pipeline", "resources/testing/ignore_file_testing.txt")
        self.ignore_terms = clinvar_to_evidence_strings.get_terms_from_file(ignore_file)

    def test_length(self):
        self.assertEqual(len(self.ignore_terms), 218)

    def test_head(self):
        self.assertEqual(self.ignore_terms[0], "http://purl.obolibrary.org/obo/HP_0011677")

    def test_tail(self):
        self.assertEqual(self.ignore_terms[-1], "http://www.orpha.net/ORDO/Orphanet_120795")

    def test_no_file(self):
        self.assertEqual(clinvar_to_evidence_strings.get_terms_from_file(None), [])
