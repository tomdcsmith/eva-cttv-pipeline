import unittest

from eva_cttv_pipeline import clinvar_to_evidence_strings
from eva_cttv_pipeline import utilities


class GetTermsFromFileTest(unittest.TestCase):
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


class GetUnmappedUrlTest(unittest.TestCase):
    def test_orphanet(self):
        url = "http://www.orpha.net/ORDO/Orphanet_2670"
        self.assertEqual(clinvar_to_evidence_strings.get_unmapped_url(url), "http://purl.bioontology.org/ORDO/Orphanet_2670")

    def test_hp(self):
        url = "http://purl.obolibrary.org/obo/HP_0000545"
        self.assertEqual(clinvar_to_evidence_strings.get_unmapped_url(url), "http://purl.bioontology.org/obo/HP_0000545")

    def test_bad_url(self):
        #TODO currently the function exits execution with bad urls
        pass


class GetCttvVariantTypeTest(unittest.TestCase):
    def setUp(self):
        self.record_single_a = ({"reference": "AGAGACGTACGTACGTACGTACGTACGTACGTACGTACG", "alternate": "C"}, "snp single")
        self.record_single_b = ({"reference": "A", "alternate": "C"}, "snp single")
        self.record_single_c = ({"reference": "AGAGACGTACGTACGTACGTACGTACGTACGTACGTACG", "alternate": "AGAGACGTACGTACGTACGTACGTACGTACGTACGTACG"}, "snp single")

        self.test_records_singles = [self.record_single_a, self.record_single_b, self.record_single_c]

        self.record_structurals_a = ({"reference": "AGAGACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT", "alternate": "C"},
                                "structural variant")
        self.record_structurals_b = ({"reference": "A", "alternate": "AGAGACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"},
                                "structural variant")
        self.record_structurals_c = ({"reference": "AGAGACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT",
                                 "alternate": "AGAGACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"},
                                "structural variant")

        self.test_records_structurals = [self.record_structurals_a, self.record_structurals_b, self.record_structurals_c]

    def test_get_cttv_variant_type_singles(self):
        for record in self.test_records_singles:
            self.assertEqual(clinvar_to_evidence_strings.get_cttv_variant_type(record[0]), record[1])

    def test_get_cttv_variant_type_structurals(self):
        for record in self.test_records_structurals:
            self.assertEqual(clinvar_to_evidence_strings.get_cttv_variant_type(record[0]), record[1])
