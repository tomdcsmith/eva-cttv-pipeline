from datetime import datetime
import unittest

import eva_cttv_pipeline.evidence_strings as ES


class CTTVGeneticsEvidenceStringTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.test_ges = ES.CTTVGeneticsEvidenceString()

    def test_set_db_xref_url(self):
        url = "http://identifiers.org/clinvar.record/RCV000128628"
        self.test_ges.set_db_xref_url(url)
        self.assertEqual(self.test_ges['evidence']['gene2variant']['provenance_type']['database']['dbxref']['url'], url)
        self.assertEqual(self.test_ges['evidence']['variant2disease']['provenance_type']['database']['dbxref']['url'], url)

    def test_set_url(self):
        url = "http://www.ncbi.nlm.nih.gov/clinvar/RCV000128628"
        self.test_ges.set_url(url)
        self.assertEqual(self.test_ges['evidence']['gene2variant']['urls'][0]['url'], url)
        self.assertEqual(self.test_ges['evidence']['variant2disease']['urls'][0]['url'], url)

    def test_set_gene_2_var_ev_codes(self):
        ev_codes = ['http://identifiers.org/eco/cttv_mapping_pipeline']
        self.test_ges.set_gene_2_var_ev_codes(ev_codes)
        self.assertEqual(self.test_ges['evidence']['gene2variant']['evidence_codes'], ev_codes)

    def test_set_gene_2_var_func_consequence(self):
        functional_consequence = 'http://purl.obolibrary.org/obo/SO_0001583'
        self.test_ges.set_gene_2_var_func_consequence(functional_consequence)
        self.assertEqual(self.test_ges['evidence']['gene2variant']['functional_consequence'], functional_consequence)

    def test_set_var_2_disease_literature_a(self):
        self.test_ges['evidence']['variant2disease']['provenance_type']['literature'] = {}

        literature_1 = "PMCID12345"
        self.test_ges.set_var_2_disease_literature([literature_1])
        self.assertEqual(self.test_ges['evidence']['variant2disease']['provenance_type']['literature']['references'], [{"lit_id": literature_1}])

        literature_2 = "PMCID9876"
        literature_3 = "PMCID7654"
        literature_list = [literature_2, literature_3]
        self.test_ges.set_var_2_disease_literature(literature_list)
        self.assertEqual(self.test_ges['evidence']['variant2disease']['provenance_type']['literature']['references'], [{"lit_id": literature_id} for literature_id in literature_list])

    def test_set_var_2_disease_literature_b(self):
        literature_1 = "PMCID12345"
        self.test_ges.set_var_2_disease_literature([literature_1])
        self.assertEqual(self.test_ges['evidence']['variant2disease']['provenance_type']['literature']['references'], [{"lit_id": literature_1}])

        literature_2 = "PMCID9876"
        literature_3 = "PMCID7654"
        literature_list = [literature_2, literature_3]
        self.test_ges.set_var_2_disease_literature(literature_list)
        self.assertEqual(self.test_ges['evidence']['variant2disease']['provenance_type']['literature']['references'], [{"lit_id": literature_id} for literature_id in literature_list])

    def test_set_association(self):
        self.test_ges.set_association(True)
        self.assertTrue(self.test_ges['evidence']['gene2variant']['is_associated'])
        self.assertTrue(self.test_ges['evidence']['variant2disease']['is_associated'])

        self.test_ges.set_association(False)
        self.assertFalse(self.test_ges['evidence']['gene2variant']['is_associated'])
        self.assertFalse(self.test_ges['evidence']['variant2disease']['is_associated'])

    def test_set_variant(self):
        test_id = "http://identifiers.org/dbsnp/rs193922494"
        test_type = "snp single"
        self.test_ges.set_variant(test_id, test_type)
        self.assertEqual(self.test_ges['variant']['id'], [test_id])
        self.assertEqual(self.test_ges['variant']['type'], test_type)

    def test_set_unique_reference(self):
        date_string = datetime.fromtimestamp(1412982000000 / 1000).isoformat()
        self.test_ges.set_date(date_string)
        self.assertEqual(self.test_ges['evidence']['gene2variant']['date_asserted'], date_string)
        self.assertEqual(self.test_ges['evidence']['variant2disease']['date_asserted'], date_string)



# class CTTVSomaticEvidenceStringTest(unittest.TestCase):
#     def setUp(self):
#         self.test_ses = ES.CTTVSomaticEvidenceString()
#
#     def test_set_db_xref_url(self):
#         url = "http://identifiers.org/clinvar.record/RCV000128628"
#         self.test_ges.set_db_xref_url(url)
#         self.assertEqual(self.test_ges['evidence']['gene2variant']['provenance_type']['database']['dbxref']['url'], url)
#         self.assertEqual(self.test_ges['evidence']['variant2disease']['provenance_type']['database']['dbxref']['url'], url)
#
#     def test_set_url(self):
#         url = "http://www.ncbi.nlm.nih.gov/clinvar/RCV000128628"
#         self.test_ges.set_url(url)
#         self.assertEqual(self.test_ges['evidence']['gene2variant']['urls'][0]['url'], url)
#         self.assertEqual(self.test_ges['evidence']['variant2disease']['urls'][0]['url'], url)
#
#     def test_set_gene_2_var_ev_codes(self):
#         ev_codes = ['http://identifiers.org/eco/cttv_mapping_pipeline']
#         self.test_ges.set_gene_2_var_ev_codes(ev_codes)
#         self.assertEqual(self.test_ges['evidence']['gene2variant']['evidence_codes'], ev_codes)
#
#     def test_set_gene_2_var_func_consequence(self):
#         functional_consequence = 'http://purl.obolibrary.org/obo/SO_0001583'
#         self.test_ges.set_gene_2_var_func_consequence(functional_consequence)
#         self.assertEqual(self.test_ges['evidence']['gene2variant']['functional_consequence'], functional_consequence)
#
#     def test_set_var_2_disease_literature_a(self):
#         self.test_ges['evidence']['variant2disease']['provenance_type']['literature'] = {}
#
#         literature_1 = "PMCID12345"
#         self.test_ges.set_var_2_disease_literature([literature_1])
#         self.assertEqual(self.test_ges['evidence']['variant2disease']['provenance_type']['literature']['references'], [{"lit_id": literature_1}])
#
#         literature_2 = "PMCID9876"
#         literature_3 = "PMCID7654"
#         literature_list = [literature_2, literature_3]
#         self.test_ges.set_var_2_disease_literature(literature_list)
#         self.assertEqual(self.test_ges['evidence']['variant2disease']['provenance_type']['literature']['references'], [{"lit_id": literature_id} for literature_id in literature_list])
#
#     def test_set_var_2_disease_literature_b(self):
#         literature_1 = "PMCID12345"
#         self.test_ges.set_var_2_disease_literature([literature_1])
#         self.assertEqual(self.test_ges['evidence']['variant2disease']['provenance_type']['literature']['references'], [{"lit_id": literature_1}])
#
#         literature_2 = "PMCID9876"
#         literature_3 = "PMCID7654"
#         literature_list = [literature_2, literature_3]
#         self.test_ges.set_var_2_disease_literature(literature_list)
#         self.assertEqual(self.test_ges['evidence']['variant2disease']['provenance_type']['literature']['references'], [{"lit_id": literature_id} for literature_id in literature_list])
#
#     def test_set_association(self):
#         self.test_ges.set_association(True)
#         self.assertTrue(self.test_ges['evidence']['gene2variant']['is_associated'])
#         self.assertTrue(self.test_ges['evidence']['variant2disease']['is_associated'])
#
#         self.test_ges.set_association(False)
#         self.assertFalse(self.test_ges['evidence']['gene2variant']['is_associated'])
#         self.assertFalse(self.test_ges['evidence']['variant2disease']['is_associated'])
#
#     def test_set_variant(self):
#         test_id = "http://identifiers.org/dbsnp/rs193922494"
#         test_type = "snp single"
#         self.test_ges.set_variant(test_id, test_type)
#         self.assertEqual(self.test_ges['variant']['id'], [test_id])
#         self.assertEqual(self.test_ges['variant']['type'], test_type)
#
#     def test_set_unique_reference(self):
#         date_string = datetime.fromtimestamp(1412982000000 / 1000).isoformat()
#         self.test_ges.set_date(date_string)
#         self.assertEqual(self.test_ges['evidence']['gene2variant']['date_asserted'], date_string)
#         self.assertEqual(self.test_ges['evidence']['variant2disease']['date_asserted'], date_string)
