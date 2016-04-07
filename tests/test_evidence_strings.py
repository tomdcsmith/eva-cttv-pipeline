from datetime import datetime
import unittest

import eva_cttv_pipeline.evidence_strings as ES


class CTTVGeneticsEvidenceStringTest(unittest.TestCase):
    def setUp(self):
        self.test_ges = ES.CTTVGeneticsEvidenceString()

    def test_db_xref_url(self):
        url = "http://identifiers.org/clinvar.record/RCV000128628"
        self.test_ges.db_xref_url = url
        self.assertEqual(self.test_ges['evidence']['gene2variant']['provenance_type']['database']['dbxref']['url'], url)
        self.assertEqual(self.test_ges['evidence']['variant2disease']['provenance_type']['database']['dbxref']['url'], url)
        self.assertEqual(self.test_ges.db_xref_url, url)

    def test_url(self):
        url = "http://www.ncbi.nlm.nih.gov/clinvar/RCV000128628"
        self.test_ges.url = url
        self.assertEqual(self.test_ges['evidence']['gene2variant']['urls'][0]['url'], url)
        self.assertEqual(self.test_ges['evidence']['variant2disease']['urls'][0]['url'], url)
        self.assertEqual(self.test_ges.url, url)

    def test_gene_2_var_ev_codes(self):
        ev_codes = ['http://identifiers.org/eco/cttv_mapping_pipeline']
        self.test_ges.gene_2_var_ev_codes = ev_codes
        self.assertEqual(self.test_ges['evidence']['gene2variant']['evidence_codes'], ev_codes)
        self.assertEqual(self.test_ges.gene_2_var_ev_codes, ev_codes)

    def test_gene_2_var_func_consequence(self):
        functional_consequence = 'http://purl.obolibrary.org/obo/SO_0001583'
        self.test_ges.gene_2_var_func_consequence = functional_consequence
        self.assertEqual(self.test_ges['evidence']['gene2variant']['functional_consequence'], functional_consequence)
        self.assertEqual(self.test_ges.gene_2_var_func_consequence, functional_consequence)

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

    def test_association(self):
        self.test_ges.association = True
        self.assertTrue(self.test_ges['evidence']['gene2variant']['is_associated'])
        self.assertTrue(self.test_ges['evidence']['variant2disease']['is_associated'])
        self.assertTrue(self.test_ges.association)

        self.test_ges.association = False
        self.assertFalse(self.test_ges['evidence']['gene2variant']['is_associated'])
        self.assertFalse(self.test_ges['evidence']['variant2disease']['is_associated'])
        self.assertFalse(self.test_ges.association)

    def test_set_variant(self):
        test_id = "http://identifiers.org/dbsnp/rs193922494"
        test_type = "snp single"
        self.test_ges.set_variant(test_id, test_type)
        self.assertEqual(self.test_ges['variant']['id'], [test_id])
        self.assertEqual(self.test_ges['variant']['type'], test_type)

    def test_unique_reference(self):
        unique_reference = "http://europepmc.org/abstract/MED/0"
        self.test_ges.unique_reference = unique_reference
        self.assertEqual(self.test_ges['evidence']['variant2disease']['unique_experiment_reference'], unique_reference)
        self.assertEqual(self.test_ges.unique_reference, unique_reference)

    def test_date(self):
        date_string = datetime.fromtimestamp(1412982000000 / 1000).isoformat()
        self.test_ges.date = date_string
        self.assertEqual(self.test_ges['evidence']['gene2variant']['date_asserted'], date_string)
        self.assertEqual(self.test_ges['evidence']['variant2disease']['date_asserted'], date_string)
        self.assertEqual(self.test_ges.date, date_string)



class CTTVSomaticEvidenceStringTest(unittest.TestCase):
    def setUp(self):
        self.test_ses = ES.CTTVSomaticEvidenceString()

    def test_db_xref_url(self):
        url = "http://identifiers.org/clinvar.record/RCV000128628"
        self.test_ses.db_xref_url = url
        self.assertEqual(self.test_ses['evidence']['provenance_type']['database']['dbxref']['url'], url)
        self.assertEqual(self.test_ses.db_xref_url, url)

    def test_url(self):
        url = "http://www.ncbi.nlm.nih.gov/clinvar/RCV000128628"
        self.test_ses.url = url
        self.assertEqual(self.test_ses['evidence']['urls'][0]['url'], url)
        self.assertEqual(self.test_ses.url, url)

    def test_evidence_literature(self):
        literature_1 = "PMCID12345"
        self.test_ses.evidence_literature = [literature_1]
        self.assertEqual(self.test_ses['evidence']['provenance_type']['literature']['references'], [{"lit_id": literature_1}])
        self.assertEqual(self.test_ses.evidence_literature, [{"lit_id": literature_1}])

        literature_2 = "PMCID9876"
        literature_3 = "PMCID7654"
        literature_list = [literature_2, literature_3]
        self.test_ses.evidence_literature = literature_list
        self.assertEqual(self.test_ses['evidence']['provenance_type']['literature']['references'], [{"lit_id": literature_id} for literature_id in literature_list])
        self.assertEqual(self.test_ses.evidence_literature, [{"lit_id": literature_id} for literature_id in literature_list])

    def test_association(self):
        self.test_ses.association = True
        self.assertTrue(self.test_ses['evidence']['is_associated'])
        self.assertTrue(self.test_ses.association)

        self.test_ses.association = False
        self.assertFalse(self.test_ses['evidence']['is_associated'])
        self.assertFalse(self.test_ses.association)

    def test_date(self):
        date_string = datetime.fromtimestamp(1412982000000 / 1000).isoformat()
        self.test_ses.date = date_string
        self.assertEqual(self.test_ses['evidence']['date_asserted'], date_string)
        self.assertEqual(self.test_ses.date, date_string)


