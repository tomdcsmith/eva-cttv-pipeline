import unittest

import eva_cttv_pipeline.evidence_strings as ES


class CTTVGeneticsEvidenceStringTest(unittest.TestCase):
    def setUp(self):
        self.gen_ev_string = ES.CTTVGeneticsEvidenceString()

    def test_set_db_xref_url(self):
        test_urls = ["http://purl.obolibrary.org/obo/GO_0044691", "http://purl.obolibrary.org/obo/HP_0000113", "http://www.orpha.net/ORDO/Orphanet_71517"]
        for url in test_urls:
            self.gen_ev_string.set_db_xref_url(url)
            ev_url1 = self.gen_ev_string['evidence']['gene2variant']['provenance_type']['database']['dbxref']['url']
            self.assertEqual(url, ev_url1)
            ev_url2 = self.gen_ev_string['evidence']['variant2disease']['provenance_type']['database']['dbxref']['url']
            self.assertEqual(url, ev_url2)
            self.assertNotEquals(url, "not a url")
