import unittest

import eva_cttv_pipeline.efo_term as EFOT


# TODO test get_available_terms, execute_query, get_obsolete_terms


class EFOTermTest(unittest.TestCase):
    # def setUpClass(cls):
    #     cls.obsolete_terms = ["ob_term_1", "ob_term_2"]
    #     cls.cttv_available_terms = ["cttv_av_term_1", "cttv_av_term_2"]

    def setUp(self):
        self.test_efoterm_id = EFOT.EFOTerm("test_id")

        self.test_efoterm_ob_a = EFOT.EFOTerm("ob_term_1")
        self.test_efoterm_ob_b = EFOT.EFOTerm("not_ob")

        self.test_efoterm_cttv_a = EFOT.EFOTerm("cttv_av_term_1")
        self.test_efoterm_cttv_b = EFOT.EFOTerm("not_cttv_av")

    def test_efoid(self):
        self.assertEqual(self.test_efoterm_id.efoid, "test_id")
        self.test_efoterm_id.efoid = "new_id"
        self.assertNotEquals(self.test_efoterm_id.efoid, "test_id")
        self.assertEqual(self.test_efoterm_id.efoid, "new_id")

    # def test_is_obsolete(self):
    #     self.assertRaises(EFOT.EFOTerm.IsObsoleteException, self.test_efoterm_ob_a.is_obsolete())
