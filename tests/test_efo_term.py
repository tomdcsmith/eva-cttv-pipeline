import unittest

import eva_cttv_pipeline.efo_term as EFOT


# TODO test get_available_terms, execute_query, get_obsolete_terms


class EFOTermTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        EFOT.EFOTerm.obsolete_terms = {"ob_id_1": "ob_term_1", "ob_id_2": "ob_term_2"}
        EFOT.EFOTerm.cttv_available_terms = \
            {"cttv_id_1": "cttv_av_term_1", "cttv_id_2": "cttv_av_term_2"}

        cls.test_efoterm_id = EFOT.EFOTerm("test_id")

        cls.test_efoterm_ob_a = EFOT.EFOTerm("ob_id_1")
        cls.test_efoterm_ob_b = EFOT.EFOTerm("not_ob")

        cls.test_efoterm_cttv_a = EFOT.EFOTerm("cttv_id_1")
        cls.test_efoterm_cttv_b = EFOT.EFOTerm("not_cttv_av")

    def test_efoid(self):
        self.assertEqual(self.test_efoterm_id.efoid, "test_id")
        self.test_efoterm_id.efoid = "new_id"
        self.assertNotEquals(self.test_efoterm_id.efoid, "test_id")
        self.assertEqual(self.test_efoterm_id.efoid, "new_id")

    def test_is_obsolete(self):
        self.assertRaises(EFOT.EFOTerm.IsObsoleteException, self.test_efoterm_ob_a.is_obsolete)
        self.assertFalse(self.test_efoterm_ob_b.is_obsolete())

    def test_is_cttv_available(self):
        self.assertRaises(EFOT.EFOTerm.NotCttvAvailableException,
                          self.test_efoterm_cttv_b.is_cttv_available)
        self.assertTrue(self.test_efoterm_cttv_a.is_cttv_available())
