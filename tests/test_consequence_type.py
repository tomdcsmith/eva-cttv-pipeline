from collections import defaultdict
import os
import unittest

from eva_cttv_pipeline import consequence_type as CT


class ProcessGeneTest(unittest.TestCase):
    def test__process_gene(self):
        test_consequence_type_dict = defaultdict(list)
        test_rs_id = "rs121912888"
        test_ensembl_gene_id = "ENSG00000139219"
        test_so_name = "missense_variant"

        test_consequence_type = CT.ConsequenceType(test_ensembl_gene_id, test_so_name)

        CT.process_gene(test_consequence_type_dict, test_rs_id, test_ensembl_gene_id, test_so_name)

        self.assertEqual(test_consequence_type_dict["rs121912888"][0], test_consequence_type)


class ProcessConsequenceTypeFileTsvTest(unittest.TestCase):
    def test__process_consequence_type_file_tsv(self):
        test_consequence_type = CT.ConsequenceType("ENSG00000021488", "missense_variant")
        snp_2_gene_file_path = os.path.join(os.path.dirname(__file__), 'resources',
                                      'coords_20170117_out_extract.tsv')
        consequence_type_dict, one_rs_multiple_genes = \
            CT.process_consequence_type_file_tsv(snp_2_gene_file_path)
        self.assertEqual(consequence_type_dict["rs121908485"][0], test_consequence_type)


class SoTermTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.test_so_term_a = CT.SoTerm("stop_gained")
        cls.test_so_term_b = CT.SoTerm("not_real_term")

    def test_accession(self):
        self.assertEqual(self.test_so_term_a.accession, "SO:0001587")
        self.assertIsNone(self.test_so_term_b.accession)

    def test_rank(self):
        self.assertEqual(self.test_so_term_a.rank, 3)
        self.assertEqual(self.test_so_term_b.rank, 34)
