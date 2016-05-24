import os
import unittest

from eva_cttv_pipeline import consequence_type as CT


class ProcessGeneTest(unittest.TestCase):
    def test__process_gene(self):
        test_consequence_type_dict = {}
        test_rs_id = "rs180177129"
        test_ensembl_gene_id = "ENSG00000083093"
        test_so_name = "intron_variant"

        test_consequence_type = CT.ConsequenceType(ensembl_gene_ids=[test_ensembl_gene_id], so_names=[test_so_name])

        CT._process_gene(test_consequence_type_dict, test_rs_id, test_ensembl_gene_id, test_so_name)

        self.assertEqual(test_consequence_type_dict["rs180177129"], test_consequence_type)


class ProcessConsequenceTypeFileTsvTest(unittest.TestCase):
    def test__process_consequence_type_file_tsv(self):
        test_consequence_type = CT.ConsequenceType(ensembl_gene_ids=["ENSG00000083093"], so_names=["intron_variant"])
        snp_2_gene_xls = os.path.join(os.path.dirname(__file__), 'resources',
                                      'cttv012_snp2gene_20160222_test_extract.tsv')
        consequence_type_dict, one_rs_multiple_genes = CT._process_consequence_type_file_tsv(snp_2_gene_xls)
        self.assertEqual(consequence_type_dict["rs180177129"], test_consequence_type)


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


class ConsequenceTypeTest(unittest.TestCase):
    def setUp(self):
        self.test_ensembl_gene_ids_a = {"ENSG00000083093"}
        self.test_so_name_a = "intron_variant"
        self.test_so_name_b = "transcript_ablation"
        self.test_so_term_a = CT.SoTerm(self.test_so_name_a)
        self.test_so_term_b = CT.SoTerm(self.test_so_name_b)
        self.test_consequence_type_a = CT.ConsequenceType(ensembl_gene_ids=self.test_ensembl_gene_ids_a,
                                                          so_names={self.test_so_name_a})
        self.test_consequence_type_b = CT.ConsequenceType()
        self.test_consequence_type_c = CT.ConsequenceType(ensembl_gene_ids=self.test_ensembl_gene_ids_a,
                                                          so_names={self.test_so_name_a, self.test_so_name_b})

    def test_ensembl_gene_ids(self):
        self.assertEqual(self.test_consequence_type_a.ensembl_gene_ids, self.test_ensembl_gene_ids_a)
        new_id = "new_id"
        self.test_ensembl_gene_ids_a.add(new_id)
        self.test_consequence_type_a.ensembl_gene_ids.add(new_id)
        self.assertEqual(self.test_consequence_type_a.ensembl_gene_ids, self.test_ensembl_gene_ids_a)

        new_ids = {"id1", "id2"}
        self.test_consequence_type_a.ensembl_gene_ids = new_ids
        self.assertEqual(self.test_consequence_type_a.ensembl_gene_ids, new_ids)

        self.assertEqual(self.test_consequence_type_b.ensembl_gene_ids, set())
        self.test_consequence_type_b.ensembl_gene_ids = new_ids
        self.assertEqual(self.test_consequence_type_b.ensembl_gene_ids, new_ids)

    def test_add_so_term(self):
        self.assertEqual(self.test_consequence_type_a.so_terms, {self.test_so_term_a})
        self.test_consequence_type_a.add_so_term(self.test_so_name_b)
        self.assertEqual(self.test_consequence_type_a.so_terms, {self.test_so_term_a, self.test_so_term_b})

        self.assertEqual(self.test_consequence_type_b.so_terms, set())
        self.test_consequence_type_b.add_so_term(self.test_so_name_b)
        self.assertEqual(self.test_consequence_type_b.so_terms, {self.test_so_term_b})

    def test_most_severe_so(self):
        self.assertEqual(self.test_consequence_type_a.most_severe_so, self.test_so_term_a)
        self.test_consequence_type_a.add_so_term(self.test_so_name_b)
        self.assertEqual(self.test_consequence_type_a.most_severe_so, self.test_so_term_b)

        self.assertEqual(self.test_consequence_type_c.most_severe_so, self.test_so_term_b)
