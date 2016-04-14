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
        snp_2_gene_xls = os.path.join(os.path.dirname(__file__), 'resources', 'cttv012_snp2gene_20160222_test_extract.tsv')
        consequence_type_dict, one_rs_multiple_genes = CT._process_consequence_type_file_tsv(snp_2_gene_xls)
        self.assertEqual(consequence_type_dict["rs180177129"], test_consequence_type)

#
# class SoTermTest(unittest.TestCase):
#


