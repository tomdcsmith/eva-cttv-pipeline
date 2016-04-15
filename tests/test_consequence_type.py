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

    def test_get_ranked_so_names(self):
        self.assertEqual(CT.SoTerm.get_ranked_so_names(), ['transcript_ablation', 'splice_acceptor_variant', 'splice_donor_variant', 'stop_gained', 'frameshift_variant', 'stop_lost', 'initiator_codon_variant', 'transcript_amplification', 'inframe_insertion', 'inframe_deletion', 'missense_variant', 'splice_region_variant', 'incomplete_terminal_codon_variant', 'stop_retained_variant', 'synonymous_variant', 'coding_sequence_variant', 'mature_miRNA_variant', '5_prime_UTR_variant', '3_prime_UTR_variant', 'non_coding_transcript_exon_variant', 'intron_variant', 'NMD_transcript_variant', 'non_coding_transcript_variant', 'upstream_gene_variant', 'downstream_gene_variant', 'TFBS_ablation', 'TFBS_amplification', 'TF_binding_site_variant', 'regulatory_region_ablation', 'regulatory_region_amplification', 'regulatory_region_variant', 'feature_elongation', 'feature_truncation', 'intergenic_variant'])


# class ConsequenceTypeTest(unittest.TestCase):
#
