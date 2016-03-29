import unittest

import eva_cttv_pipeline.utilities as util


class ArgParserTest(unittest.TestCase):
    clin_sig = 'pathogenic,likely pathogenic'
    ignore = '/path/to/ignore/file'
    out = '/path/to/out/file'
    efo_map_file = '/path/to/efo/file'
    snp_2_gene_file = '/path/to/snp/to/gene/file'
    variant_summary_file = '/path/to/variant/summary/file'

    def setUp(self):
        argv = ['clinvar_to_evidence_strings.py', '--clinSig', self.clin_sig, '--ignore', self.ignore,
                '--out', self.out, '-e', self.efo_map_file, '-g', self.snp_2_gene_file, '-v', self.variant_summary_file]
        self.argparser = util.ArgParser(argv)

    def test_clin_sig(self):
        self.assertEquals(self.argparser.clinical_significance, self.clin_sig)

    def test_ignore(self):
        self.assertEquals(self.argparser.ignore_terms_file, self.ignore)

    def test_out(self):
        self.assertEquals(self.argparser.out, self.out)

    def test_efo_map_file(self):
        self.assertEquals(self.argparser.efo_mapping_file, self.efo_map_file)

    def test_snp_2_gene_file(self):
        self.assertEquals(self.argparser.snp_2_gene_file, self.snp_2_gene_file)

    def test_variant_summary_file(self):
        self.assertEquals(self.argparser.variant_summary_file, self.variant_summary_file)


class GetResourceFileTest(unittest.TestCase):
    def test_get_resource_file_existant(self):
        pass

    def test_get_resource_file_nonexistant(self):
        self.assertEqual(util.get_resource_file("not_a_real_package_39146", "not_a_real_file"), None)
