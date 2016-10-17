import itertools
import os
import unittest

import sys

from eva_cttv_pipeline import clinvar_to_evidence_strings
from tests import test_clinvar


def _get_mappings():
    efo_mapping_file = os.path.join(os.path.dirname(__file__), 'resources',
                                    'feb16_jul16_combined_trait_to_url.tsv')
    ignore_file = os.path.join(os.path.dirname(__file__), 'resources', 'ignore_file.txt')
    snp_2_gene_file = os.path.join(os.path.dirname(__file__), 'resources',
                                   'snp2gene_assignment_jul2016_extract.tsv')
    variant_summary_file = os.path.join(os.path.dirname(__file__), 'resources',
                                        'variant_summary_2016-05_test_extract.txt')

    mappings = clinvar_to_evidence_strings.get_mappings(efo_mapping_file, ignore_file, None,
                                                        snp_2_gene_file, variant_summary_file)

    return mappings


MAPPINGS = _get_mappings()


class GetMappingsTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.mappings = MAPPINGS

    def test_efo_mapping(self):
        self.assertEqual(len(self.mappings.trait_2_efo), 5055)
        self.assertEqual(len(self.mappings.unavailable_efo_dict), 0)

        self.assertEqual(self.mappings.trait_2_efo["renal-hepatic-pancreatic dysplasia 2"],
                         ['http://www.orpha.net/ORDO/Orphanet_294415'])
        self.assertEqual(self.mappings.trait_2_efo["frontotemporal dementia"],
                         ['http://purl.obolibrary.org/obo/HP_0000713'])
        self.assertEqual(
            self.mappings.trait_2_efo["3 beta-hydroxysteroid dehydrogenase deficiency"],
            ['http://www.orpha.net/ORDO/Orphanet_90791'])

    def test_consequence_type_dict(self):
        self.assertEqual(len(self.mappings.consequence_type_dict), 56)

        self.assertTrue("rs121908485" in self.mappings.consequence_type_dict)
        self.assertTrue("rs121912888" in self.mappings.consequence_type_dict)
        self.assertTrue("rs137852558" in self.mappings.consequence_type_dict)
        self.assertTrue("rs137853008" in self.mappings.consequence_type_dict)

        self.assertFalse("rs0" in self.mappings.consequence_type_dict)
        self.assertFalse("rs5" in self.mappings.consequence_type_dict)
        self.assertFalse("rs9" in self.mappings.consequence_type_dict)

    def test_rcv_to_rs_nsv(self):
        self.assertEqual(len(self.mappings.rcv_to_rs), 21)
        self.assertEqual(len(self.mappings.rcv_to_nsv), 6)

        self.assertEqual(self.mappings.rcv_to_nsv["RCV000020147"], "nsv1067916")
        self.assertEqual(self.mappings.rcv_to_nsv["RCV000004182"], "nsv1067860")
        self.assertEqual(self.mappings.rcv_to_nsv["RCV000004183"], "nsv1067861")

        self.assertEqual(self.mappings.rcv_to_rs["RCV000000012"], "rs397704705")
        self.assertEqual(self.mappings.rcv_to_rs["RCV000000204"], "rs121965059")
        self.assertEqual(self.mappings.rcv_to_rs["RCV000000381"], "rs137854556")


class CreateTraitTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.trait = clinvar_to_evidence_strings.create_trait(9, ["Ciliary dyskinesia, primary, 7"],
                                                             MAPPINGS.trait_2_efo)

    def test_clinvar_trait_list(self):
        self.assertEqual(self.trait.clinvar_trait_list, ['Ciliary dyskinesia, primary, 7'])

    def test_efo_list(self):
        self.assertEqual(self.trait.efo_list, ['http://www.ebi.ac.uk/efo/EFO_0003900'])

    def test_return_none(self):
        none_trait = \
            clinvar_to_evidence_strings.create_trait(9, ["not a real trait"], MAPPINGS.trait_2_efo)
        self.assertIsNone(none_trait)


class SkipRecordTest(unittest.TestCase):

    def setUp(self):
        self.clinvar_record = test_clinvar.get_test_record()
        report = clinvar_to_evidence_strings.Report()
        # skip_record(clinvarRecord, cellbase_record, allowed_clinical_significance, counters)
        self.args = [self.clinvar_record, {"reference": "A", "alternate": "T"},
                     ["not provided"], report.counters]
        # allowed clin sig changed to just "non provided" to match that in the test record

    def test_return_true(self):
        self.assertTrue(clinvar_to_evidence_strings.skip_record(*self.args))

    def test_ref_eq_alt(self):
        self.args[1] = {"reference": "A", "alternate": "A"}
        self.assertTrue(clinvar_to_evidence_strings.skip_record(*self.args))

    def test_rs_is_none(self):
        self.clinvar_record.rs = None
        self.assertTrue(clinvar_to_evidence_strings.skip_record(*self.args))

    def test_con_type_is_none(self):
        self.clinvar_record.consequence_type = None
        self.assertTrue(clinvar_to_evidence_strings.skip_record(*self.args))


class LoadEfoMappingTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        ignore_file = os.path.join(os.path.dirname(__file__), 'resources', 'ignore_file.txt')
        efo_file = \
            os.path.join(os.path.dirname(__file__), 'resources', 'feb16_jul16_combined_trait_to_url.tsv')

        cls.trait_2_efo, cls.unavailable_efo = \
            clinvar_to_evidence_strings.load_efo_mapping(efo_file)
        cls.trait_2_efo_w_ignore, cls.unavailable_efo_w_ignore = \
            clinvar_to_evidence_strings.load_efo_mapping(efo_file, ignore_terms_file=ignore_file)

    def test_just_mapping_trait_2_efo(self):
        self.assertEqual(len(self.trait_2_efo), 5283)

    def test_w_ignore_trait_2_efo(self):
        self.assertEqual(len(self.trait_2_efo_w_ignore), 5055)


class GetUnmappedUrlTest(unittest.TestCase):
    def test_orphanet(self):
        url = "http://www.orpha.net/ORDO/Orphanet_2670"
        self.assertEqual(clinvar_to_evidence_strings.get_unmapped_url(url),
                         "http://purl.bioontology.org/ORDO/Orphanet_2670")

    def test_hp(self):
        url = "http://purl.obolibrary.org/obo/HP_0000545"
        self.assertEqual(clinvar_to_evidence_strings.get_unmapped_url(url),
                         "http://purl.bioontology.org/obo/HP_0000545")

    def test_bad_url(self):
        #TODO currently the function exits execution with bad urls
        pass


class GetTermsFromFileTest(unittest.TestCase):
    #TODO do the same for adapt terms file?
    @classmethod
    def setUpClass(cls):
        ignore_file = os.path.join(os.path.dirname(__file__), 'resources', 'ignore_file.txt')
        cls.ignore_terms = clinvar_to_evidence_strings.get_terms_from_file(ignore_file)

    def test_with_file(self):
        self.assertEqual(len(self.ignore_terms), 218)
        self.assertEqual(self.ignore_terms[0], "http://purl.obolibrary.org/obo/HP_0011677")
        self.assertEqual(self.ignore_terms[-1], "http://www.orpha.net/ORDO/Orphanet_120795")

    def test_no_file(self):
        self.assertEqual(clinvar_to_evidence_strings.get_terms_from_file(None), [])


# class Clinvar2OTAlleleOrigin(unittest.TestCase):
#     @classmethod
#     def setUpClass(cls):
#         cls.germline_allele_origins = ["germline", "inherited", "maternal", "paternal",
#                                        "biparental", "uniparental"]
#         cls.somatic_allele_origins = ["somatic", "de novo"]
#         cls.unused_allele_origins = ["unknown", "not provided", "tested-inconclusive",
#                                      "not applicable"]
#
#         cls.germline_combinations = []
#         cls.somatic_combinations = []
#         cls.unused_combinations = []
#         cls.populated_germline_combinations = []
#         cls.populated_somatic_combinations = []
#         cls.populated_unused_combinations = []
#
#         for L in range(0, len(cls.germline_allele_origins) + 1):
#             combinations = list(itertools.combinations(cls.germline_allele_origins, L))
#             cls.germline_combinations += combinations
#             if len(combinations) > 1:
#                 cls.populated_germline_combinations += combinations
#         #
#         # print(cls.germline_combinations)
#         # print()
#         # print(cls.populated_germline_combinations)
#         # sys.exit(1)
#
#         for L in range(0, len(cls.somatic_allele_origins) + 1):
#             combinations = list(itertools.combinations(cls.somatic_allele_origins, L))
#             cls.somatic_combinations += combinations
#             if len(combinations) > 1:
#                 cls.populated_somatic_combinations += combinations
#
#         for L in range(0, len(cls.unused_allele_origins) + 1):
#             combinations = list(itertools.combinations(cls.unused_allele_origins, L))
#             cls.unused_combinations += combinations
#             if len(combinations) > 1:
#                 cls.populated_unused_combinations += combinations
#
#     def test_just_germline(self):
#         output_test_list = ["germline"]
#
#         input_test_lists = []
#         for germ in self.populated_germline_combinations:
#             for unu in self.unused_combinations:
#                 input_test_lists.append([list(germ) + list(unu)])
#
#         # input_test_lists = [germ + unu
#         #                     for germ in self.populated_germline_combinations
#         #                     for unu in self.unused_combinations]
#
#         for test_list in input_test_lists:
#             self.assertEqual(clinvar_to_evidence_strings.clinvar_2_ot_allele_origin(test_list), output_test_list)
#
#     def test_just_somatic(self):
#         output_test_list = ["somatic"]
#
#         input_test_lists = []
#         for som in self.populated_somatic_combinations:
#             for unu in self.unused_combinations:
#                 input_test_lists.append([list(som) + list(unu)])
#
#         # input_test_lists = [som + unu
#         #                     for som in self.populated_somatic_combinations
#         #                     for unu in self.unused_combinations]
#
#         for test_list in input_test_lists:
#             self.assertEqual(clinvar_to_evidence_strings.clinvar_2_ot_allele_origin(test_list), output_test_list)
#
#     def test_just_unused(self):
#         output_test_list = []
#
#         input_test_lists = self.unused_combinations
#
#         for test_list in input_test_lists:
#             self.assertEqual(clinvar_to_evidence_strings.clinvar_2_ot_allele_origin(test_list), output_test_list)
#
#     def test_germline_and_somatic(self):
#         output_test_list = ["somatic", "germline"]
#
#         input_test_lists = []
#         for som in self.populated_somatic_combinations:
#             for germ in self.populated_germline_combinations:
#                 input_test_lists.append([list(som) + list(germ)])
#
#         # print(input_test_lists)
#
#         # input_test_lists = [som + germ
#         #                     for som in self.populated_somatic_combinations
#         #                     for germ in self.populated_germline_combinations]
#
#         for test_list in input_test_lists:
#             self.assertEqual(sorted(clinvar_to_evidence_strings.clinvar_2_ot_allele_origin(test_list)), sorted(output_test_list))


class TestGetDefaultAllowedClincalSignificance(unittest.TestCase):
    def test_get_default_allowed_clincal_significance(self):
        correct_list = ['unknown', 'untested', 'non-pathogenic', 'probable-non-pathogenic',
         'probable-pathogenic', 'pathogenic', 'drug-response', 'drug response',
         'histocompatibility', 'other', 'benign', 'protective', 'not provided',
         'likely benign', 'confers sensitivity', 'uncertain significance',
         'likely pathogenic', 'conflicting data from submitters', 'risk factor',
         'association']
        self.assertEqual(clinvar_to_evidence_strings.get_default_allowed_clinical_significance(),
                         correct_list)
