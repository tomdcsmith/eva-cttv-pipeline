import unittest

import eva_cttv_pipeline.trait_mapping.trait as trait
import eva_cttv_pipeline.trait_mapping.zooma as zooma


class TestTrait(unittest.TestCase):
    def test_is_finished_true(self):
        test_trait = trait.Trait('aprt deficiency, japanese type', 1)
        test_trait.finished_mapping_set.add(trait.OntologyEntry('http://www.orpha.net/ORDO/Orphanet_976',
                                                                'Adenine phosphoribosyltransferase deficiency'))
        self.assertTrue(test_trait.is_finished)

    def test_is_finished_false(self):
        test_trait = trait.Trait('aprt deficiency, japanese type', 1)
        self.assertFalse(test_trait.is_finished)

    def test_process_zooma_result(self):
        test_trait = trait.Trait('aprt deficiency, japanese type', 1)

        test_zooma_result = zooma.ZoomaResult(['http://www.orpha.net/ORDO/Orphanet_976'],
                                              'Adenine phosphoribosyltransferase deficiency',
                                              'HIGH', 'eva-clinvar')
        entry = test_zooma_result.mapping_list[0]
        entry.confidence = zooma.ZoomaConfidence.HIGH
        entry.in_efo = True
        entry.is_current = True
        entry.ontology_label = "Adenine phosphoribosyltransferase deficiency"
        entry.source = 'eva-clinvar'
        entry.uri = 'http://www.orpha.net/ORDO/Orphanet_976'

        test_trait.zooma_result_list.append(test_zooma_result)

        test_trait.process_zooma_results()
        self.assertEqual(1, len(test_trait.finished_mapping_set))




