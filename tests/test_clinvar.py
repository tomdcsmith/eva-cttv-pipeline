from datetime import datetime
import json
import os
import unittest

from eva_cttv_pipeline import clinvar
from eva_cttv_pipeline import consequence_type as CT
from eva_cttv_pipeline import utilities

from tests import config


class TestClinvarRecord(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.test_clinvar_record = get_test_record()

    def test_date(self):
        self.assertEqual(self.test_clinvar_record.date,
                         datetime.fromtimestamp(1435359600000/1000).isoformat())

    def test_score(self):
        self.assertEqual(self.test_clinvar_record.score, 1)

    def test_acc(self):
        self.assertEqual(self.test_clinvar_record.accession, "RCV000002127")

    def test_traits(self):
        self.assertEqual(self.test_clinvar_record.traits, [['Usher syndrome, type 3',
                                                            'Usher Syndrome, Type III',
                                                            'USHER SYNDROME, TYPE IIIA',
                                                            'Usher syndrome, type 3A']])

    def test_trait_pubmed_refs(self):
        self.assertEqual(self.test_clinvar_record.trait_pubmed_refs, [[21697857]])

    def test_observed_pubmed_refs(self):
        self.assertEqual(self.test_clinvar_record.observed_pubmed_refs, [11524702, 12145752])

    def test_clinical_significance(self):
        self.assertEqual(self.test_clinvar_record.clinical_significance, "Pathogenic")

    def test_allele_origins(self):
        self.assertEqual(self.test_clinvar_record.allele_origins, ['germline'])


class TestClinvarRecordMeasure(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.test_crm = get_test_record().measures[0]
        cls.consequence_type_dict = CT.process_consequence_type_file(config.snp_2_gene_file)

    def test_hgvs(self):
        self.assertEqual(self.test_crm.hgvs,
                         ['NM_174878.2:c.528T>G',
                          'NM_001256819.1:c.*142T>G',
                          'NM_052995.2:c.300T>G',
                          'NM_001195794.1:c.567T>G',
                          'NG_009168.1:g.49893T>G',
                          'NC_000003.12:g.150928107A>C',
                          'NC_000003.11:g.150645894A>C',
                          'NR_046380.2:n.1009T>G',
                          'NR_046380.1:n.1010T>G',
                          'p.Tyr176X',
                          'NP_443721.1:p.Tyr100Ter',
                          'NP_777367.1:p.Tyr176Ter',
                          'NP_001182723.1:p.Tyr189Ter'])

    def test_rs(self):
        self.assertEqual(self.test_crm.rs_id, "rs121908140")

    def test_nsv(self):
        self.assertEqual(self.test_crm.nsv_id, "nsv123456")

    def test_variant_type(self):
        self.assertEqual(self.test_crm.variant_type, "single nucleotide variant")

    def test_measure_set_pubmed_refs(self):
        self.assertEqual(self.test_crm.pubmed_refs, [])


def get_test_record():
    test_clinvar_record_filepath = os.path.join(os.path.dirname(__file__), 'resources',
                                              'test_clinvar_record.json')
    with utilities.open_file(test_clinvar_record_filepath, "rt") as f:
        test_record_dict = json.load(f)
    test_record = clinvar.ClinvarRecord(test_record_dict)
    return test_record

