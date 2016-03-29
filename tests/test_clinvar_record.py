from datetime import datetime
import json
import unittest

from eva_cttv_pipeline import clinvar_record
from eva_cttv_pipeline import utilities


class TestClinvarRecord(unittest.TestCase):
    def setUp(self):
        self.clinvar_record = get_test_record()

    #TODO get second clinvar record to test other return point of get_gene_id
    def test_get_gene_id(self):
        self.assertEqual(self.clinvar_record.get_gene_id(), "NM_000548")

    def test_get_ensembl_id(self):
        self.assertEqual(self.clinvar_record.get_ensembl_id(), "ENSG00000008710")

    def test_get_date(self):
        self.assertEqual(self.clinvar_record.get_date(), datetime.fromtimestamp(1412982000000/1000).isoformat())

    def test_get_score(self):
        self.assertEqual(self.clinvar_record.get_score(), None)

    def test_get_acc(self):
        self.assertEqual(self.clinvar_record.get_acc(), "RCV000055062")

    def test_get_traits(self):
        self.assertEqual(self.clinvar_record.get_traits(), [['Tuberous sclerosis syndrome']])

    def test_get_trait_pubmed_refs(self):
        self.assertEqual(self.clinvar_record.get_trait_pubmed_refs(), [[20301399]])

    def test_get_observed_pubmed_refs(self):
        self.assertEqual(self.clinvar_record.get_observed_pubmed_refs(), [])

    def test_get_measure_set_pubmed_refs(self):
        self.assertEqual(self.clinvar_record.get_measure_set_pubmed_refs(), [])

    def test_get_hgvs(self):
        self.assertEqual(self.clinvar_record.get_hgvs(), ['NM_000548.3:c.*154dup', 'NM_001009944.2:c.*963dupC', 'NG_005895.1:g.44459dupG', 'NC_000016.10:g.2088764dupG', 'NC_000016.9:g.2138765dupG', 'p.(=)'])

    def test_get_clinical_significance(self):
        self.assertEqual(self.clinvar_record.get_clinical_significance(), "not provided")

    #TODO needs rcv_to_rs dict
    def test_get_rs(self):
        pass

    #TODO needs rcv_to_nsv dict
    def test_get_nsv(self):
        pass

    #TODO needs consequence_type_dict and rsv_to_rs dicts
    def test_get_main_consequence_types(self):
        pass

    def test_get_variant_type(self):
        self.assertEqual(self.clinvar_record.get_variant_type(), "Duplication")

    def test_get_allele_origins(self):
        self.assertEqual(self.clinvar_record.get_allele_origins(), ['germline'])


class TestGetRcvToRSNSVMapping(unittest.TestCase):
    variant_summary_file_path = utilities.get_resource_file("eva_cttv_pipeline", "resources/variant_summary_2015-05_testing.txt")
    rcv_to_rs, rcv_to_nsv = clinvar_record.get_rcv_to_rsnsv_mapping(variant_summary_file_path)

    def test_rcv_to_rs(self):
        self.assertEqual(self.rcv_to_rs["RCV000000012"], "rs397704705")
        self.assertEqual(self.rcv_to_rs["RCV000000381"], "rs137854556")
        self.assertEqual(self.rcv_to_rs["RCV000000204"], "rs121965059")

    def test_rcv_to_nsv(self):
        self.assertEqual(self.rcv_to_nsv["RCV000004182"], "nsv1067860")
        self.assertEqual(self.rcv_to_nsv["RCV000004183"], "nsv1067861")
        self.assertEqual(self.rcv_to_nsv["RCV000004554"], "nsv1067916")


def get_test_record():
    test_record = clinvar_record.ClinvarRecord(
        {"recordStatus": "current", "title": "NM_000548.3(TSC2):c.*154dup AND Tuberous sclerosis syndrome",
         "referenceClinVarAssertion": {
             "clinVarAccession": {"acc": "RCV000055062", "version": 1, "type": "RCV", "dateUpdated": 1412982000000},
             "recordStatus": "current",
             "clinicalSignificance": {"reviewStatus": "NOT_CLASSIFIED_BY_SUBMITTER", "description": "not provided"},
             "assertion": {"type": "VARIATION_TO_DISEASE"}, "observedIn": [{"sample": {"origin": "germline",
                                                                                       "species": {"value": "human",
                                                                                                   "taxonomyId": 9606},
                                                                                       "affectedStatus": "yes"},
                                                                            "method": [
                                                                                {"methodType": "LITERATURE_ONLY"}],
                                                                            "observedData": [{"attribute": {
                                                                                "integerValue": 1,
                                                                                "type": "VariantAlleles"},
                                                                                              "id": 3619513}]}],
             "measureSet": {"measure": [
                 {"name": [{"elementValue": {"value": "NM_000548.3(TSC2):c.*154dup", "type": "Preferred"}}],
                  "attributeSet": [{"attribute": {"value": "NM_000548.3:c.*154dup", "type": "HGVS, coding, RefSeq",
                                                  "change": "c.*154dup"}}, {
                                       "attribute": {"value": "NM_001009944.2:c.*963dupC",
                                                     "type": "HGVS, coding, RefSeq", "change": "c.*963dupC"}}, {
                                       "attribute": {"value": "NG_005895.1:g.44459dupG",
                                                     "type": "HGVS, genomic, RefSeqGene", "change": "g.44459dupG"}}, {
                                       "attribute": {"value": "NC_000016.10:g.2088764dupG", "integerValue": 38,
                                                     "type": "HGVS, genomic, top level", "change": "g.2088764dupG"}}, {
                                       "attribute": {"value": "NC_000016.9:g.2138765dupG", "integerValue": 37,
                                                     "type": "HGVS, genomic, top level, previous",
                                                     "change": "g.2138765dupG"}},
                                   {"attribute": {"value": "p.(=)", "type": "HGVS, non-validated"}},
                                   {"attribute": {"value": "Exon 41", "type": "Location"}},
                                   {"attribute": {"value": "3 prime UTR variant", "type": "MolecularConsequence"},
                                    "xref": [{"db": "Sequence Ontology", "id": "SO:0001624", "status": "CURRENT"},
                                             {"db": "RefSeq", "id": "NM_001009944.2:c.*962_*963insC",
                                              "status": "CURRENT"}]},
                                   {"attribute": {"value": "500B downstream variant", "type": "MolecularConsequence"},
                                    "xref": [{"db": "Sequence Ontology", "id": "SO:0001634", "status": "CURRENT"},
                                             {"db": "RefSeq", "id": "NM_000548.3:c.*154_*155insG",
                                              "status": "CURRENT"}]}], "cytogeneticLocation": ["16p13.3"],
                  "sequenceLocation": [
                      {"assembly": "GRCh37", "chr": "16", "accession": "NC_000016.9", "start": 2138765, "stop": 2138765,
                       "variantLength": 2, "referenceAllele": "G", "alternateAllele": "GG"},
                      {"assembly": "GRCh38", "chr": "16", "accession": "NC_000016.10", "start": 2088764,
                       "stop": 2088764, "variantLength": 2, "referenceAllele": "G", "alternateAllele": "GG"}],
                  "measureRelationship": [{"name": [{"elementValue": {
                      "value": "polycystic kidney disease 1 (autosomal dominant)", "type": "Preferred"}}],
                                           "symbol": [{"elementValue": {"value": "PKD1", "type": "Preferred"}}],
                                           "attributeSet": [{"attribute": {
                                               "value": "Sufficient evidence for dosage pathogenicity",
                                               "dateValue": 1329868800000, "type": "Haploinsufficiency"}, "citation": [{
                                                                                                                           "url": "http://www.ncbi.nlm.nih.gov/projects/dbvar/ISCA/isca_gene.cgi?sym=PKD1"}]},
                                                            {"attribute": {"value": "No evidence available",
                                                                           "dateValue": 1329868800000,
                                                                           "type": "Triplosensitivity"}, "citation": [{
                                                                                                                          "url": "http://www.ncbi.nlm.nih.gov/projects/dbvar/ISCA/isca_gene.cgi?sym=PKD1"}]}],
                                           "sequenceLocation": [
                                               {"assembly": "GRCh37", "chr": "16", "accession": "NC_000016.9",
                                                "start": 2138710, "stop": 2185898, "strand": "-"},
                                               {"assembly": "GRCh38", "chr": "16", "accession": "NC_000016.10",
                                                "start": 2088707, "stop": 2135897, "strand": "-"}],
                                           "type": "variant in gene",
                                           "xref": [{"db": "Gene", "id": "5310", "status": "CURRENT"},
                                                    {"db": "OMIM", "id": "601313", "type": "MIM",
                                                     "status": "CURRENT"}]}, {"name": [
                      {"elementValue": {"value": "tuberous sclerosis 2", "type": "Preferred"}}], "symbol": [
                      {"elementValue": {"value": "TSC2", "type": "Preferred"}}], "attributeSet": [{"attribute": {
                      "value": "Sufficient evidence for dosage pathogenicity", "dateValue": 1329868800000,
                      "type": "Haploinsufficiency"}, "citation": [
                      {"url": "http://www.ncbi.nlm.nih.gov/projects/dbvar/ISCA/isca_gene.cgi?sym=TSC2"}]}, {
                                                                                                      "attribute": {
                                                                                                          "value": "No evidence available",
                                                                                                          "dateValue": 1329868800000,
                                                                                                          "type": "Triplosensitivity"},
                                                                                                      "citation": [{
                                                                                                                       "url": "http://www.ncbi.nlm.nih.gov/projects/dbvar/ISCA/isca_gene.cgi?sym=TSC2"}]}],
                                                                              "sequenceLocation": [
                                                                                  {"assembly": "GRCh37", "chr": "16",
                                                                                   "accession": "NC_000016.9",
                                                                                   "start": 2097989, "stop": 2138712,
                                                                                   "strand": "+"},
                                                                                  {"assembly": "GRCh38", "chr": "16",
                                                                                   "accession": "NC_000016.10",
                                                                                   "start": 2047801, "stop": 2088711,
                                                                                   "strand": "+"}], "comment": [{
                                                                                                                    "value": "This gene is cited in the ACMG recommendations of 2013 (PubMed 23788249) for reporting incidental findings in exons.",
                                                                                                                    "dataSource": "NCBI curation"}],
                                                                              "type": "variant in gene", "xref": [
                          {"db": "Gene", "id": "7249", "status": "CURRENT"},
                          {"db": "OMIM", "id": "191092", "type": "MIM", "status": "CURRENT"}]}], "type": "Duplication",
                  "id": 75791,
                  "xref": [{"db": "Tuberous sclerosis database (TSC2)", "id": "TSC2_02318", "status": "CURRENT"},
                           {"db": "dbSNP", "id": "397514891", "type": "rs", "status": "CURRENT"}]}], "name": [
                 {"elementValue": {"value": "NM_000548.3(TSC2):c.*154dup", "type": "preferred name"}}],
                            "type": "Variant", "id": 64862}, "traitSet": {"trait": [{"name": [
                 {"elementValue": {"value": "Tuberous sclerosis syndrome", "type": "Preferred"},
                  "xref": [{"db": "SNOMED CT", "id": "7199000", "status": "CURRENT"}]}], "symbol": [
                 {"elementValue": {"value": "TSC", "type": "Preferred"},
                  "xref": [{"db": "OMIM", "id": "191100", "type": "MIM", "status": "CURRENT"}]},
                 {"elementValue": {"value": "TS", "type": "Alternate"},
                  "xref": [{"db": "OMIM", "id": "191100", "type": "MIM", "status": "CURRENT"}]}], "attributeSet": [{
                                                                                                                       "attribute": {
                                                                                                                           "value": "Tuberous sclerosis complex (TSC) involves abnormalities of the skin (hypomelanotic macules, facial angiofibromas, shagreen patches, fibrous facial plaques, ungual fibromas); brain (cortical tubers, subependymal nodules [SENs] and subependymal giant cell astrocytomas [SEGAs], seizures, intellectual disability/developmental delay); kidney (angiomyolipomas, cysts, renal cell carcinomas); heart (rhabdomyomas, arrhythmias); and lungs (lymphangioleiomyomatosis [LAM]). CNS tumors are the leading cause of morbidity and mortality; renal disease is the second leading cause of early death.",
                                                                                                                           "type": "public definition"},
                                                                                                                       "xref": [
                                                                                                                           {
                                                                                                                               "db": "GeneReviews",
                                                                                                                               "id": "NBK1220",
                                                                                                                               "status": "CURRENT"}]},
                                                                                                                   {
                                                                                                                       "attribute": {
                                                                                                                           "value": "Neoplasm",
                                                                                                                           "type": "keyword"}},
                                                                                                                   {
                                                                                                                       "attribute": {
                                                                                                                           "value": "Hereditary cancer syndrome",
                                                                                                                           "type": "keyword"}}],
                                                                                     "citation": [{
                                                                                                      "url": "https://www.orpha.net/data/patho/Pro/en/Emergency_TuberousSclerosis.pdf",
                                                                                                      "citationText": "Orphanet, Tuberous sclerosis, 2007",
                                                                                                      "type": "practice guideline",
                                                                                                      "abbrev": "Orphanet, 2007"},
                                                                                                  {"id": {
                                                                                                      "value": "20301399",
                                                                                                      "source": "PubMed"},
                                                                                                   "type": "review",
                                                                                                   "abbrev": "GeneReviews"}],
                                                                                     "type": "Disease", "id": 15993,
                                                                                     "xref": [{"db": "MedGen",
                                                                                               "id": "C0041341",
                                                                                               "status": "CURRENT"}]}],
                                                                          "type": "Disease", "id": 8139},
             "dateCreated": 1379286000000, "dateLastUpdated": 1412982000000, "id": 144533}, "clinVarAssertion": [{
                                                                                                                     "clinVarSubmissionID": {
                                                                                                                         "submitter": "Tuberous sclerosis database (TSC2)",
                                                                                                                         "title": "NM_000548.3:c.*154dup AND TSC",
                                                                                                                         "localKey": "TSC2_02318_TSC",
                                                                                                                         "submitterDate": 1376002800000},
                                                                                                                     "clinVarAccession": {
                                                                                                                         "acc": "SCV000083280",
                                                                                                                         "version": 1,
                                                                                                                         "type": "SCV",
                                                                                                                         "orgID": 500074,
                                                                                                                         "dateUpdated": 1403478000000},
                                                                                                                     "recordStatus": "current",
                                                                                                                     "clinicalSignificance": {
                                                                                                                         "reviewStatus": "CLASSIFIED_BY_SINGLE_SUBMITTER",
                                                                                                                         "description": [
                                                                                                                             "not provided"]},
                                                                                                                     "assertion": {
                                                                                                                         "type": "variation to disease"},
                                                                                                                     "externalID": {
                                                                                                                         "db": "Tuberous sclerosis database (TSC2)",
                                                                                                                         "id": "TSC2_02318",
                                                                                                                         "status": "CURRENT"},
                                                                                                                     "observedIn": [
                                                                                                                         {
                                                                                                                             "sample": {
                                                                                                                                 "origin": "germline",
                                                                                                                                 "species": {
                                                                                                                                     "value": "human"},
                                                                                                                                 "affectedStatus": "yes"},
                                                                                                                             "method": [
                                                                                                                                 {
                                                                                                                                     "methodType": "LITERATURE_ONLY"}],
                                                                                                                             "observedData": [
                                                                                                                                 {
                                                                                                                                     "attribute": {
                                                                                                                                         "value": "1",
                                                                                                                                         "type": "VariantAlleles"}}]}],
                                                                                                                     "measureSet": {
                                                                                                                         "measure": [
                                                                                                                             {
                                                                                                                                 "attributeSet": [
                                                                                                                                     {
                                                                                                                                         "attribute": {
                                                                                                                                             "value": "NM_000548.3:c.*154dup",
                                                                                                                                             "type": "HGVS"}},
                                                                                                                                     {
                                                                                                                                         "attribute": {
                                                                                                                                             "value": "p.(=)",
                                                                                                                                             "type": "HGVS"}},
                                                                                                                                     {
                                                                                                                                         "attribute": {
                                                                                                                                             "value": "Exon 41",
                                                                                                                                             "type": "Location"}}],
                                                                                                                                 "measureRelationship": [
                                                                                                                                     {
                                                                                                                                         "symbol": [
                                                                                                                                             {
                                                                                                                                                 "elementValue": {
                                                                                                                                                     "value": "TSC2",
                                                                                                                                                     "type": "Preferred"}}],
                                                                                                                                         "type": "variant in gene"}],
                                                                                                                                 "type": "Variation"}],
                                                                                                                         "type": "Variant"},
                                                                                                                     "traitSet": {
                                                                                                                         "trait": [
                                                                                                                             {
                                                                                                                                 "name": [
                                                                                                                                     {
                                                                                                                                         "elementValue": {
                                                                                                                                             "value": "TSC",
                                                                                                                                             "type": "Preferred"}}],
                                                                                                                                 "type": "Disease"}],
                                                                                                                         "type": "Disease"},
                                                                                                                     "id": 143786}],
         "id": 3756609})
    # record_string = json.load(test_record)
    return test_record


print(get_test_record().get_allele_origins())