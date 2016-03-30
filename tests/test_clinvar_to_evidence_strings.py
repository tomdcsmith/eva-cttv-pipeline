import json
import unittest

from eva_cttv_pipeline import clinvar_to_evidence_strings
from eva_cttv_pipeline import utilities
from eva_cttv_pipeline import consequence_type
from eva_cttv_pipeline import clinvar_record


def get_args_GetCttvGeneticsEvidenceStringTest():
    efo_list = ['http://www.orpha.net/ORDO/Orphanet_88991']
    clin_sig = "likely pathogenic"
    clin_sig_2_activity = {'association': 'http://identifiers.org/cttv.activity/damaging_to_target',
                           'likely pathogenic': 'http://identifiers.org/cttv.activity/predicted_damaging',
                           'untested': 'http://identifiers.org/cttv.activity/unknown',
                           'histocompatibility': 'http://identifiers.org/cttv.activity/unknown',
                           'pathogenic': 'http://identifiers.org/cttv.activity/damaging_to_target',
                           'not provided': 'http://identifiers.org/cttv.activity/unknown',
                           'non-pathogenic': 'http://identifiers.org/cttv.activity/tolerated_by_target',
                           'uncertain significance': 'http://identifiers.org/cttv.activity/unknown',
                           'likely benign': 'http://identifiers.org/cttv.activity/predicted_tolerated',
                           'protective': 'http://identifiers.org/cttv.activity/tolerated_by_target',
                           'benign': 'http://identifiers.org/cttv.activity/tolerated_by_target',
                           'risk factor': 'http://identifiers.org/cttv.activity/predicted_damaging',
                           'probable-non-pathogenic': 'http://identifiers.org/cttv.activity/predicted_tolerated',
                           'confers sensitivity': 'http://identifiers.org/cttv.activity/predicted_damaging',
                           'unknown': 'http://identifiers.org/cttv.activity/unknown',
                           'other': 'http://identifiers.org/cttv.activity/unknown',
                           'drug-response': 'http://identifiers.org/cttv.activity/unknown',
                           'conflicting data from submitters': 'http://identifiers.org/cttv.activity/unknown',
                           'probable-pathogenic': 'http://identifiers.org/cttv.activity/predicted_damaging'}
    clinvarRecord = clinvar_record.ClinvarRecord({'referenceClinVarAssertion': {'dateLastUpdated': 1412982000000, 'dateCreated': 1405551600000, 'id': 300575, 'traitSet': {'trait': [{'xref': [{'id': 'C0018798', 'db': 'MedGen', 'status': 'CURRENT'}, {'id': '140500', 'db': 'OMIM', 'status': 'CURRENT', 'type': 'MIM'}, {'id': '234750', 'db': 'OMIM', 'status': 'CURRENT', 'type': 'MIM'}], 'name': [{'elementValue': {'value': 'Malformation of the heart', 'type': 'Preferred'}}, {'xref': [{'id': '140500', 'db': 'OMIM', 'status': 'CURRENT', 'type': 'MIM'}, {'id': '234750', 'db': 'OMIM', 'status': 'CURRENT', 'type': 'MIM'}], 'elementValue': {'value': 'HEART, MALFORMATION OF', 'type': 'Alternate'}}, {'elementValue': {'value': 'Congenital heart defect', 'type': 'Alternate'}}], 'id': 15882, 'type': 'Disease'}], 'id': 16428, 'type': 'Disease'}, 'measureSet': {'name': [{'elementValue': {'value': 'NM_002471.3(MYH6):c.2033A>G (p.Asn678Ser)', 'type': 'preferred name'}}], 'id': 139663, 'measure': [{'id': 150127, 'cytogeneticLocation': ['14q11.2'], 'measureRelationship': [{'symbol': [{'elementValue': {'value': 'MYH6', 'type': 'Preferred'}}], 'xref': [{'id': '4624', 'db': 'Gene', 'status': 'CURRENT'}, {'id': '160710', 'db': 'OMIM', 'status': 'CURRENT', 'type': 'MIM'}], 'type': 'variant in gene', 'name': [{'elementValue': {'value': 'myosin, heavy chain 6, cardiac muscle, alpha', 'type': 'Preferred'}}], 'sequenceLocation': [{'stop': 23877485, 'accession': 'NC_000014.8', 'assembly': 'GRCh37', 'chr': '14', 'strand': '-', 'start': 23851198}, {'stop': 23409621, 'accession': 'NC_000014.9', 'assembly': 'GRCh38', 'chr': '14', 'strand': '-', 'start': 23380732}]}], 'sequenceLocation': [{'stop': 23866396, 'accession': 'NC_000014.8', 'variantLength': 1, 'assembly': 'GRCh37', 'referenceAllele': 'T', 'chr': '14', 'start': 23866396, 'alternateAllele': 'C'}, {'stop': 23397187, 'accession': 'NC_000014.9', 'variantLength': 1, 'assembly': 'GRCh38', 'referenceAllele': 'T', 'chr': '14', 'start': 23397187, 'alternateAllele': 'C'}], 'attributeSet': [{'attribute': {'change': 'c.2033A>G', 'value': 'LRG_389t1:c.2033A>G', 'type': 'HGVS, coding, LRG'}}, {'attribute': {'change': 'c.2033A>G', 'value': 'NM_002471.3:c.2033A>G', 'type': 'HGVS, coding, RefSeq'}}, {'attribute': {'change': 'g.16091A>G', 'value': 'LRG_389:g.16091A>G', 'type': 'HGVS, genomic, LRG'}}, {'attribute': {'change': 'g.16091A>G', 'value': 'NG_023444.1:g.16091A>G', 'type': 'HGVS, genomic, RefSeqGene'}}, {'attribute': {'change': 'g.23397187T>C', 'value': 'NC_000014.9:g.23397187T>C', 'integerValue': 38, 'type': 'HGVS, genomic, top level'}}, {'attribute': {'change': 'g.23866396T>C', 'value': 'NC_000014.8:g.23866396T>C', 'integerValue': 37, 'type': 'HGVS, genomic, top level, previous'}}, {'attribute': {'change': 'p.Asn678Ser', 'value': 'LRG_389p1:p.Asn678Ser', 'type': 'HGVS, protein'}}, {'attribute': {'change': 'p.Asn678Ser', 'value': 'NP_002462.2:p.Asn678Ser', 'type': 'HGVS, protein, RefSeq'}}, {'xref': [{'id': 'SO:0001583', 'db': 'Sequence Ontology', 'status': 'CURRENT'}, {'id': 'NM_002471.3:c.2033A>G', 'db': 'RefSeq', 'status': 'CURRENT'}], 'attribute': {'value': 'missense variant', 'type': 'MolecularConsequence'}}, {'attribute': {'value': 'N678S', 'type': 'ProteinChange1LetterCode'}}], 'xref': [{'id': '515726230', 'db': 'dbSNP', 'status': 'CURRENT', 'type': 'rs'}], 'name': [{'elementValue': {'value': 'NM_002471.3(MYH6):c.2033A>G (p.Asn678Ser)', 'type': 'Preferred'}}], 'type': 'single nucleotide variant'}], 'type': 'Variant'}, 'assertion': {'type': 'VARIATION_TO_DISEASE'}, 'clinicalSignificance': {'reviewStatus': 'CLASSIFIED_BY_SINGLE_SUBMITTER', 'description': 'Likely pathogenic'}, 'observedIn': [{'method': [{'methodType': 'NOT_PROVIDED'}], 'sample': {'affectedStatus': 'not provided', 'species': {'taxonomyId': 9606, 'value': 'human'}, 'origin': 'germline', 'numberTested': 1}, 'observedData': [{'id': 3592924, 'attribute': {'value': 'not provided', 'type': 'Description'}}]}, {'method': [{'methodType': 'NOT_PROVIDED'}], 'sample': {'affectedStatus': 'not provided', 'species': {'taxonomyId': 9606, 'value': 'human'}, 'origin': 'inherited', 'numberTested': 1}, 'observedData': [{'id': 3642470, 'attribute': {'value': 'not provided', 'type': 'Description'}}]}], 'clinVarAccession': {'acc': 'RCV000128628', 'dateUpdated': 1412982000000, 'version': 1, 'type': 'RCV'}, 'recordStatus': 'current'}, 'title': 'NM_002471.3(MYH6):c.2033A>G (p.Asn678Ser) AND Malformation of the heart', 'recordStatus': 'current', 'id': 3829163, 'clinVarAssertion': [{'id': 299693, 'traitSet': {'trait': [{'xref': [{'id': 'C14.240.400', 'db': 'MESH', 'status': 'UNDER_REVIEW'}], 'name': [{'elementValue': {'value': 'not provided', 'type': 'Preferred'}}], 'type': 'Disease'}, {'xref': [{'id': 'C14.280.400', 'db': 'MESH', 'status': 'UNDER_REVIEW'}], 'name': [{'elementValue': {'value': 'not provided', 'type': 'Preferred'}}], 'type': 'Disease'}, {'xref': [{'id': 'C16.131.240.400', 'db': 'MESH', 'status': 'UNDER_REVIEW'}], 'name': [{'elementValue': {'value': 'not provided', 'type': 'Preferred'}}], 'type': 'Disease'}], 'type': 'Disease'}, 'measureSet': {'measure': [{'attributeSet': [{'attribute': {'value': 'NM_002471.3:c.2033A>G', 'type': 'HGVS'}}], 'type': 'Variation'}], 'type': 'Variant'}, 'comment': [{'value': 'this variant was identified in a 4 generation family with multiple members affected with congenital heart defects (multiple types)'}], 'assertion': {'type': 'variation to disease'}, 'submissionName': 'cardiac_targeted_resequencing', 'clinVarSubmissionID': {'submitterDate': 1405465200000, 'submitter': 'Laboratory for Genetics of Human Development Center for Human Genetics, Catholic University of Leuven', 'localKey': 'NM_002471.3:c.2033A>G|LABGENHUDEVKULEUVEN'}, 'observedIn': [{'method': [{'methodType': 'NOT_PROVIDED'}], 'sample': {'affectedStatus': 'not provided', 'species': {'value': 'human'}, 'origin': 'germline', 'numberTested': 1}, 'observedData': [{'attribute': {'value': 'not provided', 'type': 'Description'}}]}, {'method': [{'methodType': 'NOT_PROVIDED'}], 'sample': {'affectedStatus': 'not provided', 'species': {'value': 'human'}, 'origin': 'inherited', 'numberTested': 1}, 'observedData': [{'attribute': {'value': 'not provided', 'type': 'Description'}}]}], 'recordStatus': 'current', 'clinVarAccession': {'orgID': 500109, 'acc': 'SCV000172246', 'dateUpdated': 1405638000000, 'version': 1, 'type': 'SCV'}, 'clinicalSignificance': {'reviewStatus': 'CLASSIFIED_BY_SINGLE_SUBMITTER', 'comment': [{'value': 'Converted during submission to Likely pathogenic.'}], 'description': ['probable-pathogenic']}}]})
    consequenceType = consequence_type.ConsequenceType({'ENSG00000197616'}, {"missense_variant"})
    ensembl_gene_id = "ENSG00000197616"
    ensembl_gene_id_uri = "http://identifiers.org/ensembl/ENSG00000197616"
    ensembl_gene_id_uris = {'http://identifiers.org/ensembl/ENSG00000135486'}
    measure_set_refs_list = []
    n_more_than_one_efo_term = 0
    observed_refs_list = []
    rcv_to_gene_evidence_codes = ['http://identifiers.org/eco/cttv_mapping_pipeline']
    record = {'end': 23866396, 'clinvarSet': {'referenceClinVarAssertion': {'dateLastUpdated': 1412982000000, 'dateCreated': 1405551600000, 'id': 300575, 'traitSet': {'trait': [{'xref': [{'id': 'C0018798', 'db': 'MedGen', 'status': 'CURRENT'}, {'id': '140500', 'db': 'OMIM', 'status': 'CURRENT', 'type': 'MIM'}, {'id': '234750', 'db': 'OMIM', 'status': 'CURRENT', 'type': 'MIM'}], 'name': [{'elementValue': {'value': 'Malformation of the heart', 'type': 'Preferred'}}, {'xref': [{'id': '140500', 'db': 'OMIM', 'status': 'CURRENT', 'type': 'MIM'}, {'id': '234750', 'db': 'OMIM', 'status': 'CURRENT', 'type': 'MIM'}], 'elementValue': {'value': 'HEART, MALFORMATION OF', 'type': 'Alternate'}}, {'elementValue': {'value': 'Congenital heart defect', 'type': 'Alternate'}}], 'id': 15882, 'type': 'Disease'}], 'id': 16428, 'type': 'Disease'}, 'measureSet': {'name': [{'elementValue': {'value': 'NM_002471.3(MYH6):c.2033A>G (p.Asn678Ser)', 'type': 'preferred name'}}], 'id': 139663, 'measure': [{'id': 150127, 'cytogeneticLocation': ['14q11.2'], 'measureRelationship': [{'symbol': [{'elementValue': {'value': 'MYH6', 'type': 'Preferred'}}], 'xref': [{'id': '4624', 'db': 'Gene', 'status': 'CURRENT'}, {'id': '160710', 'db': 'OMIM', 'status': 'CURRENT', 'type': 'MIM'}], 'type': 'variant in gene', 'name': [{'elementValue': {'value': 'myosin, heavy chain 6, cardiac muscle, alpha', 'type': 'Preferred'}}], 'sequenceLocation': [{'stop': 23877485, 'accession': 'NC_000014.8', 'assembly': 'GRCh37', 'chr': '14', 'strand': '-', 'start': 23851198}, {'stop': 23409621, 'accession': 'NC_000014.9', 'assembly': 'GRCh38', 'chr': '14', 'strand': '-', 'start': 23380732}]}], 'sequenceLocation': [{'stop': 23866396, 'accession': 'NC_000014.8', 'variantLength': 1, 'assembly': 'GRCh37', 'referenceAllele': 'T', 'chr': '14', 'start': 23866396, 'alternateAllele': 'C'}, {'stop': 23397187, 'accession': 'NC_000014.9', 'variantLength': 1, 'assembly': 'GRCh38', 'referenceAllele': 'T', 'chr': '14', 'start': 23397187, 'alternateAllele': 'C'}], 'attributeSet': [{'attribute': {'change': 'c.2033A>G', 'value': 'LRG_389t1:c.2033A>G', 'type': 'HGVS, coding, LRG'}}, {'attribute': {'change': 'c.2033A>G', 'value': 'NM_002471.3:c.2033A>G', 'type': 'HGVS, coding, RefSeq'}}, {'attribute': {'change': 'g.16091A>G', 'value': 'LRG_389:g.16091A>G', 'type': 'HGVS, genomic, LRG'}}, {'attribute': {'change': 'g.16091A>G', 'value': 'NG_023444.1:g.16091A>G', 'type': 'HGVS, genomic, RefSeqGene'}}, {'attribute': {'change': 'g.23397187T>C', 'value': 'NC_000014.9:g.23397187T>C', 'integerValue': 38, 'type': 'HGVS, genomic, top level'}}, {'attribute': {'change': 'g.23866396T>C', 'value': 'NC_000014.8:g.23866396T>C', 'integerValue': 37, 'type': 'HGVS, genomic, top level, previous'}}, {'attribute': {'change': 'p.Asn678Ser', 'value': 'LRG_389p1:p.Asn678Ser', 'type': 'HGVS, protein'}}, {'attribute': {'change': 'p.Asn678Ser', 'value': 'NP_002462.2:p.Asn678Ser', 'type': 'HGVS, protein, RefSeq'}}, {'xref': [{'id': 'SO:0001583', 'db': 'Sequence Ontology', 'status': 'CURRENT'}, {'id': 'NM_002471.3:c.2033A>G', 'db': 'RefSeq', 'status': 'CURRENT'}], 'attribute': {'value': 'missense variant', 'type': 'MolecularConsequence'}}, {'attribute': {'value': 'N678S', 'type': 'ProteinChange1LetterCode'}}], 'xref': [{'id': '515726230', 'db': 'dbSNP', 'status': 'CURRENT', 'type': 'rs'}], 'name': [{'elementValue': {'value': 'NM_002471.3(MYH6):c.2033A>G (p.Asn678Ser)', 'type': 'Preferred'}}], 'type': 'single nucleotide variant'}], 'type': 'Variant'}, 'assertion': {'type': 'VARIATION_TO_DISEASE'}, 'clinicalSignificance': {'reviewStatus': 'CLASSIFIED_BY_SINGLE_SUBMITTER', 'description': 'Likely pathogenic'}, 'observedIn': [{'method': [{'methodType': 'NOT_PROVIDED'}], 'sample': {'affectedStatus': 'not provided', 'species': {'taxonomyId': 9606, 'value': 'human'}, 'origin': 'germline', 'numberTested': 1}, 'observedData': [{'id': 3592924, 'attribute': {'value': 'not provided', 'type': 'Description'}}]}, {'method': [{'methodType': 'NOT_PROVIDED'}], 'sample': {'affectedStatus': 'not provided', 'species': {'taxonomyId': 9606, 'value': 'human'}, 'origin': 'inherited', 'numberTested': 1}, 'observedData': [{'id': 3642470, 'attribute': {'value': 'not provided', 'type': 'Description'}}]}], 'clinVarAccession': {'acc': 'RCV000128628', 'dateUpdated': 1412982000000, 'version': 1, 'type': 'RCV'}, 'recordStatus': 'current'}, 'title': 'NM_002471.3(MYH6):c.2033A>G (p.Asn678Ser) AND Malformation of the heart', 'recordStatus': 'current', 'id': 3829163, 'clinVarAssertion': [{'id': 299693, 'traitSet': {'trait': [{'xref': [{'id': 'C14.240.400', 'db': 'MESH', 'status': 'UNDER_REVIEW'}], 'name': [{'elementValue': {'value': 'not provided', 'type': 'Preferred'}}], 'type': 'Disease'}, {'xref': [{'id': 'C14.280.400', 'db': 'MESH', 'status': 'UNDER_REVIEW'}], 'name': [{'elementValue': {'value': 'not provided', 'type': 'Preferred'}}], 'type': 'Disease'}, {'xref': [{'id': 'C16.131.240.400', 'db': 'MESH', 'status': 'UNDER_REVIEW'}], 'name': [{'elementValue': {'value': 'not provided', 'type': 'Preferred'}}], 'type': 'Disease'}], 'type': 'Disease'}, 'measureSet': {'measure': [{'attributeSet': [{'attribute': {'value': 'NM_002471.3:c.2033A>G', 'type': 'HGVS'}}], 'type': 'Variation'}], 'type': 'Variant'}, 'comment': [{'value': 'this variant was identified in a 4 generation family with multiple members affected with congenital heart defects (multiple types)'}], 'assertion': {'type': 'variation to disease'}, 'submissionName': 'cardiac_targeted_resequencing', 'clinVarSubmissionID': {'submitterDate': 1405465200000, 'submitter': 'Laboratory for Genetics of Human Development Center for Human Genetics, Catholic University of Leuven', 'localKey': 'NM_002471.3:c.2033A>G|LABGENHUDEVKULEUVEN'}, 'observedIn': [{'method': [{'methodType': 'NOT_PROVIDED'}], 'sample': {'affectedStatus': 'not provided', 'species': {'value': 'human'}, 'origin': 'germline', 'numberTested': 1}, 'observedData': [{'attribute': {'value': 'not provided', 'type': 'Description'}}]}, {'method': [{'methodType': 'NOT_PROVIDED'}], 'sample': {'affectedStatus': 'not provided', 'species': {'value': 'human'}, 'origin': 'inherited', 'numberTested': 1}, 'observedData': [{'attribute': {'value': 'not provided', 'type': 'Description'}}]}], 'recordStatus': 'current', 'clinVarAccession': {'orgID': 500109, 'acc': 'SCV000172246', 'dateUpdated': 1405638000000, 'version': 1, 'type': 'SCV'}, 'clinicalSignificance': {'reviewStatus': 'CLASSIFIED_BY_SINGLE_SUBMITTER', 'comment': [{'value': 'Converted during submission to Likely pathogenic.'}], 'description': ['probable-pathogenic']}}]}, 'reference': 'T', 'chromosome': '14', 'alternate': 'C', 'start': 23866396, 'annot': {'end': 23866396, 'hgvs': ['ENST00000405093.3:c.2033A>G', 'ENSP00000386041.3:p.Asn678Ser', 'ENST00000356287.3:c.2033A>G', 'ENSP00000348634.3:p.Asn678Ser'], 'alternativeAllele': 'C', 'chromosome': '14', 'referenceAllele': 'T', 'start': 23866396, 'consequenceTypes': [{'ensemblTranscriptId': 'ENST00000405093', 'biotype': 'protein_coding', 'soTerms': [{'soName': 'missense_variant', 'soAccession': 'SO:0001583'}], 'geneName': 'MYH6', 'codon': 'aAt/aGt', 'proteinSubstitutionScores': [{'description': 'deleterious', 'source': 'Sift', 'score': 0.0}, {'description': 'probably_damaging', 'source': 'Polyphen', 'score': 0.999}], 'cDnaPosition': 2104, 'aaPosition': 678, 'aaChange': 'N/S', 'cdsPosition': 2033, 'ensemblGeneId': 'ENSG00000197616', 'strand': '-'}, {'ensemblTranscriptId': 'ENST00000356287', 'biotype': 'protein_coding', 'soTerms': [{'soName': 'missense_variant', 'soAccession': 'SO:0001583'}], 'geneName': 'MYH6', 'codon': 'aAt/aGt', 'proteinSubstitutionScores': [{'description': 'deleterious', 'source': 'Sift', 'score': 0.0}, {'description': 'probably_damaging', 'source': 'Polyphen', 'score': 0.999}], 'cDnaPosition': 2063, 'aaPosition': 678, 'aaChange': 'N/S', 'cdsPosition': 2033, 'ensemblGeneId': 'ENSG00000197616', 'strand': '-'}]}}
    rs = "rs515726230"
    trait_counter = 0
    traits_ref_list = [[]]
    traits = {'http://www.ebi.ac.uk/efo/EFO_0003840'}
    unrecognised_clin_sigs = set()

    test_args_1 = (efo_list, clin_sig, clin_sig_2_activity, clinvarRecord, consequenceType, ensembl_gene_id,
                 ensembl_gene_id_uri, ensembl_gene_id_uris, measure_set_refs_list, n_more_than_one_efo_term,
                 observed_refs_list, rcv_to_gene_evidence_codes, record, rs, trait_counter, traits_ref_list, traits,
                 unrecognised_clin_sigs)

    return test_args_1


class GetCttvGeneticsEvidenceStringTest(unittest.TestCase):
    def setUp(self):
        test_args = get_args_GetCttvGeneticsEvidenceStringTest()
        self.evidence_string, n_more_than_one_efo_term = clinvar_to_evidence_strings.get_cttv_genetics_evidence_string(*test_args)

    def test_evidence_string(self):
        test_dict = {'type': 'genetic_association', 'variant': {'type': 'snp single', 'id': ['http://identifiers.org/dbsnp/rs515726230']}, 'unique_association_fields': {'phenotype': 'http://www.orpha.net/ORDO/Orphanet_88991', 'alleleOrigin': 'germline', 'clinvarAccession': 'RCV000128628', 'gene': 'ENSG00000197616'}, 'access_level': 'public', 'sourceID': 'eva', 'target': {'id': ['http://identifiers.org/ensembl/ENSG00000197616'], 'activity': 'http://identifiers.org/cttv.activity/predicted_damaging', 'target_type': 'http://identifiers.org/cttv.target/gene_variant'}, 'disease': {'id': ['http://www.orpha.net/ORDO/Orphanet_88991']}, 'evidence': {'variant2disease': {'resource_score': {'type': 'pvalue', 'method': {'url': '', 'description': 'Not provided by data supplier'}, 'value': 1e-07}, 'unique_experiment_reference': 'http://europepmc.org/abstract/MED/0', 'is_associated': True, 'evidence_codes': ['http://purl.obolibrary.org/obo/ECO_0000205'], 'date_asserted': '2014-10-11T00:00:00', 'provenance_type': {'database': {'id': 'EVA', 'dbxref': {'id': 'http://identifiers.org/clinvar', 'url': 'http://identifiers.org/clinvar.record/RCV000128628', 'version': '2015-04'}, 'version': '1.0'}, 'expert': {'statement': 'Primary submitter of data', 'status': True}}, 'urls': [{'url': 'http://www.ncbi.nlm.nih.gov/clinvar/RCV000128628', 'nice_name': 'Further details in ClinVar database'}]}, 'gene2variant': {'date_asserted': '2014-10-11T00:00:00', 'is_associated': True, 'evidence_codes': ['http://identifiers.org/eco/cttv_mapping_pipeline'], 'functional_consequence': 'http://purl.obolibrary.org/obo/SO_0001583', 'provenance_type': {'database': {'id': 'EVA', 'dbxref': {'id': 'http://identifiers.org/clinvar', 'url': 'http://identifiers.org/clinvar.record/RCV000128628', 'version': '2015-04'}, 'version': '1.0'}, 'expert': {'statement': 'Primary submitter of data', 'status': True}}, 'urls': [{'url': 'http://www.ncbi.nlm.nih.gov/clinvar/RCV000128628', 'nice_name': 'Further details in ClinVar database'}]}}, 'validated_against_schema_version': '1.2.2'}
        self.assertEqual(self.evidence_string, test_dict)


############
def get_args_GetCttvSomaticEvidenceStringTest():
    efo_list = ['http://www.ebi.ac.uk/efo/EFO_0003840']
    clin_sig = "likely pathogenic"
    clin_sig_2_activity = {'other': 'http://identifiers.org/cttv.activity/unknown', 'unknown': 'http://identifiers.org/cttv.activity/unknown', 'protective': 'http://identifiers.org/cttv.activity/tolerated_by_target', 'probable-pathogenic': 'http://identifiers.org/cttv.activity/predicted_damaging', 'non-pathogenic': 'http://identifiers.org/cttv.activity/tolerated_by_target', 'benign': 'http://identifiers.org/cttv.activity/tolerated_by_target', 'likely pathogenic': 'http://identifiers.org/cttv.activity/predicted_damaging', 'probable-non-pathogenic': 'http://identifiers.org/cttv.activity/predicted_tolerated', 'pathogenic': 'http://identifiers.org/cttv.activity/damaging_to_target', 'association': 'http://identifiers.org/cttv.activity/damaging_to_target', 'conflicting data from submitters': 'http://identifiers.org/cttv.activity/unknown', 'uncertain significance': 'http://identifiers.org/cttv.activity/unknown', 'likely benign': 'http://identifiers.org/cttv.activity/predicted_tolerated', 'histocompatibility': 'http://identifiers.org/cttv.activity/unknown', 'not provided': 'http://identifiers.org/cttv.activity/unknown', 'untested': 'http://identifiers.org/cttv.activity/unknown', 'confers sensitivity': 'http://identifiers.org/cttv.activity/predicted_damaging', 'drug-response': 'http://identifiers.org/cttv.activity/unknown', 'risk factor': 'http://identifiers.org/cttv.activity/predicted_damaging'}
    clinvarRecord = clinvar_record.ClinvarRecord({'title': 'NM_031157.2(HNRNPA1):c.1054C>T (p.Arg352Ter) AND Chronic progressive multiple sclerosis', 'clinVarAssertion': [{'measureSet': {'type': 'Variant', 'measure': [{'attributeSet': [{'attribute': {'type': 'HGVS', 'value': 'NM_002136.2:c.898C>T'}}], 'type': 'Variation'}]}, 'id': 286280, 'traitSet': {'type': 'Disease', 'trait': [{'name': [{'elementValue': {'type': 'Preferred', 'value': 'not provided'}}], 'xref': [{'id': 'C10.114.375.500.200', 'db': 'MESH', 'status': 'UNDER_REVIEW'}], 'type': 'Disease'}]}, 'submissionName': 'HNRNPA1_VARIATION', 'clinVarAccession': {'orgID': 504916, 'type': 'SCV', 'dateUpdated': 1409698800000, 'acc': 'SCV000154960', 'version': 1}, 'assertion': {'type': 'variation to disease'}, 'clinVarSubmissionID': {'submitter': 'Demyelinating Disease Laboratories, VA Medical Center and University of Tennessee', 'submitterDate': 1402268400000, 'localKey': 'NP_002127.1:c.898C>T'}, 'clinicalSignificance': {'description': ['probable-pathogenic'], 'reviewStatus': 'CLASSIFIED_BY_SINGLE_SUBMITTER', 'comment': [{'value': 'Converted during submission to Likely pathogenic.'}]}, 'recordStatus': 'current', 'observedIn': [{'observedData': [{'attribute': {'type': 'Description', 'value': 'not provided'}}], 'sample': {'affectedStatus': 'not provided', 'species': {'value': 'human'}, 'origin': 'somatic', 'numberTested': 1}, 'method': [{'methodType': 'NOT_PROVIDED'}]}]}], 'referenceClinVarAssertion': {'measureSet': {'name': [{'elementValue': {'type': 'preferred name', 'value': 'NM_031157.2(HNRNPA1):c.1054C>T (p.Arg352Ter)'}}], 'measure': [{'attributeSet': [{'attribute': {'type': 'HGVS, coding, RefSeq', 'change': 'c.1054C>T', 'value': 'NM_031157.2:c.1054C>T'}}, {'attribute': {'type': 'HGVS, coding, RefSeq', 'change': 'c.898C>T', 'value': 'NM_002136.2:c.898C>T'}}, {'attribute': {'type': 'HGVS, genomic, RefSeqGene', 'change': 'g.8255C>T', 'value': 'NG_033830.1:g.8255C>T'}}, {'attribute': {'type': 'HGVS, genomic, top level', 'integerValue': 38, 'change': 'g.54283958C>T', 'value': 'NC_000012.12:g.54283958C>T'}}, {'attribute': {'type': 'HGVS, genomic, top level, previous', 'integerValue': 37, 'change': 'g.54677742C>T', 'value': 'NC_000012.11:g.54677742C>T'}}, {'attribute': {'type': 'HGVS, protein, RefSeq', 'change': 'p.Arg300Ter', 'value': 'NP_002127.1:p.Arg300Ter'}}, {'attribute': {'type': 'HGVS, protein, RefSeq', 'change': 'p.Arg352Ter', 'value': 'NP_112420.1:p.Arg352Ter'}}, {'xref': [{'id': 'SO:0001587', 'db': 'Sequence Ontology', 'status': 'CURRENT'}, {'id': 'NM_031157.2:c.1054C>T', 'db': 'RefSeq', 'status': 'CURRENT'}], 'attribute': {'type': 'MolecularConsequence', 'value': 'nonsense'}}, {'attribute': {'type': 'ProteinChange1LetterCode', 'value': 'R352*'}}, {'attribute': {'type': 'ProteinChange1LetterCode', 'value': 'R300*'}}], 'cytogeneticLocation': ['12q13.13'], 'sequenceLocation': [{'variantLength': 1, 'chr': '12', 'start': 54677742, 'assembly': 'GRCh37', 'stop': 54677742, 'alternateAllele': 'T', 'referenceAllele': 'C', 'accession': 'NC_000012.11'}, {'variantLength': 1, 'chr': '12', 'start': 54283958, 'assembly': 'GRCh38', 'stop': 54283958, 'alternateAllele': 'T', 'referenceAllele': 'C', 'accession': 'NC_000012.12'}], 'name': [{'elementValue': {'type': 'Preferred', 'value': 'NM_031157.2(HNRNPA1):c.1054C>T (p.Arg352Ter)'}}], 'xref': [{'type': 'rs', 'id': '483353037', 'db': 'dbSNP', 'status': 'CURRENT'}], 'type': 'single nucleotide variant', 'measureRelationship': [{'symbol': [{'elementValue': {'type': 'Preferred', 'value': 'HNRNPA1'}}], 'name': [{'elementValue': {'type': 'Preferred', 'value': 'heterogeneous nuclear ribonucleoprotein A1'}}], 'xref': [{'id': '3178', 'db': 'Gene', 'status': 'CURRENT'}, {'type': 'MIM', 'id': '164017', 'db': 'OMIM', 'status': 'CURRENT'}], 'type': 'variant in gene', 'sequenceLocation': [{'chr': '12', 'start': 54674487, 'assembly': 'GRCh37', 'stop': 54679029, 'strand': '+', 'accession': 'NC_000012.11'}, {'chr': '12', 'start': 54280695, 'assembly': 'GRCh38', 'stop': 54287086, 'strand': '+', 'accession': 'NC_000012.12'}]}], 'id': 139332}], 'type': 'Variant', 'id': 135606}, 'dateCreated': 1402354800000, 'traitSet': {'type': 'Disease', 'trait': [{'name': [{'elementValue': {'type': 'Preferred', 'value': 'Chronic progressive multiple sclerosis'}, 'xref': [{'id': '230373008', 'db': 'SNOMED CT', 'status': 'CURRENT'}]}], 'xref': [{'id': 'C0393665', 'db': 'MedGen', 'status': 'CURRENT'}], 'type': 'Disease', 'id': 18795}], 'id': 13814}, 'dateLastUpdated': 1412982000000, 'clinVarAccession': {'type': 'RCV', 'dateUpdated': 1412982000000, 'acc': 'RCV000122455', 'version': 1}, 'assertion': {'type': 'VARIATION_TO_DISEASE'}, 'clinicalSignificance': {'description': 'Likely pathogenic', 'reviewStatus': 'CLASSIFIED_BY_SINGLE_SUBMITTER'}, 'recordStatus': 'current', 'observedIn': [{'observedData': [{'id': 3574852, 'attribute': {'type': 'Description', 'value': 'not provided'}}], 'sample': {'affectedStatus': 'not provided', 'species': {'taxonomyId': 9606, 'value': 'human'}, 'origin': 'somatic', 'numberTested': 1}, 'method': [{'methodType': 'NOT_PROVIDED'}]}], 'id': 286319}, 'id': 3823230, 'recordStatus': 'current'})
    ensembl_gene_id = "ENSG00000135486"
    ensembl_gene_id_uri = "http://identifiers.org/ensembl/ENSG00000135486"
    ensembl_gene_id_uris = set()
    measure_set_refs_list = []
    n_more_than_one_efo_term = 0
    observed_regs_list = []
    trait_counter = 0
    trait_refs_list = [[]]
    traits = set()
    unrecognised_clin_sigs = set()
    consequenceType = consequence_type.ConsequenceType(ensembl_gene_ids={'ENSG00000135486'}, so_names={"stop_gained"})

    test_args_1 = (efo_list, clin_sig, clin_sig_2_activity, clinvarRecord, ensembl_gene_id, ensembl_gene_id_uri,
                   ensembl_gene_id_uris, measure_set_refs_list, n_more_than_one_efo_term, observed_regs_list,
                   trait_counter, trait_refs_list, traits, unrecognised_clin_sigs, consequenceType)

    return test_args_1


class GetCttvSomaticEvidenceStringTest(unittest.TestCase):
    def setUp(self):
        test_args = get_args_GetCttvSomaticEvidenceStringTest()
        self.evidence_string, n_more_than_one_efo_term = clinvar_to_evidence_strings.get_cttv_somatic_evidence_string(*test_args)

    def test_evidence_string(self):
        test_dict = {'evidence': {'provenance_type': {'expert': {'status': True, 'statement': 'Primary submitter of data'}, 'database': {'dbxref': {'id': 'http://identifiers.org/clinvar', 'url': 'http://identifiers.org/clinvar.record/RCV000122455', 'version': '2015-04'}, 'id': 'EVA', 'version': '1.0'}}, 'resource_score': {'value': 1, 'type': 'probability'}, 'date_asserted': '2014-10-11T00:00:00', 'is_associated': True, 'known_mutations': [{'functional_consequence': 'http://purl.obolibrary.org/obo/SO_0001587', 'preferred_name': 'stop_gained'}], 'urls': [{'url': 'http://www.ncbi.nlm.nih.gov/clinvar/RCV000122455', 'nice_name': 'Further details in ClinVar database'}], 'evidence_codes': ['http://purl.obolibrary.org/obo/ECO_0000205']}, 'type': 'somatic_mutation', 'sourceID': 'eva_somatic', 'access_level': 'public', 'disease': {'id': ['http://www.ebi.ac.uk/efo/EFO_0003840']}, 'validated_against_schema_version': '1.2.2', 'target': {'activity': 'http://identifiers.org/cttv.activity/predicted_damaging', 'id': ['http://identifiers.org/ensembl/ENSG00000135486'], 'target_type': 'http://identifiers.org/cttv.target/gene_variant'}, 'unique_association_fields': {'phenotype': 'http://www.ebi.ac.uk/efo/EFO_0003840', 'alleleOrigin': 'somatic', 'gene': 'ENSG00000135486', 'clinvarAccession': 'RCV000122455'}}
        self.assertEqual(self.evidence_string, test_dict)

############


class GetCttvVariantTypeTest(unittest.TestCase):
    def setUp(self):
        self.record_single_a = ({"reference": "AGAGACGTACGTACGTACGTACGTACGTACGTACGTACG", "alternate": "C"}, "snp single")
        self.record_single_b = ({"reference": "A", "alternate": "C"}, "snp single")
        self.record_single_c = ({"reference": "AGAGACGTACGTACGTACGTACGTACGTACGTACGTACG", "alternate": "AGAGACGTACGTACGTACGTACGTACGTACGTACGTACG"}, "snp single")

        self.test_records_singles = [self.record_single_a, self.record_single_b, self.record_single_c]

        self.record_structurals_a = ({"reference": "AGAGACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT", "alternate": "C"},
                                "structural variant")
        self.record_structurals_b = ({"reference": "A", "alternate": "AGAGACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"},
                                "structural variant")
        self.record_structurals_c = ({"reference": "AGAGACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT",
                                 "alternate": "AGAGACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"},
                                "structural variant")

        self.test_records_structurals = [self.record_structurals_a, self.record_structurals_b, self.record_structurals_c]

    def test_get_cttv_variant_type_singles(self):
        for record in self.test_records_singles:
            self.assertEqual(clinvar_to_evidence_strings.get_cttv_variant_type(record[0]), record[1])

    def test_get_cttv_variant_type_structurals(self):
        for record in self.test_records_structurals:
            self.assertEqual(clinvar_to_evidence_strings.get_cttv_variant_type(record[0]), record[1])


class LoadEfoMappingTest(unittest.TestCase):
    def setUp(self):
        ignore_file = utilities.get_resource_file("eva_cttv_pipeline", "resources/testing/ignore_file.txt")
        efo_file = utilities.get_resource_file("eva_cttv_pipeline", "resources/testing/ClinVar_Traits_EFO_090915.xls")

        self.trait_2_efo, self.unavailable_efo = clinvar_to_evidence_strings.load_efo_mapping(efo_file)
        self.trait_2_efo_w_ignore, self.unavailable_efo_w_ignore = clinvar_to_evidence_strings.load_efo_mapping(efo_file, ignore_terms_file=ignore_file)

    def test_just_mapping_trait_2_efo(self):
        self.assertEqual(len(self.trait_2_efo), 3819)

    def test_w_ignore_trait_2_efo(self):
        self.assertEqual(len(self.trait_2_efo_w_ignore), 3528)


class GetUnmappedUrlTest(unittest.TestCase):
    def test_orphanet(self):
        url = "http://www.orpha.net/ORDO/Orphanet_2670"
        self.assertEqual(clinvar_to_evidence_strings.get_unmapped_url(url), "http://purl.bioontology.org/ORDO/Orphanet_2670")

    def test_hp(self):
        url = "http://purl.obolibrary.org/obo/HP_0000545"
        self.assertEqual(clinvar_to_evidence_strings.get_unmapped_url(url), "http://purl.bioontology.org/obo/HP_0000545")

    def test_bad_url(self):
        #TODO currently the function exits execution with bad urls
        pass


class GetTermsFromFileTest(unittest.TestCase):
    #TODO do the same for adapt terms file?
    def setUp(self):
        ignore_file = utilities.get_resource_file("eva_cttv_pipeline", "resources/testing/ignore_file.txt")
        self.ignore_terms = clinvar_to_evidence_strings.get_terms_from_file(ignore_file)

    def test_length(self):
        self.assertEqual(len(self.ignore_terms), 218)

    def test_head(self):
        self.assertEqual(self.ignore_terms[0], "http://purl.obolibrary.org/obo/HP_0011677")

    def test_tail(self):
        self.assertEqual(self.ignore_terms[-1], "http://www.orpha.net/ORDO/Orphanet_120795")

    def test_no_file(self):
        self.assertEqual(clinvar_to_evidence_strings.get_terms_from_file(None), [])


# def temp():
#     ignore_file = utilities.get_resource_file("eva_cttv_pipeline", "resources/testing/ignore_file.txt")
#     efo_file = utilities.get_resource_file("eva_cttv_pipeline", "resources/testing/ClinVar_Traits_EFO_090915.xls")
#
#     trait_2_efo, unavailable_efo = clinvar_to_evidence_strings.load_efo_mapping(efo_file)
#
#     print([])
#     print(len(unavailable_efo))
#
#
# temp()