__author__ = 'Javier Lopez: javild@gmail.com'

import json
import os

import jsonschema

from eva_cttv_pipeline.CTTVEvidenceString import CTTVEvidenceString

SCHEMA_FILE = os.path.dirname(__file__) + "/resources/schema_local/genetics.local.json"
# SCHEMA_FILE = os.path.dirname(__file__) + "/resources/schema/genetics.json"


class CTTVGeneticsEvidenceString(CTTVEvidenceString):
    schema = json.loads(open(SCHEMA_FILE, 'r').read())

    def __init__(self):
        CTTVEvidenceString.__init__(self,
                                    {'sourceID': 'eva',
                                     'validated_against_schema_version': '1.2.1',
                                     'type': 'genetic_association',
                                     'access_level': 'public',
                                     'unique_association_fields': {},
                                     'target': {
                                         'id': [],  # ENSEMBL ID
                                         'target_type': 'http://identifiers.org/cttv.target/gene_variant',
                                         'activity': None  # http://identifiers.org/cttv.activity/predicted_damaging
                                     },
                                     'variant': {
                                         'id': [],  # rs
                                         'type': None  # snp single
                                     },
                                     'evidence': {
                                         'gene2variant': {
                                             'provenance_type': {
                                                 'expert': {
                                                     'status': True,
                                                     'statement': 'Primary submitter of data'
                                                 },
                                                 'database': {
                                                     'dbxref': {
                                                         'id': 'http://identifiers.org/clinvar',
                                                         'url': None,
                                                         'version': '2015-04',
                                                     },
                                                     'id': 'EVA',
                                                     'version': '1.0'

                                                 }
                                             },
                                             'is_associated': True,
                                             'date_asserted': None,
                                             'evidence_codes': [],
                                             'functional_consequence': None,
                                             'urls': [
                                                 {
                                                     'nice_name': 'Further details in ClinVar database',
                                                     'url': None  # http://www.ncbi.nlm.nih.gov/clinvar/RCV000038347
                                                 }
                                             ]
                                         },
                                         'variant2disease': {
                                             'resource_score': {
                                                 'type': 'pvalue',
                                                 'value': 0.0000001,
                                                 'method': {
                                                     'url': '',
                                                     'description': 'Not provided by data supplier'
                                                 }
                                             },
                                             'provenance_type': {
                                                 'expert': {
                                                     'status': True,
                                                     'statement': 'Primary submitter of data'
                                                 },
                                                 'database': {
                                                     'dbxref': {
                                                         'id': 'http://identifiers.org/clinvar',
                                                         'url': None,
                                                         'version': '2015-04',
                                                     },
                                                     'id': 'EVA',
                                                     'version': '1.0'

                                                 }
                                             },
                                             'unique_experiment_reference': 'http://europepmc.org/abstract/MED/0',
                                             'is_associated': True,
                                             'date_asserted': None,
                                             'evidence_codes': [
                                                 'http://purl.obolibrary.org/obo/ECO_0000205'
                                             ],
                                             'urls': [
                                                 {
                                                     'nice_name': 'Further details in ClinVar database',
                                                     'url': None  # http://www.ncbi.nlm.nih.gov/clinvar/RCV000038347
                                                 }
                                             ]
                                         },
                                     },
                                     'disease': {
                                         'id': [],  # EFO terms
                                     }
                                     }
                                    )

    def setDbxrefUrl(self, url):
        self['evidence']['gene2variant']['provenance_type']['database']['dbxref']['url'] = url
        self['evidence']['variant2disease']['provenance_type']['database']['dbxref']['url'] = url

    def setUrl(self, url):
        self['evidence']['gene2variant']['urls'][0]['url'] = url
        self['evidence']['variant2disease']['urls'][0]['url'] = url

    def setGene2VariantEvidenceCodes(self, gene2VariantEvidenceCodes):
        self['evidence']['gene2variant']['evidence_codes'] = gene2VariantEvidenceCodes

    def setGene2VariantFunctionalConsequence(self, SOTerm):
        self['evidence']['gene2variant']['functional_consequence'] = SOTerm

    def setVariant2DiseaseLiterature(self, referenceList):
        if 'literature' not in self['evidence']['variant2disease']['provenance_type']:
            self['evidence']['variant2disease']['provenance_type']['literature'] = {
                'references': [{'lit_id': reference} for reference in referenceList]}
        else:
            self['evidence']['variant2disease']['provenance_type']['literature']['references'] = [{'lit_id': reference}
                                                                                                  for reference in
                                                                                                  referenceList]

    def setAssociation(self, isAssociated):
        self['evidence']['gene2variant']['is_associated'] = isAssociated
        self['evidence']['variant2disease']['is_associated'] = isAssociated

    def validate(self):
        jsonschema.validate(self, CTTVGeneticsEvidenceString.schema, format_checker=jsonschema.FormatChecker())
        self.getDisease().isObsolete()

    def setVariant(self, id, type):
        self['variant']['id'].append(id)
        self['variant']['type'] = type

    def setUniqueReference(self, reference):
        self['evidence']['variant2disease']['unique_experiment_reference'] = reference

    def setDate(self, dateString):
        self['evidence']['gene2variant']['date_asserted'] = dateString
        self['evidence']['variant2disease']['date_asserted'] = dateString
