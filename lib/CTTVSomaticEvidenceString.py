__author__ = 'Javier Lopez: javild@gmail.com'

from CTTVEvidenceString import CTTVEvidenceString
import json
import jsonschema
import os

SCHEMA_FILE = os.path.dirname(__file__) + "/resources/schema_local/literature_curated.local.json"
# SCHEMA_FILE = os.path.dirname(__file__) + "/resources/schema/literature_curated.json"


class CTTVSomaticEvidenceString(CTTVEvidenceString):
    schema = json.loads(open(SCHEMA_FILE, 'r').read())

    def __init__(self):

        CTTVEvidenceString.__init__(self,
                                    {'sourceID': 'eva_somatic',
                                     'validated_against_schema_version': '1.2.1',
                                     'type': 'somatic_mutation',
                                     'access_level': 'public',
                                     'unique_association_fields': {},
                                     'target': {
                                         'id': [],  # ENSEMBL ID
                                         'target_type': 'http://identifiers.org/cttv.target/gene_variant',
                                         'activity': None  # http://identifiers.org/cttv.activity/predicted_damaging
                                     },
                                     'evidence': {
                                         'resource_score': {
                                             'type': 'probability',
                                             'value': 1
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

                                     'disease': {
                                         'id': [],  # EFO terms
                                     }
                                     }
                                    )

    def setDbxrefUrl(self, url):
        self['evidence']['provenance_type']['database']['dbxref']['url'] = url
        self['evidence']['provenance_type']['database']['dbxref']['url'] = url

    def setUrl(self, url):
        self['evidence']['urls'][0]['url'] = url
        self['evidence']['urls'][0]['url'] = url

    def setLiterature(self, referenceList):
        if 'literature' not in self['evidence']['provenance_type']:
            self['evidence']['provenance_type']['literature'] = {
                'references': [{'lit_id': reference} for reference in referenceList]}
        else:
            self['evidence']['provenance_type']['literature']['references'] = [{'lit_id': reference} for reference in
                                                                               referenceList]

    def setAssociation(self, isAssociated):
        self['evidence']['is_associated'] = isAssociated

    def validate(self):
        jsonschema.validate(self, CTTVSomaticEvidenceString.schema, format_checker=jsonschema.FormatChecker())
        self.getDisease().isObsolete()

    def setDate(self, dateString):
        self['evidence']['date_asserted'] = dateString
