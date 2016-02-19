import sys

__author__ = 'Javier Lopez: javild@gmail.com'

import json

import eva_cttv_pipeline.utilities as utilities
import jsonschema

from eva_cttv_pipeline.CTTVEvidenceString import CTTVEvidenceString

SCHEMA_FILE = utilities.get_resource_file("eva_cttv_pipeline", "resources/schema_local/src/literature_curated.json")
# SCHEMA_FILE = utilities.get_resource_file("eva_cttv_pipeline", "resources/json_schema/src/literature_curated.json")


class CTTVSomaticEvidenceString(CTTVEvidenceString):
    schema = json.loads(open(SCHEMA_FILE, 'r').read())

    def __init__(self):

         CTTVEvidenceString.__init__(self,
                                    {'sourceID': 'eva_somatic',
                                     'validated_against_schema_version': '1.2.2',
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
                                         'known_mutations': [
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

    def setUrl(self, url):
        self['evidence']['urls'][0]['url'] = url

    def set_evidence_literature(self, referenceList):
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

    def add_known_mutation(self, new_functional_consequence):
        new_known_mutation = {'functional_consequence': new_functional_consequence}
        self['evidence']['known_mutations'].append(new_known_mutation)

    def set_known_mutations(self, consequence_type):
        for so_term in consequence_type.get_so_terms():
            if so_term.getAccession():
                new_functional_consequence = "http://purl.obolibrary.org/obo/" + so_term.getAccession().replace(':', '_')
            else:
                new_functional_consequence = 'http://targetvalidation.org/sequence/' + so_term.getName()
            self.add_known_mutation(new_functional_consequence)
