import json

import jsonschema

import eva_cttv_pipeline.config as config
import eva_cttv_pipeline.utilities as utilities
from eva_cttv_pipeline import efo_term

__author__ = 'Javier Lopez: javild@gmail.com'


utilities.check_for_local_schema()
GEN_SCHEMA_FILE = utilities.get_resource_file(__package__, config.local_schema + "/src/genetics.json")
SOM_SCHEMA_FILE = utilities.get_resource_file(__package__, config.local_schema + "/src/literature_curated.json")


class CTTVEvidenceString(dict):
    def __init__(self, a_dictionary):
        dict.__init__(self, a_dictionary)

    def add_unique_association_field(self, key, value):
        self['unique_association_fields'][key] = value

    def set_unique_association_field(self, key, value):
        self['unique_association_fields'][key] = value

    def set_target(self, id, activity):
        self['target']['id'].append(id)
        self['target']['activity'] = activity

    def set_pubmed_refs(self, pubmed_list):
        if pubmed_list is not None and len(pubmed_list) > 0:
            if 'literature' in self['evidence']['provenance_type']:
                self['evidence']['provenance_type']['literature']['pubmed_refs'] = pubmed_list
            else:
                self['evidence']['provenance_type']['literature'] = {'pubmed_refs': pubmed_list}

    @property
    def disease(self):
        return efo_term.EFOTerm(self['disease']['id'][0])

    @disease.setter
    def disease(self, id):
        self['disease']['id'].append(id)

    @property
    def evidence_codes(self):
        return self['evidence']['evidence_codes']

    @evidence_codes.setter
    def evidence_codes(self, ev_code_list):
        self.evidence_codes = ev_code_list

    def set_top_level_literature(self, reference_list):
        if 'literature' not in self:
            self['literature'] = {'references': [{'lit_id': reference} for reference in reference_list]}
        else:
            self['literature']['references'] = [{'lit_id': reference} for reference in reference_list]


class CTTVGeneticsEvidenceString(CTTVEvidenceString):
    schema = json.loads(open(GEN_SCHEMA_FILE, 'r').read())

    def __init__(self):

        CTTVEvidenceString.__init__(self,
                                    {'sourceID': 'eva',
                                     'validated_against_schema_version': '1.2.2',
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

    def set_db_xref_url(self, url):
        self['evidence']['gene2variant']['provenance_type']['database']['dbxref']['url'] = url
        self['evidence']['variant2disease']['provenance_type']['database']['dbxref']['url'] = url

    def set_url(self, url):
        self['evidence']['gene2variant']['urls'][0]['url'] = url
        self['evidence']['variant2disease']['urls'][0]['url'] = url

    @property
    def gene_2_var_ev_codes(self):
        return self['evidence']['gene2variant']['evidence_codes']

    @gene_2_var_ev_codes.setter
    def gene_2_var_ev_codes(self, gene_2_var_ev_codes):
        self.gene_2_var_ev_codes = gene_2_var_ev_codes

    @property
    def gene_2_var_func_consequence(self):
        return self['evidence']['gene2variant']['functional_consequence']

    @gene_2_var_func_consequence.setter
    def gene_2_var_func_consequence(self, so_term):
        self.gene_2_var_func_consequence = so_term

    def set_var_2_disease_literature(self, ref_list):
        if 'literature' not in self['evidence']['variant2disease']['provenance_type']:
            self['evidence']['variant2disease']['provenance_type']['literature'] = {
                'references': [{'lit_id': reference} for reference in ref_list]}
        else:
            self['evidence']['variant2disease']['provenance_type']['literature']['references'] = [{'lit_id': reference}
                                                                                                  for reference in
                                                                                                  ref_list]

    def set_association(self, is_associated):
        self['evidence']['gene2variant']['is_associated'] = is_associated
        self['evidence']['variant2disease']['is_associated'] = is_associated

    def validate(self):
        # try:
        jsonschema.validate(self, CTTVGeneticsEvidenceString.schema, format_checker=jsonschema.FormatChecker())
        # except Exception as e:
        #     print(str(self))
        #     print(e)
        #     traceback.print_stack()
        #     sys.exit()
        self.disease.is_obsolete()

    def set_variant(self, id, type):
        self['variant']['id'].append(id)
        self['variant']['type'] = type

    @property
    def unique_reference(self):
        return self['evidence']['variant2disease']['unique_experiment_reference']

    @unique_reference.setter
    def unique_reference(self, reference):
        self.unique_reference = reference

    def set_date(self, date_string):
        self['evidence']['gene2variant']['date_asserted'] = date_string
        self['evidence']['variant2disease']['date_asserted'] = date_string


class CTTVSomaticEvidenceString(CTTVEvidenceString):
    schema = json.loads(open(SOM_SCHEMA_FILE, 'r').read())

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

    @property
    def db_xref_url(self):
        return self['evidence']['provenance_type']['database']['dbxref']['url']

    @db_xref_url.setter
    def db_xref_url(self, url):
        self.db_xref_url = url

    @property
    def url(self):
        return self['evidence']['urls'][0]['url']

    @url.setter
    def url(self, url):
        self.url = url

    def set_evidence_literature(self, ref_list):
        if 'literature' not in self['evidence']['provenance_type']:
            self['evidence']['provenance_type']['literature'] = {
                'references': [{'lit_id': reference} for reference in ref_list]}
        else:
            self['evidence']['provenance_type']['literature']['references'] = [{'lit_id': reference} for reference in
                                                                               ref_list]

    @property
    def association(self):
        return self['evidence']['is_associated']

    @association.setter
    def association(self, is_associated):
        self.association = is_associated

    def validate(self):
        jsonschema.validate(self, CTTVSomaticEvidenceString.schema, format_checker=jsonschema.FormatChecker())
        self.disease.is_obsolete()

    @property
    def date(self):
        return self['evidence']['date_asserted']

    @date.setter
    def date(self, date_string):
        self.date = date_string

    def add_known_mutation(self, new_functional_consequence, so_name):
        new_known_mutation = {'functional_consequence': new_functional_consequence, 'preferred_name': so_name}
        self['evidence']['known_mutations'].append(new_known_mutation)

    def set_known_mutations(self, consequenceType):
        for so_term in consequenceType.get_so_terms():
            so_name = so_term.get_name()
            if so_term.get_accession():
                new_functional_consequence = "http://purl.obolibrary.org/obo/" + so_term.get_accession().replace(':', '_')
            else:
                new_functional_consequence = 'http://targetvalidation.org/sequence/' + so_term.get_name()
            self.add_known_mutation(new_functional_consequence, so_name)
