from collections import UserDict
import json

import jsonschema

import eva_cttv_pipeline.config as config
import eva_cttv_pipeline.utilities as utilities
from eva_cttv_pipeline import efo_term

__author__ = 'Javier Lopez: javild@gmail.com'


utilities.check_for_local_schema()
GEN_SCHEMA_FILE = utilities.get_resource_file(__package__, config.local_schema + "/src/genetics.json")
SOM_SCHEMA_FILE = utilities.get_resource_file(__package__, config.local_schema + "/src/literature_curated.json")


clin_sig_2_activity = {'other': 'http://identifiers.org/cttv.activity/unknown', 'unknown': 'http://identifiers.org/cttv.activity/unknown', 'protective': 'http://identifiers.org/cttv.activity/tolerated_by_target', 'probable-pathogenic': 'http://identifiers.org/cttv.activity/predicted_damaging', 'non-pathogenic': 'http://identifiers.org/cttv.activity/tolerated_by_target', 'benign': 'http://identifiers.org/cttv.activity/tolerated_by_target', 'likely pathogenic': 'http://identifiers.org/cttv.activity/predicted_damaging', 'probable-non-pathogenic': 'http://identifiers.org/cttv.activity/predicted_tolerated', 'pathogenic': 'http://identifiers.org/cttv.activity/damaging_to_target', 'association': 'http://identifiers.org/cttv.activity/damaging_to_target', 'conflicting data from submitters': 'http://identifiers.org/cttv.activity/unknown', 'uncertain significance': 'http://identifiers.org/cttv.activity/unknown', 'likely benign': 'http://identifiers.org/cttv.activity/predicted_tolerated', 'histocompatibility': 'http://identifiers.org/cttv.activity/unknown', 'not provided': 'http://identifiers.org/cttv.activity/unknown', 'untested': 'http://identifiers.org/cttv.activity/unknown', 'confers sensitivity': 'http://identifiers.org/cttv.activity/predicted_damaging', 'drug-response': 'http://identifiers.org/cttv.activity/unknown', 'risk factor': 'http://identifiers.org/cttv.activity/predicted_damaging'}


def get_cttv_variant_type(ref, alt):
    if len(ref) < 2 and len(alt) < 2:
        cttv_variant_type = 'snp single'
    elif len(ref) > 50 or len(alt) > 50:
        cttv_variant_type = 'structural variant'
    else:
        cttv_variant_type = 'snp single'  # Sam asked for this in his email 21/05/2015
        # cttv_variant_type = 'snp multiple'

    return cttv_variant_type


class CTTVEvidenceString(UserDict):
    def __init__(self, a_dictionary):
        super().__init__(a_dictionary)
        # dict.__init__(a_dictionary)

    def add_unique_association_field(self, key, value):
        self['unique_association_fields'][key] = value

    def _clear_target(self):
        self['target']['id'] = []

    def set_target(self, id, activity):
        self['target']['id'].append(id)
        self['target']['activity'] = activity

    @property
    def disease(self):
        return efo_term.EFOTerm(self['disease']['id'][0])

    @disease.setter
    def disease(self, value):
        self['disease']['id'] = [value]

    @property
    def evidence_codes(self):
        return self['evidence']['evidence_codes']

    @evidence_codes.setter
    def evidence_codes(self, ev_code_list):
        self['evidence']['evidence_codes'] = ev_code_list

    @property
    def top_level_literature(self):
        return self['literature']['references']

    @top_level_literature.setter
    def top_level_literature(self, reference_list):
        self['literature'] = {'references': [{'lit_id': reference} for reference in reference_list]}


class CTTVGeneticsEvidenceString(CTTVEvidenceString):
    schema = json.loads(open(GEN_SCHEMA_FILE, 'r').read())

    def __init__(self, efo_list, clin_sig, clinvarRecord, consequenceType, ensembl_gene_id,
                 ensembl_gene_id_uri, measure_set_refs_list, observed_refs_list, rcv_to_gene_evidence_codes, record, rs,
                 trait_counter, traits_ref_list, unrecognised_clin_sigs):

        a_dictionary = {'sourceID': 'eva',
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

        # CTTVEvidenceString.__init__(self, a_dictionary)
        super().__init__(a_dictionary)

        self.add_unique_association_field('gene', ensembl_gene_id)
        self.add_unique_association_field('clinvarAccession', clinvarRecord.acc)
        self.add_unique_association_field('alleleOrigin', 'germline')
        try:
            self.set_target(ensembl_gene_id_uri, clin_sig_2_activity[clin_sig])
        except KeyError:
            unrecognised_clin_sigs.add(clin_sig)
            self.set_target(ensembl_gene_id_uri, 'http://identifiers.org/cttv.activity/unknown')
        self.set_variant('http://identifiers.org/dbsnp/' + rs, get_cttv_variant_type(record['reference'], record['alternate']))
        self.date = clinvarRecord.date
        self.db_xref_url = 'http://identifiers.org/clinvar.record/' + clinvarRecord.acc
        self.url = 'http://www.ncbi.nlm.nih.gov/clinvar/' + clinvarRecord.acc
        self.association = clin_sig != 'non-pathogenic' and clin_sig != 'probable-non-pathogenic' and clin_sig != 'likely benign' and clin_sig != 'benign'
        self.gene_2_var_ev_codes = rcv_to_gene_evidence_codes
        most_severe_so_term = consequenceType.most_severe_so
        if most_severe_so_term.accession is None:
            self.gene_2_var_func_consequence = 'http://targetvalidation.org/sequence/' + most_severe_so_term.so_name
        else:
            self.gene_2_var_func_consequence = 'http://purl.obolibrary.org/obo/' + most_severe_so_term.accession.replace(':', '_')

        ref_list = list(set(traits_ref_list[trait_counter] + observed_refs_list + measure_set_refs_list))
        if len(ref_list) > 0:
            self.set_var_2_disease_literature(ref_list)
            # Arbitrarily select only one reference among all
            self.unique_reference = ref_list[0]
            self.top_level_literature = ref_list
        efo_list.sort()
        # Just (arbitrarily) adding one of the potentially multiple EFO terms because of schema constraints
        self.disease = efo_list[0]
        self.add_unique_association_field('phenotype', efo_list[0])

    @property
    def db_xref_url(self):
        if self['evidence']['gene2variant']['provenance_type']['database']['dbxref']['url'] \
                == self['evidence']['variant2disease']['provenance_type']['database']['dbxref']['url']:
            return self['evidence']['variant2disease']['provenance_type']['database']['dbxref']['url']
        else:
            raise Exception("db_xref_url attributes different")

    @db_xref_url.setter
    def db_xref_url(self, url):
        self['evidence']['gene2variant']['provenance_type']['database']['dbxref']['url'] = url
        self['evidence']['variant2disease']['provenance_type']['database']['dbxref']['url'] = url

    @property
    def url(self):
        if self['evidence']['gene2variant']['urls'][0]['url'] \
                == self['evidence']['variant2disease']['urls'][0]['url']:
            return self['evidence']['gene2variant']['urls'][0]['url']
        else:
            raise Exception("url attributes different")

    @url.setter
    def url(self, url):
        self['evidence']['gene2variant']['urls'][0]['url'] = url
        self['evidence']['variant2disease']['urls'][0]['url'] = url

    @property
    def gene_2_var_ev_codes(self):
        return self['evidence']['gene2variant']['evidence_codes']

    @gene_2_var_ev_codes.setter
    def gene_2_var_ev_codes(self, gene_2_var_ev_codes):
        self['evidence']['gene2variant']['evidence_codes'] = gene_2_var_ev_codes

    @property
    def gene_2_var_func_consequence(self):
        return self['evidence']['gene2variant']['functional_consequence']

    @gene_2_var_func_consequence.setter
    def gene_2_var_func_consequence(self, so_term):
        self['evidence']['gene2variant']['functional_consequence'] = so_term

    def set_var_2_disease_literature(self, ref_list):
        self['evidence']['variant2disease']['provenance_type']['literature'] = {'references': [{'lit_id': reference} for reference in ref_list]}

    @property
    def association(self):
        if self['evidence']['gene2variant']['is_associated'] \
                == self['evidence']['variant2disease']['is_associated']:
            return self['evidence']['gene2variant']['is_associated']
        else:
            raise Exception("association attributes different")

    @association.setter
    def association(self, is_associated):
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

    def _clear_variant(self):
        self['variant']['id'] = []

    def set_variant(self, id, type):
        self['variant']['id'].append(id)
        self['variant']['type'] = type

    @property
    def unique_reference(self):
        return self['evidence']['variant2disease']['unique_experiment_reference']

    @unique_reference.setter
    def unique_reference(self, reference):
        self['evidence']['variant2disease']['unique_experiment_reference'] = reference

    @property
    def date(self):
        if self['evidence']['gene2variant']['date_asserted'] == self['evidence']['variant2disease']['date_asserted']:
            return self['evidence']['gene2variant']['date_asserted']
        else:
            raise Exception("date attributes have different values")

    @date.setter
    def date(self, date_string):
        self['evidence']['gene2variant']['date_asserted'] = date_string
        self['evidence']['variant2disease']['date_asserted'] = date_string


class CTTVSomaticEvidenceString(CTTVEvidenceString):
    schema = json.loads(open(SOM_SCHEMA_FILE, 'r').read())

    def __init__(self, efo_list, clin_sig, clinvarRecord,
                                     ensembl_gene_id, ensembl_gene_id_uri, measure_set_refs_list,
                                     observed_refs_list, trait_counter, trait_refs_list,
                                     unrecognised_clin_sigs, consequenceType):
        a_dictionary = {'sourceID': 'eva_somatic',
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
        # CTTVEvidenceString.__init__(self,a_dictionary)

        super().__init__(a_dictionary)

        self.add_unique_association_field('gene', ensembl_gene_id)
        self.add_unique_association_field('clinvarAccession', clinvarRecord.acc)
        self.add_unique_association_field('alleleOrigin', 'somatic')
        try:
            self.set_target(ensembl_gene_id_uri, clin_sig_2_activity[clin_sig])
        except KeyError:
            unrecognised_clin_sigs.add(clin_sig)
            self.set_target(ensembl_gene_id_uri, 'http://identifiers.org/cttv.activity/unknown')

        self.date = clinvarRecord.date
        self.db_xref_url = 'http://identifiers.org/clinvar.record/' + clinvarRecord.acc
        self.url = 'http://www.ncbi.nlm.nih.gov/clinvar/' + clinvarRecord.acc
        self.association = (clin_sig != 'non-pathogenic' and clin_sig != 'probable-non-pathogenic' and clin_sig != 'likely benign' and clin_sig != 'benign')

        self.set_known_mutations(consequenceType)

        ref_list = list(set(trait_refs_list[trait_counter] + observed_refs_list + measure_set_refs_list))
        if len(ref_list) > 0:
            self.evidence_literature = ref_list
            self.top_level_literature = ref_list

        efo_list.sort()
        # Just (arbitrarily) adding one of the potentially multiple EFO terms because of schema constraints
        self.disease = efo_list[0]
        self.add_unique_association_field('phenotype', efo_list[0])

    @property
    def db_xref_url(self):
        return self['evidence']['provenance_type']['database']['dbxref']['url']

    @db_xref_url.setter
    def db_xref_url(self, url):
        self['evidence']['provenance_type']['database']['dbxref']['url'] = url

    @property
    def url(self):
        return self['evidence']['urls'][0]['url']

    @url.setter
    def url(self, url):
        self['evidence']['urls'][0]['url'] = url

    @property
    def evidence_literature(self):
        return self['evidence']['provenance_type']['literature']['references']

    @evidence_literature.setter
    def evidence_literature(self, ref_list):
        self['evidence']['provenance_type']['literature'] = {'references': [{'lit_id': reference} for reference in ref_list]}

    @property
    def association(self):
        return self['evidence']['is_associated']

    @association.setter
    def association(self, is_associated):
        self['evidence']['is_associated'] = is_associated

    def validate(self):
        jsonschema.validate(self, CTTVSomaticEvidenceString.schema, format_checker=jsonschema.FormatChecker())
        self.disease.is_obsolete()

    @property
    def date(self):
        return self['evidence']['date_asserted']

    @date.setter
    def date(self, date_string):
        self['evidence']['date_asserted'] = date_string

    def _clear_known_mutations(self):
        self['evidence']['known_mutations'] = []

    def add_known_mutation(self, new_functional_consequence, so_name):
        new_known_mutation = {'functional_consequence': new_functional_consequence, 'preferred_name': so_name}
        self['evidence']['known_mutations'].append(new_known_mutation)

    def set_known_mutations(self, consequence_type):
        for so_term in consequence_type.so_terms:
            so_name = so_term.so_name
            if so_term.accession:
                new_functional_consequence = "http://purl.obolibrary.org/obo/" + so_term.accession.replace(':', '_')
            else:
                new_functional_consequence = 'http://targetvalidation.org/sequence/' + so_term.so_name
            self.add_known_mutation(new_functional_consequence, so_name)
