from eva_cttv_pipeline import EFOTerm
import json
import eva_cttv_pipeline.utilities as utilities
import jsonschema
import eva_cttv_pipeline.config as config

__author__ = 'Javier Lopez: javild@gmail.com'


utilities.check_for_local_schema()
GEN_SCHEMA_FILE = utilities.get_resource_file(__package__, config.local_schema + "/src/genetics.json")
SOM_SCHEMA_FILE = utilities.get_resource_file(__package__, config.local_schema + "/src/literature_curated.json")


class CTTVEvidenceString(dict):
    def __init__(self, a_dictionary):
        dict.__init__(self, a_dictionary)

    def addUniqueAssociationField(self, key, value):
        self['unique_association_fields'][key] = value

    def setUniqueAssociationField(self, key, value):
        self['unique_association_fields'][key] = value

    def setTarget(self, id, activity):
        self['target']['id'].append(id)
        self['target']['activity'] = activity

    def setPubmedrefs(self, pubmedList):
        if pubmedList is not None and len(pubmedList)>0:
            if 'literature' in self['evidence']['provenance_type']:
                self['evidence']['provenance_type']['literature']['pubmed_refs'] = pubmedList
            else:
                self['evidence']['provenance_type']['literature'] = {'pubmed_refs' : pubmedList}

    def setDisease(self, id):
        self['disease']['id'].append(id)

    def get_disease(self):
        return EFOTerm.EFOTerm(self['disease']['id'][0])

    def setEvidenceCodes(self, evidenceCodeList):
        self['evidence']['evidence_codes'] = evidenceCodeList

    def getEvidenceCodes(self):
        return self['evidence']['evidence_codes']

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
        # try:
        jsonschema.validate(self, CTTVGeneticsEvidenceString.schema, format_checker=jsonschema.FormatChecker())
        # except Exception as e:
        #     print(str(self))
        #     print(e)
        #     traceback.print_stack()
        #     sys.exit()
        self.get_disease().isObsolete()

    def setVariant(self, id, type):
        self['variant']['id'].append(id)
        self['variant']['type'] = type

    def setUniqueReference(self, reference):
        self['evidence']['variant2disease']['unique_experiment_reference'] = reference

    def setDate(self, dateString):
        self['evidence']['gene2variant']['date_asserted'] = dateString
        self['evidence']['variant2disease']['date_asserted'] = dateString


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
        self.get_disease().isObsolete()

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
