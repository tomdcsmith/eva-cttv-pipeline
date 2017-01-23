import json
import gzip
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from datetime import datetime
import http.client
from collections import UserDict

from eva_cttv_pipeline import utilities


def get_rcv_to_rsnsv_mapping(variant_summary_file):

    rcv_to_rs = {}
    rcv_to_nsv = {}

    print('Loading mapping RCV->rs/nsv')
    fdr = utilities.open_file(variant_summary_file, "rt")
    fdr.readline()
    for line in fdr:
        parts = line.split('\t')
        rcv_list = parts[8].split(';')
        if parts[6] != '-' and parts[6] != '-1' and parts[6] != '':
            for rcv_id in rcv_list:
                rcv_to_rs[rcv_id] = 'rs' + parts[6]
        if parts[7] != '-' and parts[7] != '-1' and parts[7] != '':
            for rcv_id in rcv_list:
                rcv_to_nsv[rcv_id] = parts[7]
    fdr.close()
    print(' Done.')

    return rcv_to_rs, rcv_to_nsv


class ClinvarRecord(UserDict):
    """
    Class of which instances hold data on individual clinvar records. Subclass of UserDict rather
    than dict in order to use attributes
    """

    score_map = {
        "CLASSIFIED_BY_SINGLE_SUBMITTER": 1,
        "NOT_CLASSIFIED_BY_SUBMITTER": None,
        "CLASSIFIED_BY_MULTIPLE_SUBMITTERS": 2,
        "REVIEWED_BY_EXPERT_PANEL": 3,
        "REVIEWED_BY_PROFESSIONAL_SOCIETY": 4
    }

    def __init__(self, cellbase_dict, consequence_type_dict):
        UserDict.__init__(self, cellbase_dict)
        self.measures = [ClinvarRecordMeasure(measure_dict, self, consequence_type_dict)
                         for measure_dict in self.data['referenceClinVarAssertion']["measureSet"]["measure"]]

    @property
    def date(self):
        return datetime.fromtimestamp(
            self.data['referenceClinVarAssertion']['dateLastUpdated'] / 1000).isoformat()

    @property
    def score(self):
        return self.score_map[
            self.data['referenceClinVarAssertion']['clinicalSignificance']['reviewStatus']
        ]

    @property
    def accession(self):
        return self.data['referenceClinVarAssertion']['clinVarAccession']['acc']

    @property
    def traits(self):
        trait_list = []
        for trait in self.data['referenceClinVarAssertion']['traitSet']['trait']:
            trait_list.append([])
            for name in trait['name']:
                # First trait name in the list will always be the "Preferred" one
                if name['elementValue']['type'] == 'Preferred':
                    trait_list[-1] = [name['elementValue']['value']] + trait_list[-1]
                elif name['elementValue']['type'] in ["EFO URL", "EFO id", "EFO name"]:
                    continue  # if the trait name not originally from clinvar
                else:
                    trait_list[-1].append(name['elementValue']['value'])

        return trait_list

    @property
    def trait_pubmed_refs(self):
        pubmed_refs_list = []
        for trait in self.data['referenceClinVarAssertion']['traitSet']['trait']:
            pubmed_refs_list.append([])
            if 'citation' in trait:
                for citation in trait['citation']:
                    if ('id' in citation) and citation['id'] is not None:
                        for citation_id in citation['id']:
                            if citation_id['source'] == 'PubMed':
                                pubmed_refs_list[-1].append(int(citation_id['value']))

        return pubmed_refs_list

    @property
    def observed_pubmed_refs(self):
        pubmed_refs_list = []
        if 'observedIn' in self.data['referenceClinVarAssertion']:
            for observed_in in self.data['referenceClinVarAssertion']['observedIn']:
                for observed_data in observed_in['observedData']:
                    if 'citation' in observed_data:
                        for citation in observed_data['citation']:
                            if ('id' in citation) and citation['id'] is not None:
                                for citation_id in citation['id']:
                                    if citation_id['source'] == 'PubMed':
                                        pubmed_refs_list.append(int(citation_id['value']))
        return pubmed_refs_list

    @property
    def trait_refs_list(self):
        return [['http://europepmc.org/abstract/MED/' + str(ref) for ref in ref_list]
                for ref_list in self.trait_pubmed_refs]

    @property
    def observed_refs_list(self):
        return ['http://europepmc.org/abstract/MED/' + str(ref)
                for ref in self.observed_pubmed_refs]

    @property
    def clinical_significance(self):
        return \
            self.data['referenceClinVarAssertion']['clinicalSignificance']['description'].lower()

    @property
    def allele_origins(self):
        allele_origins = set()
        for clinvar_assetion_document in self.data['clinVarAssertion']:
            for observed_in_document in clinvar_assetion_document['observedIn']:
                allele_origins.add(observed_in_document['sample']['origin'].lower())

        return list(allele_origins)


class ClinvarRecordMeasure(UserDict):

    def __init__(self, clinvar_measure_dict, clinvar_record, consequence_type_dict):
        UserDict.__init__(self, clinvar_measure_dict)
        self.clinvar_record = clinvar_record
        self.consequence_type = self.__get_main_consequence_types(consequence_type_dict)

    @property
    def rs_id(self):
        for xref in self.data["xref"]:
            if xref["db"].lower() == "dbsnp":
                return "rs{}".format(xref["id"])
        return None

    @property
    def nsv_id(self):
        for xref in self.data["xref"]:
            if xref["db"].lower() == "dbvar":
                return xref["id"]
        return None

    def __get_main_consequence_types(self, consequence_type_dict):

        if self.rs_id is not None and self.rs_id in consequence_type_dict:
            return consequence_type_dict[self.rs_id]
        elif self.nsv_id is not None and self.nsv_id in consequence_type_dict:
            return consequence_type_dict[self.nsv_id]
        elif self.clinvar_record.accession in consequence_type_dict:
            return consequence_type_dict[self.clinvar_record.accession]  # todo change this depending upon OT gene mapping file
        else:
            return None

    @property
    def hgvs(self):
        hgvs_list = []
        for attribute_set in self.data['attributeSet']:
            if attribute_set['attribute']['type'].startswith('HGVS'):
                hgvs_list.append(attribute_set['attribute']['value'])

        return hgvs_list

    @property
    def variant_type(self):
        return self.data['type']

    @property
    def measure_set_pubmed_refs(self):
        pubmed_refs_list = []
        if 'citation' in self.data:
            for citation in self.data['citation']:
                if 'id' in citation and citation['id'] is not None:
                    for citation_id in citation['id']:
                        if citation_id['source'] == 'PubMed':
                            pubmed_refs_list.append(int(citation_id['value']))
        return pubmed_refs_list

    @property
    def measure_set_refs_list(self):
        return ['http://europepmc.org/abstract/MED/' + str(ref)
                for ref in self.measure_set_pubmed_refs]

    @property
    def chr(self):
        return self.sequence_location_helper("chr")

    @property
    def start(self):
        return self.sequence_location_helper("start")

    @property
    def stop(self):
        return self.sequence_location_helper("stop")

    @property
    def ref(self):
        return self.sequence_location_helper("referenceAllele")

    @property
    def alt(self):
        return self.sequence_location_helper("alternateAllele")

    def sequence_location_helper(self, attr):
        for sequence_location in self.data["sequenceLocation"]:
            if sequence_location["assembly"].lower() == "grch38":
                if attr in sequence_location:
                    return sequence_location[attr]
                else:
                    return None

