import json
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from datetime import datetime
import http.client

import xlrd

import eva_cttv_pipeline.config as config
import eva_cttv_pipeline.utilities as utilities
from eva_cttv_pipeline import consequence_type

__author__ = 'Javier Lopez: javild@gmail.com'

CTMAPPINGFILE = utilities.get_resource_file(__package__, config.con_type_file)
RCVTORSFILE = utilities.get_resource_file(__package__, config.variant_summary)


def _process_con_type_file_xls():

    one_rs_multiple_genes = set()
    consequence_type_dict = {}

    ct_mapping_read_book = xlrd.open_workbook(CTMAPPINGFILE, formatting_info=True)
    ct_mapping_read_sheet = ct_mapping_read_book.sheet_by_index(0)
    for i in range(1, ct_mapping_read_sheet.nrows):
        if ct_mapping_read_sheet.cell_value(rowx=i, colx=2) != 'Not found':

            rs_id = ct_mapping_read_sheet.cell_value(rowx=i, colx=0)
            ensembl_gene_id = ct_mapping_read_sheet.cell_value(rowx=i, colx=2)
            so_term = ct_mapping_read_sheet.cell_value(rowx=i, colx=1)

            if rs_id in consequence_type_dict:
                if ensembl_gene_id not in consequence_type_dict[rs_id].get_ensembl_gene_ids():
                    print('WARNING (clinvar_record.py): different genes and annotations found for a given gene.')
                    print(' Variant id: ' + rs_id + ', ENSG: ' + so_term + ', ENSG: ' + str(consequence_type_dict[rs_id].get_ensembl_gene_ids()))
                    print('Skipping')
                    one_rs_multiple_genes.add(rs_id)
                else:
                    consequence_type_dict[rs_id].add_so_term(so_term)
            else:
                consequence_type_dict[rs_id] = consequence_type.ConsequenceType(ensembl_gene_id, [so_term])

    return one_rs_multiple_genes, consequence_type_dict


def _process_gene(consequence_type_dict, rs_id, ensembl_gene_id, so_term):
    if rs_id in consequence_type_dict:
        consequence_type_dict[rs_id].add_ensembl_gene_id(ensembl_gene_id)
        consequence_type_dict[rs_id].add_so_term(so_term)
    else:
        consequence_type_dict[rs_id] = consequence_type.ConsequenceType([ensembl_gene_id], [so_term])


def _process_con_type_file_tsv():
    one_rs_multiple_genes = set()
    consequence_type_dict = {}

    with open(CTMAPPINGFILE, "rt") as f:
        for line in f:
            line = line.rstrip()
            line_list = line.split("\t")

            rs_id = line_list[0]
            ensembl_gene_id = line_list[2]
            if not ensembl_gene_id or rs_id == "rs":
                continue
            so_term = line_list[4]

            if "," in ensembl_gene_id:
                ensembl_gene_ids = ensembl_gene_id.split(",")
                for ensembl_gene_id in ensembl_gene_ids:
                    _process_gene(consequence_type_dict, rs_id, ensembl_gene_id, so_term)
            else:
                _process_gene(consequence_type_dict, rs_id, ensembl_gene_id, so_term)

    return one_rs_multiple_genes, consequence_type_dict


def _process_con_type_file():

    print('Loading mappintg rs->ENSG/SOterms')

    if config.con_type_file.endswith(".xls"):
        one_rs_multiple_genes, consequence_type_dict = _process_con_type_file_xls()
    else:
        one_rs_multiple_genes, consequence_type_dict = _process_con_type_file_tsv()

    print(str(len(consequence_type_dict)) + ' rs->ENSG/SOterms mappings loaded')
    print(str(len(one_rs_multiple_genes)) + ' rsIds with multiple gene associations')
    print('Done.')

    return consequence_type_dict


class ClinvarRecord(dict):
    cached_symbol_2_ensembl = {}

    consequence_type_dict = _process_con_type_file()

    print('Loading mapping RCV->rs/nsv')
    rcv_to_rs = {}
    rcv_to_nsv = {}
    fdr = open(RCVTORSFILE, "r")
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

    def __init__(self, a_dictionary=None):
        if a_dictionary is None:
            dict.__init__(self)
        else:
            dict.__init__(self, a_dictionary)

        self.score_map = {
            "CLASSIFIED_BY_SINGLE_SUBMITTER": 1,
            "NOT_CLASSIFIED_BY_SUBMITTER": None,
            "CLASSIFIED_BY_MULTIPLE_SUBMITTERS": 2,
            "REVIEWED_BY_EXPERT_PANEL": 3,
            "REVIEWED_BY_PROFESSIONAL_SOCIETY": 4
        }

    def get_gene_id(self):
        j = 0
        measure = self['referenceClinVarAssertion']['measureSet']['measure']
        found = False
        while j < len(measure) and not found:
            attribute_set = measure[j]['attributeSet']
            i = 0
            while i < len(attribute_set) and not attribute_set[i]['attribute']['type'].startswith('HGVS'):
                i += 1
            found = (i < len(attribute_set))
            j += 1

        if found:
            return attribute_set[i]['attribute']['value'][:9]
        else:
            return self['referenceClinVarAssertion']['measureSet']['measure'][0]['name'][0]['elementValue']['value']

    def get_ensembl_id(self):
        global ensembl_json
        j = 0
        measure = self['referenceClinVarAssertion']['measureSet']['measure']
        while j < len(measure):
            if 'measureRelationship' in measure[j]:
                measureRelationship = measure[j]['measureRelationship']
                l = 0
                while l < len(measureRelationship):
                    if 'symbol' in measureRelationship[l]:
                        symbol = measureRelationship[l]['symbol']
                        i = 0
                        while i < len(symbol):
                            try:
                                return ClinvarRecord.cached_symbol_2_ensembl[symbol[i]['elementValue']['value']]
                            except KeyError:
                                not_solved = True
                                while not_solved:
                                    try:
                                        ensembl_json = json.loads(urllib.request.urlopen(
                                            'http://rest.ensembl.org/lookup/symbol/homo_sapiens/' +
                                            symbol[i]['elementValue'][
                                                'value'] + '?content-type=application/json').read())
                                        not_solved = False
                                    except urllib.error.HTTPError as e:
                                        if e.code == 400:
                                            print('WARNING: Bad request code returned from ENSEMBL rest.')
                                            print(' ClinVar accession: ' + self.get_acc())
                                            print(' Gene symbol: ' + symbol[i]['elementValue']['value'])
                                            print(' Error: ')
                                            print(e)
                                            print(' Returning None.')
                                            ClinvarRecord.cached_symbol_2_ensembl[
                                                symbol[i]['elementValue']['value']] = None
                                            return None
                                        else:
                                            time.sleep(0.05)
                                    except http.client.BadStatusLine:
                                        time.sleep(3)

                                if len(ensembl_json) > 0:
                                    if (len(ensembl_json) == 14) and (ensembl_json['object_type'] == 'Gene'):
                                        ClinvarRecord.cached_symbol_2_ensembl[symbol[i]['elementValue']['value']] = \
                                            ensembl_json['id']
                                        return ensembl_json['id']
                                    else:
                                        print("WARNING at clinvar_record.py: ENSEMBL's REST API returned an unexpected json object")
                                        print(" Queried gene symbol: " + symbol[i]['elementValue']['value'])
                                        print(" ENSEMBL's API response:")
                                        print(ensembl_json)
                                        print("Exiting.")
                                        sys.exit(1)
                                else:
                                    ClinvarRecord.cached_symbol_2_ensembl[symbol[i]['elementValue']['value']] = None
                            i += 1
                    l += 1
            j += 1

        return None

    def get_date(self):
        return datetime.fromtimestamp(self['referenceClinVarAssertion']['dateLastUpdated'] / 1000).isoformat()

    def get_score(self):
        return self.score_map[self['referenceClinVarAssertion']['clinicalSignificance']['reviewStatus']]

    def get_acc(self):
        return self['referenceClinVarAssertion']['clinVarAccession']['acc']

    def get_traits(self):
        trait_list = []
        for trait_record in self['referenceClinVarAssertion']['traitSet']['trait']:
            trait_list.append([])
            for name_record in trait_record['name']:
                if (name_record['elementValue'][
                        'type'] == 'Preferred'):  # First trait name in the list will always be the "Preferred" one
                    trait_list[-1] = [name_record['elementValue']['value']] + trait_list[-1]
                else:
                    trait_list[-1].append(name_record['elementValue']['value'])

        return trait_list

    def get_trait_pubmed_refs(self):
        pubmed_refs_list = []
        for trait_record in self['referenceClinVarAssertion']['traitSet']['trait']:
            pubmed_refs_list.append([])
            if 'citation' in trait_record:
                for citation_record in trait_record['citation']:
                    if ('id' in citation_record) and citation_record['id'] is not None:
                        if citation_record['id']['source'] == 'PubMed':
                            pubmed_refs_list[-1].append(int(citation_record['id']['value']))

        return pubmed_refs_list

    def get_observed_pubmed_refs(self):
        pubmedrefsList = []
        if 'observedIn' in self['referenceClinVarAssertion']:
            for observedInRecord in self['referenceClinVarAssertion']['observedIn']:
                for observedDataRecord in observedInRecord['observedData']:
                    if 'citation' in observedDataRecord:
                        for citationRecord in observedDataRecord['citation']:
                            if ('id' in citationRecord) and citationRecord['id'] is not None:
                                if citationRecord['id']['source'] == 'PubMed':
                                    pubmedrefsList.append(int(citationRecord['id']['value']))
        return pubmedrefsList

    def get_measure_set_pubmed_refs(self):
        pubmed_refs_list = []
        for measure_record in self['referenceClinVarAssertion']['measureSet']['measure']:
            if 'citation' in measure_record:
                for ciration_record in measure_record['citation']:
                    if ('id' in ciration_record) and ciration_record['id'] is not None:
                        if ciration_record['id']['source'] == 'PubMed':
                            pubmed_refs_list.append(int(ciration_record['id']['value']))
        return pubmed_refs_list

    def get_hgvs(self):
        hgvs_list = []
        for measure_record in self['referenceClinVarAssertion']['measureSet']['measure']:
            for attribute_set_record in measure_record['attributeSet']:
                if attribute_set_record['attribute']['type'].startswith('HGVS'):
                    hgvs_list.append(attribute_set_record['attribute']['value'])

        return hgvs_list

    def get_clinical_significance(self):
        return self['referenceClinVarAssertion']['clinicalSignificance']['description']

    def get_rs(self):
        try:
            return ClinvarRecord.rcv_to_rs[self.get_acc()]
        except KeyError:
            return None

    def get_nsv(self):
        try:
            return ClinvarRecord.rcv_to_nsv[self.get_acc()]
        except KeyError:
            return None

    def get_main_consequence_types(self):
        new_rs_id = self.get_rs()
        if new_rs_id is not None and (new_rs_id in ClinvarRecord.consequence_type_dict):
            return ClinvarRecord.consequence_type_dict[new_rs_id]
        else:
            return None

    def get_variant_type(self):
        return self['referenceClinVarAssertion']['measureSet']['measure'][0]['type']

    def get_allele_origins(self):
        allele_origins = set()
        for clinvar_assetion_document in self['clinVarAssertion']:
            for observed_in_document in clinvar_assetion_document['observedIn']:
                allele_origins.add(observed_in_document['sample']['origin'])

        return list(allele_origins)
