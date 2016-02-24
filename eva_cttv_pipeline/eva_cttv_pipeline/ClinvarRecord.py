import http.client
import json
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from datetime import datetime
import xlrd
import eva_cttv_pipeline.utilities as utilities
import eva_cttv_pipeline.config as config
from eva_cttv_pipeline import ConsequenceType

__author__ = 'Javier Lopez: javild@gmail.com'

CTMAPPINGFILE = utilities.get_resource_file(__package__, config.con_type_file)
RCVTORSFILE = utilities.get_resource_file(__package__, config.variant_summary)


def _process_con_type_file_xls():

    oneRsMultipleGenes = set()
    consequenceTypeDict = {}

    CTMappingReadBook = xlrd.open_workbook(CTMAPPINGFILE, formatting_info=True)
    CTMappingReadSheet = CTMappingReadBook.sheet_by_index(0)
    for i in range(1, CTMappingReadSheet.nrows):
        if CTMappingReadSheet.cell_value(rowx=i, colx=2) != 'Not found':

            rs_id = CTMappingReadSheet.cell_value(rowx=i, colx=0)
            ensemblGeneId = CTMappingReadSheet.cell_value(rowx=i, colx=2)
            so_term = CTMappingReadSheet.cell_value(rowx=i, colx=1)

            if rs_id in consequenceTypeDict:
                if ensemblGeneId != consequenceTypeDict[rs_id].getEnsemblGeneId():
                    print('WARNING (ClinvarRecord.py): different genes and annotations found for a given gene.')
                    print(' Variant id: ' + rs_id + ', ENSG: ' + so_term + ', ENSG: ' + consequenceTypeDict[rs_id].getEnsemblGeneId())
                    print('Skipping')
                    oneRsMultipleGenes.add(rs_id)
                else:
                    consequenceTypeDict[rs_id].addSoTerm(so_term)
            else:
                consequenceTypeDict[rs_id] = ConsequenceType.ConsequenceType(ensemblGeneId, [so_term])

    return oneRsMultipleGenes, consequenceTypeDict


def _process_gene(consequenceTypeDict, rs_id, ensembl_gene_id, so_term):
    if rs_id in consequenceTypeDict:
        consequenceTypeDict[rs_id].add_ensembl_gene_id(ensembl_gene_id)
        consequenceTypeDict[rs_id].addSoTerm(so_term)
    else:
        consequenceTypeDict[rs_id] = ConsequenceType.ConsequenceType([ensembl_gene_id], [so_term])


def _process_con_type_file_tsv():
    oneRsMultipleGenes = set()
    consequenceTypeDict = {}

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
                    _process_gene(consequenceTypeDict, rs_id, ensembl_gene_id, so_term)
            else:
                _process_gene(consequenceTypeDict, rs_id, ensembl_gene_id, so_term)

    return oneRsMultipleGenes, consequenceTypeDict


def _process_con_type_file():

    print('Loading mappintg rs->ENSG/SOterms')

    if config.con_type_file.endswith(".xls"):
        oneRsMultipleGenes, consequenceTypeDict = _process_con_type_file_xls()
    else:
        oneRsMultipleGenes, consequenceTypeDict = _process_con_type_file_tsv()

    print(str(len(consequenceTypeDict)) + ' rs->ENSG/SOterms mappings loaded')
    print(str(len(oneRsMultipleGenes)) + ' rsIds with multiple gene associations')
    print('Done.')

    return consequenceTypeDict


class ClinvarRecord(dict):
    cachedSymbol2Ensembl = {}

    consequenceTypeDict = _process_con_type_file()

    print('Loading mapping RCV->rs/nsv')
    rcvToRs = {}
    rcvToNsv = {}
    fdr = open(RCVTORSFILE, "r")
    fdr.readline()
    for line in fdr:
        parts = line.split('\t')
        rcvList = parts[8].split(';')
        if parts[6] != '-' and parts[6] != '-1' and parts[6] != '':
            for rcvId in rcvList:
                rcvToRs[rcvId] = 'rs' + parts[6]
        if parts[7] != '-' and parts[7] != '-1' and parts[7] != '':
            for rcvId in rcvList:
                rcvToNsv[rcvId] = parts[7]
    fdr.close()
    print(' Done.')

    def __init__(self, aDictionary=None):
        if aDictionary is None:
            dict.__init__(self)
        else:
            dict.__init__(self, aDictionary)

        self.scoreMap = {
            "CLASSIFIED_BY_SINGLE_SUBMITTER": 1,
            "NOT_CLASSIFIED_BY_SUBMITTER": None,
            "CLASSIFIED_BY_MULTIPLE_SUBMITTERS": 2,
            "REVIEWED_BY_EXPERT_PANEL": 3,
            "REVIEWED_BY_PROFESSIONAL_SOCIETY": 4
        }

    def getGeneId(self):
        j = 0
        measure = self['referenceClinVarAssertion']['measureSet']['measure']
        found = False
        while j < len(measure) and not found:
            attributeSet = measure[j]['attributeSet']
            i = 0
            while i < len(attributeSet) and not attributeSet[i]['attribute']['type'].startswith('HGVS'):
                i += 1
            found = (i < len(attributeSet))
            j += 1

        if found:
            return attributeSet[i]['attribute']['value'][:9]
        else:
            return self['referenceClinVarAssertion']['measureSet']['measure'][0]['name'][0]['elementValue']['value']

    def getEnsemblId(self):
        global ensemblJson
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
                                return ClinvarRecord.cachedSymbol2Ensembl[symbol[i]['elementValue']['value']]
                            except KeyError:
                                notSolved = True
                                while notSolved:
                                    try:
                                        ensemblJson = json.loads(urllib.request.urlopen(
                                            'http://rest.ensembl.org/lookup/symbol/homo_sapiens/' +
                                            symbol[i]['elementValue'][
                                                'value'] + '?content-type=application/json').read())
                                        notSolved = False
                                    except urllib.error.HTTPError as e:
                                        if e.code == 400:
                                            print('WARNING: Bad request code returned from ENSEMBL rest.')
                                            print(' ClinVar accession: ' + self.getAcc())
                                            print(' Gene symbol: ' + symbol[i]['elementValue']['value'])
                                            print(' Error: ')
                                            print(e)
                                            print(' Returning None.')
                                            ClinvarRecord.cachedSymbol2Ensembl[
                                                symbol[i]['elementValue']['value']] = None
                                            return None
                                        else:
                                            time.sleep(0.05)
                                    except http.client.BadStatusLine:
                                        time.sleep(3)

                                if len(ensemblJson) > 0:
                                    if (len(ensemblJson) == 14) and (ensemblJson['object_type'] == 'Gene'):
                                        ClinvarRecord.cachedSymbol2Ensembl[symbol[i]['elementValue']['value']] = \
                                            ensemblJson['id']
                                        return ensemblJson['id']
                                    else:
                                        print("WARNING at ClinvarRecord.py: ENSEMBL's REST API returned an unexpected json object")
                                        print(" Queried gene symbol: " + symbol[i]['elementValue']['value'])
                                        print(" ENSEMBL's API response:")
                                        print(ensemblJson)
                                        print("Exiting.")
                                        sys.exit(1)
                                else:
                                    ClinvarRecord.cachedSymbol2Ensembl[symbol[i]['elementValue']['value']] = None
                            i += 1
                    l += 1
            j += 1

        return None

    def getDate(self):
        return datetime.fromtimestamp(self['referenceClinVarAssertion']['dateLastUpdated'] / 1000).isoformat()

    def getScore(self):
        return self.scoreMap[self['referenceClinVarAssertion']['clinicalSignificance']['reviewStatus']]

    def getAcc(self):
        return self['referenceClinVarAssertion']['clinVarAccession']['acc']

    def getTraits(self):
        traitList = []
        for traitRecord in self['referenceClinVarAssertion']['traitSet']['trait']:
            traitList.append([])
            for nameRecord in traitRecord['name']:
                if (nameRecord['elementValue'][
                        'type'] == 'Preferred'):  # First trait name in the list will always be the "Preferred" one
                    traitList[-1] = [nameRecord['elementValue']['value']] + traitList[-1]
                else:
                    traitList[-1].append(nameRecord['elementValue']['value'])

        return traitList

    def getTraitPubmedrefs(self):
        pubmedrefsList = []
        for traitRecord in self['referenceClinVarAssertion']['traitSet']['trait']:
            pubmedrefsList.append([])
            if 'citation' in traitRecord:
                for citationRecord in traitRecord['citation']:
                    if ('id' in citationRecord) and citationRecord['id'] is not None:
                        if citationRecord['id']['source'] == 'PubMed':
                            pubmedrefsList[-1].append(int(citationRecord['id']['value']))

        return pubmedrefsList

    def getObservedPubmedrefs(self):
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

    def getMeasureSetPubmedrefs(self):
        pubmedrefsList = []
        for measureRecord in self['referenceClinVarAssertion']['measureSet']['measure']:
            if 'citation' in measureRecord:
                for citationRecord in measureRecord['citation']:
                    if ('id' in citationRecord) and citationRecord['id'] is not None:
                        if citationRecord['id']['source'] == 'PubMed':
                            pubmedrefsList.append(int(citationRecord['id']['value']))
        return pubmedrefsList

    def getHGVS(self):
        hgvsList = []
        for measureRecord in self['referenceClinVarAssertion']['measureSet']['measure']:
            for attributeSetRecord in measureRecord['attributeSet']:
                if attributeSetRecord['attribute']['type'].startswith('HGVS'):
                    hgvsList.append(attributeSetRecord['attribute']['value'])

        return hgvsList

    def getClinicalSignificance(self):
        return self['referenceClinVarAssertion']['clinicalSignificance']['description']

    def getRs(self):
        try:
            return ClinvarRecord.rcvToRs[self.getAcc()]
        except KeyError:
            return None

    def getNsv(self):
        try:
            return ClinvarRecord.rcvToNsv[self.getAcc()]
        except KeyError:
            return None

    def getMainConsequenceTypes(self):
        newRsId = self.getRs()
        if newRsId is not None and (newRsId in ClinvarRecord.consequenceTypeDict):
            return ClinvarRecord.consequenceTypeDict[newRsId]
        else:
            return None

    def getVariantType(self):
        return self['referenceClinVarAssertion']['measureSet']['measure'][0]['type']

    def getAlleleOrigins(self):
        alleleOrigins = set()
        for clinvarAssertionDocument in self['clinVarAssertion']:
            for observedInDocument in clinvarAssertionDocument['observedIn']:
                alleleOrigins.add(observedInDocument['sample']['origin'])

        return list(alleleOrigins)
