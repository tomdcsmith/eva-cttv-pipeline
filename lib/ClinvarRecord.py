__author__ = 'Javier Lopez: javild@gmail.com'


from datetime import datetime
import json
import urllib2
import httplib
import time
import sys
import xlrd
import os

import ConsequenceType


CTMAPPINGFILE = os.path.dirname(__file__)+"/resources/eva_cttv_snp2gene_mapping_20150512.xls"
RCVTORSFILE = os.path.dirname(__file__)+"/resources/variant_summary.txt"


class ClinvarRecord(dict):

    cachedSymbol2Ensembl = {}
    oneRsMultipleGenes = set()

    print 'Loading mappintg rs->ENSG/SOterms '
    CTMappingReadBook = xlrd.open_workbook(CTMAPPINGFILE, formatting_info=True)
    CTMappingReadSheet = CTMappingReadBook.sheet_by_index(0)
    consequenceTypeDict = {}
    for i in range(1, CTMappingReadSheet.nrows):
        if CTMappingReadSheet.cell_value(rowx=i, colx=2) != 'Not found':
            ensemblGeneId = CTMappingReadSheet.cell_value(rowx=i, colx=2)
            if CTMappingReadSheet.cell_value(rowx=i, colx=0) in consequenceTypeDict:
                if ensemblGeneId != consequenceTypeDict[CTMappingReadSheet.cell_value(rowx=i, colx=0)].getEnsemblGeneId():
                    print 'WARNING (ClinvarRecord.py): different genes and annotations found for a given gene.'
                    print ' Variant id: '+CTMappingReadSheet.cell_value(rowx=i, colx=0)+', ENSG: '+CTMappingReadSheet.cell_value(rowx=i, colx=1)+', ENSG: '+consequenceTypeDict[CTMappingReadSheet.cell_value(rowx=i, colx=0)].getEnsemblGeneId()
                    print 'Skipping'
                    oneRsMultipleGenes.add(CTMappingReadSheet.cell_value(rowx=i, colx=0))
                else:
                    consequenceTypeDict[CTMappingReadSheet.cell_value(rowx=i, colx=0)].addSoTerm(CTMappingReadSheet.cell_value(rowx=i, colx=1))
            else:
                consequenceTypeDict[CTMappingReadSheet.cell_value(rowx=i, colx=0)] = ConsequenceType.ConsequenceType(CTMappingReadSheet.cell_value(rowx=i, colx=2), [CTMappingReadSheet.cell_value(rowx=i, colx=1)])
    print str(len(consequenceTypeDict)) + ' rs->ENSG/SOterms mappings loaded'
    print str(len(oneRsMultipleGenes))+' rsIds with multiple gene associations'
    print ' Done.'

    print 'Loading mapping RCV->rs/nsv'
    rcvToRs = {}
    rcvToNsv = {}
    fdr = file(RCVTORSFILE)
    fdr.readline()
    for line in fdr:
        parts = line.split('\t')
        rcvList = parts[8].split(';')
        if parts[6] != '-' and parts[6] != '-1' and parts[6] != '':
            for rcvId in rcvList:
                rcvToRs[rcvId] = 'rs'+parts[6]
        if parts[7] != '-' and parts[7] != '-1' and parts[7] != '':
            for rcvId in rcvList:
                rcvToNsv[rcvId] = parts[7]
    fdr.close()
    print ' Done.'

    def __init__(self, aDictionary=None):
        if aDictionary is None:
            dict.__init__(self)
        else:
            dict.__init__(self, aDictionary)

        self.scoreMap = {
                            "CLASSIFIED_BY_SINGLE_SUBMITTER" : 1,
                            "NOT_CLASSIFIED_BY_SUBMITTER" : None,
                            "CLASSIFIED_BY_MULTIPLE_SUBMITTERS" : 2,
                            "REVIEWED_BY_EXPERT_PANEL" : 3,
                            "REVIEWED_BY_PROFESSIONAL_SOCIETY" : 4
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
            found = (i<len(attributeSet))
            j += 1

        if(found):
            return attributeSet[i]['attribute']['value'][:9]
        else:
            return self['referenceClinVarAssertion']['measureSet']['measure'][0]['name'][0]['elementValue']['value']

    def getEnsemblId(self):
        j = 0
        measure = self['referenceClinVarAssertion']['measureSet']['measure']
        while(j<len(measure)):
            if('measureRelationship' in measure[j]):
                measureRelationship = measure[j]['measureRelationship']
                l=0
                while(l<len(measureRelationship)):
                    if('symbol' in measureRelationship[l]):
                        symbol = measureRelationship[l]['symbol']
                        i = 0
                        while(i<len(symbol)):
                            try:
                                return ClinvarRecord.cachedSymbol2Ensembl[symbol[i]['elementValue']['value']]
                            except KeyError:
                                notSolved=True
                                while(notSolved):
                                    try:
                                        ensemblJson = json.loads(urllib2.urlopen('http://rest.ensembl.org/lookup/symbol/homo_sapiens/'+symbol[i]['elementValue']['value']+'?content-type=application/json').read())
                                        notSolved = False
                                    except urllib2.HTTPError, e:
                                        if(e.code == 400):
                                            print 'WARNING: Bad request code returned from ENSEMBL rest.'
                                            print ' ClinVar accession: '+self.getAcc()
                                            print ' Gene symbol: '+symbol[i]['elementValue']['value']
                                            print ' Error: '
                                            print e
                                            print ' Returning None.'
                                            ClinvarRecord.cachedSymbol2Ensembl[symbol[i]['elementValue']['value']] = None
                                            return None
                                        else:
                                            time.sleep(0.05)
                                    except httplib.BadStatusLine:
                                        time.sleep(3)

                                if(len(ensemblJson)>0):
                                    if((len(ensemblJson)==14) and (ensemblJson['object_type']=='Gene')):
                                        ClinvarRecord.cachedSymbol2Ensembl[symbol[i]['elementValue']['value']] = ensemblJson['id']
                                        return ensemblJson['id']
                                    else:
                                        print "WARNING at ClinvarRecord.py: ENSEMBL's REST API returned an unexpected json object"
                                        print " Queried gene symbol: "+symbol[i]['elementValue']['value']
                                        print " ENSEMBL's API response:"
                                        print ensemblJson
                                        print "Exiting."
                                        sys.exit(1)
                                else:
                                    ClinvarRecord.cachedSymbol2Ensembl[symbol[i]['elementValue']['value']] = None
                            i+=1
                    l += 1
            j+=1

        return None

    def getDate(self):
        return datetime.fromtimestamp(self['referenceClinVarAssertion']['dateLastUpdated']/1000).isoformat()

    def getScore(self):
        return self.scoreMap[self['referenceClinVarAssertion']['clinicalSignificance']['reviewStatus']]

    def getAcc(self):
        return self['referenceClinVarAssertion']['clinVarAccession']['acc']

    def getTraits(self):
        traitList = []
        for traitRecord in self['referenceClinVarAssertion']['traitSet']['trait']:
            traitList.append([])
            for nameRecord in traitRecord['name']:
                if(nameRecord['elementValue']['type'] == 'Preferred') :  # First trait name in the list will always be the "Preferred" one
                    traitList[-1] = [nameRecord['elementValue']['value']] + traitList[-1]
                else:
                    traitList[-1].append(nameRecord['elementValue']['value'])

        return traitList

    def getTraitPubmedrefs(self):
        pubmedrefsList = []
        for traitRecord in self['referenceClinVarAssertion']['traitSet']['trait']:
            pubmedrefsList.append([])
            if('citation' in traitRecord):
                for citationRecord in traitRecord['citation']:
                    if(('id' in citationRecord) and citationRecord['id']!=None):
                        if(citationRecord['id']['source'] == 'PubMed'):
                            pubmedrefsList[-1].append(int(citationRecord['id']['value']))

        return pubmedrefsList

    def getObservedPubmedrefs(self):
        pubmedrefsList = []
        if('observedIn' in self['referenceClinVarAssertion']):
            for observedInRecord in self['referenceClinVarAssertion']['observedIn']:
                for observedDataRecord in observedInRecord['observedData']:
                    if('citation' in observedDataRecord):
                        for citationRecord in observedDataRecord['citation']:
                            if(('id' in citationRecord) and citationRecord['id']!=None):
                                if(citationRecord['id']['source'] == 'PubMed'):
                                    pubmedrefsList.append(int(citationRecord['id']['value']))
        return pubmedrefsList

    def getMeasureSetPubmedrefs(self):
        pubmedrefsList = []
        for measureRecord in self['referenceClinVarAssertion']['measureSet']['measure']:
            if('citation' in measureRecord):
                for citationRecord in measureRecord['citation']:
                    if(('id' in citationRecord) and citationRecord['id']!=None):
                        if(citationRecord['id']['source'] == 'PubMed'):
                            pubmedrefsList.append(int(citationRecord['id']['value']))
        return pubmedrefsList

    def getHGVS(self):
        hgvsList = []
        for measureRecord in self['referenceClinVarAssertion']['measureSet']['measure']:
            for attributeSetRecord in measureRecord['attributeSet']:
                if(attributeSetRecord['attribute']['type'].startswith('HGVS')):
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
        if(newRsId != None and (newRsId in ClinvarRecord.consequenceTypeDict)):
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









