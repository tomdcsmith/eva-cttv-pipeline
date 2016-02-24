#!/usr/bin/python

import argparse
import codecs
import json
import os
import sys
import urllib.error
import urllib.parse
import urllib.request
import jsonschema
import xlrd
from eva_cttv_pipeline import EFOTerm, clinvar_record, evidence_strings, utilities

__author__ = 'Javier Lopez: javild@gmail.com'

# sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.realpath(__file__))) + "/lib")
# print(os.path.dirname(os.path.dirname(__file__)) + "/lib")
# print(os.path.dirname(os.path.dirname(os.path.realpath(__file__))) + "/lib")

# sys.path.append(os.path.dirname(os.path.dirname(__file__)) + "/lib")  # Adds eva_cttv_pipeline root dir to the
#  PYTHONPATH

BATCH_SIZE = 200
# HOST = 'localhost:8080'
HOST = 'www.ebi.ac.uk'
EFOMAPPINGFILE = os.path.dirname(clinvar_record.__file__) + "/resources/ClinVar_Traits_EFO_090915.xls"
EVIDENCESTRINGSFILENAME = 'evidence_strings.json'
EVIDENCERECORDSFILENAME = 'evidence_records.tsv'
UNMAPPEDTRAITSFILENAME = 'unmappedTraits.tsv'
UNAVAILABLEEFOFILENAME = 'unavailableefo.tsv'
NSVLISTFILE = 'nsvlist.txt'
TMPDIR = '/tmp/'


def clinvarToEvidenceStrings(dirOut, allowedClinicalSignificance=None, ignoreTermsFile=None, adaptTermsFile=None):
    trait2EFO, unavailableEFODict = loadEFOMapping(ignoreTermsFile, adaptTermsFile)

    clinicalSignificance2Activity = {'unknown': 'http://identifiers.org/cttv.activity/unknown',
                                     'untested': 'http://identifiers.org/cttv.activity/unknown',
                                     'non-pathogenic': 'http://identifiers.org/cttv.activity/tolerated_by_target',
                                     'probable-non-pathogenic': 'http://identifiers.org/cttv.activity/predicted_tolerated',
                                     'probable-pathogenic': 'http://identifiers.org/cttv.activity/predicted_damaging',
                                     'pathogenic': 'http://identifiers.org/cttv.activity/damaging_to_target',
                                     'drug-response': 'http://identifiers.org/cttv.activity/unknown',
                                     'histocompatibility': 'http://identifiers.org/cttv.activity/unknown',
                                     'other': 'http://identifiers.org/cttv.activity/unknown',
                                     'benign': 'http://identifiers.org/cttv.activity/tolerated_by_target',
                                     'protective': 'http://identifiers.org/cttv.activity/tolerated_by_target',
                                     'not provided': 'http://identifiers.org/cttv.activity/unknown',
                                     'likely benign': 'http://identifiers.org/cttv.activity/predicted_tolerated',
                                     'confers sensitivity': 'http://identifiers.org/cttv.activity/predicted_damaging',
                                     'uncertain significance': 'http://identifiers.org/cttv.activity/unknown',
                                     'likely pathogenic': 'http://identifiers.org/cttv.activity/predicted_damaging',
                                     'conflicting data from submitters': 'http://identifiers.org/cttv.activity/unknown',
                                     'risk factor': 'http://identifiers.org/cttv.activity/predicted_damaging',
                                     'association': 'http://identifiers.org/cttv.activity/damaging_to_target'
                                     }

    if allowedClinicalSignificance is None:
        allowedClinicalSignificance = ['unknown', 'untested', 'non-pathogenic', 'probable-non-pathogenic',
                                       'probable-pathogenic', 'pathogenic', 'drug-response', 'drug response',
                                       'histocompatibility', 'other', 'benign', 'protective', 'not provided',
                                       'likely benign', 'confers sensitivity', 'uncertain significance',
                                       'likely pathogenic', 'conflicting data from submitters', 'risk factor',
                                       'association']
    nsvList = []
    evidenceStringList = []
    nProcessedClinvarRecords = 0
    nPathogenicNoRs = 0
    nMultipleEvidenceStrings = 0
    nMultipleAlleleOrigin = 0
    nGermlineSomatic = 0
    nRecordsWoRecognizedAlleleOrigin = 0
    noVariantToENSGMapping = 0
    nMoreThanOneEfoTerm = 0
    nSameRefAlt = 0
    nMissedStringsUnmappedTraits = 0
    nNsvs = 0
    nValidRsAndNsv = 0
    nNsvSkippedClinicalSignificance = 0
    nNsvSkippedWrongRefAlt = 0
    unmappedTraits = {}
    nUnrecognizedAlleleOrigin = {}
    ensemblGeneIdUris = set()
    traits = set()
    unrecognizedClinicalSignificances = set()
    evidenceList = []  # To store Hellen Parkinson records of the form [RCV, rs, ClinVar trait, EFO url]
    recordCounter = 0
    nTotalClinvarRecords = 0
    skip = 0
    limit = BATCH_SIZE

    answer = urllib.request.urlopen('http://' + HOST + '/cellbase/webservices/rest/v3/hsapiens/feature/clinical/all?source=clinvar&skip=' + str(skip) + '&limit=' + str(limit))
    reader = codecs.getreader("utf-8")
    currResponse = json.load(reader(answer))['response'][0]
    currResultList = currResponse['result']

    # A progress bar is initialized
    # widgets = ['Loading evidence strings: ', progressbar.Percentage(), ' ', progressbar.Bar(marker=progressbar.RotatingMarker()), ' ', progressbar.ETA()]
    # pbar = progressbar.ProgressBar(widgets=widgets, maxval=currResponse['numTotalResults']).start()
    print(str(currResponse['numTotalResults']) + ' ClinVar records in total.')
    while len(currResultList) > 0:
        nTotalClinvarRecords += len(currResultList)
        for record in currResultList:
            nEvidenceStringsPerRecord = 0
            clinvarRecord = clinvar_record.ClinvarRecord(record['clinvarSet'])
            clinicalSignificance = clinvarRecord.get_clinical_significance().lower()
            nNsvs += (clinvarRecord.get_nsv() is not None)
            if clinicalSignificance in allowedClinicalSignificance:
                if record['reference'] != record['alternate']:
                    nsvList = appendNsv(nsvList, clinvarRecord)
                    rs = clinvarRecord.get_rs()
                    if rs is not None:
                        consequenceType = clinvarRecord.get_main_consequence_types()
                        # Mapping rs->Gene was found at Mick's file and therefore ensemblGeneId will never be None
                        if consequenceType is not None:

                            for ensemblGeneId in consequenceType.get_ensembl_gene_ids():

                                rcvToGeneEvidenceCodes = ['http://identifiers.org/eco/cttv_mapping_pipeline']  # Evidence codes provided by Mick
                                ensemblGeneIdUri = 'http://identifiers.org/ensembl/' + ensemblGeneId
                                traitRefsList = [['http://europepmc.org/abstract/MED/' + str(ref) for ref in refList] for refList in clinvarRecord.get_trait_pubmed_refs()]
                                observedRefsList = ['http://europepmc.org/abstract/MED/' + str(ref) for ref in clinvarRecord.get_observed_pubmed_refs()]
                                measureSetRefsList = ['http://europepmc.org/abstract/MED/' + str(ref) for ref in clinvarRecord.get_measure_set_pubmed_refs()]
                                for traitCounter, traitList in enumerate(clinvarRecord.get_traits()):
                                    clinvarTraitList, EFOList = mapEFO(trait2EFO, traitList)
                                    # Only ClinVar records associated to a trait with mapped EFO term will generate evidence_strings
                                    if len(EFOList) > 0:
                                        clinvaRecordAlleleOrigins = clinvarRecord.get_allele_origins()
                                        nMultipleAlleleOrigin += (len(clinvaRecordAlleleOrigins) > 1)
                                        nGermlineSomatic += (('germline' in clinvaRecordAlleleOrigins) and (
                                        'somatic' in clinvaRecordAlleleOrigins))
                                        nRecordsWoRecognizedAlleleOrigin += (
                                        ('germline' not in clinvaRecordAlleleOrigins) and (
                                        'somatic' not in clinvaRecordAlleleOrigins))
                                        for alleleOriginCounter, alleleOrigin in enumerate(clinvaRecordAlleleOrigins):
                                            if alleleOrigin == 'germline':
                                                evidenceString, nMoreThanOneEfoTerm = getCTTVGeneticsEvidenceString(EFOList,
                                                                                                                    clinicalSignificance,
                                                                                                                    clinicalSignificance2Activity,
                                                                                                                    clinvarRecord,
                                                                                                                    consequenceType,
                                                                                                                    ensemblGeneId,
                                                                                                                    ensemblGeneIdUri,
                                                                                                                    ensemblGeneIdUris,
                                                                                                                    measureSetRefsList,
                                                                                                                    nMoreThanOneEfoTerm,
                                                                                                                    observedRefsList,
                                                                                                                    rcvToGeneEvidenceCodes,
                                                                                                                    record,
                                                                                                                    rs,
                                                                                                                    traitCounter,
                                                                                                                    traitRefsList,
                                                                                                                    traits,
                                                                                                                    unrecognizedClinicalSignificances)
                                                nEvidenceStringsPerRecord = addEvidenceString(clinvarRecord, evidenceString,
                                                                                              evidenceStringList,
                                                                                              nEvidenceStringsPerRecord)
                                                evidenceList.append(
                                                    [clinvarRecord.get_acc(), rs, ','.join(clinvarTraitList),
                                                     ','.join(EFOList)])
                                                nValidRsAndNsv += (clinvarRecord.get_nsv() is not None)
                                            elif alleleOrigin == 'somatic':
                                                evidenceString, nMoreThanOneEfoTerm = getCTTVSomaticEvidenceString(EFOList,
                                                                                                                   clinicalSignificance,
                                                                                                                   clinicalSignificance2Activity,
                                                                                                                   clinvarRecord,
                                                                                                                   ensemblGeneId,
                                                                                                                   ensemblGeneIdUri,
                                                                                                                   ensemblGeneIdUris,
                                                                                                                   measureSetRefsList,
                                                                                                                   nMoreThanOneEfoTerm,
                                                                                                                   observedRefsList,
                                                                                                                   traitCounter,
                                                                                                                   traitRefsList,
                                                                                                                   traits,
                                                                                                                   unrecognizedClinicalSignificances,
                                                                                                                   consequenceType)
                                                nEvidenceStringsPerRecord = addEvidenceString(clinvarRecord, evidenceString,
                                                                                              evidenceStringList,
                                                                                              nEvidenceStringsPerRecord)
                                                evidenceList.append(
                                                    [clinvarRecord.get_acc(), rs, ','.join(clinvarTraitList),
                                                     ','.join(EFOList)])
                                                nValidRsAndNsv += (clinvarRecord.get_nsv() is not None)
                                            elif alleleOrigin not in nUnrecognizedAlleleOrigin:
                                                nUnrecognizedAlleleOrigin[alleleOrigin] = 1
                                            else:
                                                nUnrecognizedAlleleOrigin[alleleOrigin] += 1
                                    else:
                                        nMissedStringsUnmappedTraits += 1
                                        if traitList[0] in unmappedTraits:
                                            unmappedTraits[traitList[0]] += 1
                                        else:
                                            unmappedTraits[traitList[0]] = 1

                                if nEvidenceStringsPerRecord > 0:
                                    nProcessedClinvarRecords += 1
                                    if nEvidenceStringsPerRecord > 1:
                                        nMultipleEvidenceStrings += 1
                        else:
                            noVariantToENSGMapping += 1
                    else:
                        nPathogenicNoRs += 1

                else:
                    nSameRefAlt += 1
                    if clinvarRecord.get_nsv() is not None:
                        nNsvSkippedWrongRefAlt += 1
            else:
                if clinvarRecord.get_nsv() is not None:
                    nNsvSkippedClinicalSignificance += 1

            # pbar.update(recordCounter)
            recordCounter += 1
        skip += BATCH_SIZE

        answer = urllib.request.urlopen('http://' + HOST + '/cellbase/webservices/rest/v3/hsapiens/feature/clinical/all?source=clinvar&skip=' + str(skip) + '&limit=' + str(limit))
        reader = codecs.getreader("utf-8")
        currResponse = json.load(reader(answer))['response'][0]
        currResultList = currResponse['result']
    # pbar.finish()

    writeStringListToFile(nsvList, dirOut + '/' + NSVLISTFILE)

    fdw = open(dirOut + '/' + UNMAPPEDTRAITSFILENAME, 'w')  # Contains traits without a mapping in Gary's xls
    fdw.write('Trait\tCount\n')
    for traitList in unmappedTraits:
        fdw.write(str(traitList.encode('utf8')) + '\t' + str(unmappedTraits[traitList]) + '\n')
    fdw.close()

    fdw = open(dirOut + '/' + UNAVAILABLEEFOFILENAME,
               'w')  # Contains urls provided by Gary which are not yet included within EFO
    fdw.write('Trait\tCount\n')
    for url in unavailableEFODict:
        fdw.write(url.encode('utf8') + '\t' + str(unavailableEFODict[url]) + '\n')
    fdw.close()

    fdw = open(dirOut + '/' + EVIDENCESTRINGSFILENAME, 'w')
    for evidenceString in evidenceStringList:
        fdw.write(json.dumps(evidenceString) + '\n')
    fdw.close()

    fdw = open(dirOut + '/' + EVIDENCERECORDSFILENAME, 'w')
    for evidenceRecord in evidenceList:
        fdw.write('\t'.join(evidenceRecord) + '\n')
    fdw.close()

    print(str(nTotalClinvarRecords) + ' ClinVar records in total')
    print(str(len(evidenceStringList)) + ' evidence string jsons generated')
    print(str(nProcessedClinvarRecords) + ' ClinVar records generated at least one evidence string')
    print(str(len(
        unrecognizedClinicalSignificances)) + " Clinical significance string(s) not found among those described in ClinVar documentation:")
    print(str(unrecognizedClinicalSignificances))
    print(str(
        nSameRefAlt) + ' ClinVar records with allowed clinical significance did present the same reference and alternate and were skipped')
    print('Activities of those ClinVar records with unrecognized clinical significances were set to "unknown".')
    print(str(len(ensemblGeneIdUris)) + ' distinct ensembl gene ids appear in generated evidence string json objects')
    print(str(len(traits)) + ' distinct trait names found to include in generated evidence string json objects')
    print(str(nPathogenicNoRs) + ' ClinVar records with allowed clinical significance DO NOT have an rs id')
    print(str(nMultipleEvidenceStrings) + ' ClinVar records generated more than one evidence_string')
    print(str(nGermlineSomatic) + ' ClinVar records with germline and somatic origins')
    print(str(nMultipleAlleleOrigin) + ' ClinVar records with more than one allele origin')
    print('Number valid ClinVar records with unprocessed allele origins:')
    for alleleOrigin in nUnrecognizedAlleleOrigin:
        print(' ' + alleleOrigin + ': ' + str(nUnrecognizedAlleleOrigin[alleleOrigin]))
    print(str(
        noVariantToENSGMapping) + ' ClinVar records with allowed clinical significance and valid rs id were skipped due to a lack of Variant->ENSG mapping.')
    print(str(
        nMissedStringsUnmappedTraits) + ' ClinVar records with allowed clinical significance, valid rs id and Variant->ENSG mapping were skipped due to a lack of EFO mapping (see ' + UNMAPPEDTRAITSFILENAME + ').')
    print(str(
        nRecordsWoRecognizedAlleleOrigin) + ' ClinVar records with allowed clinical significance, valid rs id, valid Variant->ENSG mapping and valid EFO mapping were skipped due to a lack of a valid alleleOrigin.')
    print(str(nMoreThanOneEfoTerm) + ' evidence strings with more than one trait mapped to EFO terms')
    print(str(len(unavailableEFODict)) + ' evidence strings were generated with traits without EFO correspondence')
    print(str(nValidRsAndNsv) + ' evidence strings were generated from ClinVar records with rs and nsv ids')
    print(str(nNsvs) + ' total nsvs found')
    print(str(
        nNsvSkippedClinicalSignificance) + ' ClinVar nsvs were skipped because of a different clinical significance')
    print(str(nNsvSkippedWrongRefAlt) + ' ClinVar nsvs were skipped because of same ref and alt')


def getCTTVGeneticsEvidenceString(EFOList, clinicalSignificance, clinicalSignificance2Activity, clinvarRecord,
                                  consequenceType, ensemblGeneId, ensemblGeneIdUri, ensemblGeneIdUris,
                                  measureSetRefsList, nMoreThanOneEfoTerm, observedRefsList, rcvToGeneEvidenceCodes,
                                  record, rs, traitCounter, traitRefsList, traits,
                                  unrecognizedClinicalSignificances):
    evidenceString = evidence_strings.CTTVGeneticsEvidenceString()
    evidenceString.addUniqueAssociationField('gene', ensemblGeneId)
    evidenceString.addUniqueAssociationField('clinvarAccession', clinvarRecord.get_acc())
    evidenceString.addUniqueAssociationField('alleleOrigin', 'germline')
    try:
        evidenceString.setTarget(ensemblGeneIdUri, clinicalSignificance2Activity[clinicalSignificance])
    except KeyError:
        unrecognizedClinicalSignificances.add(clinicalSignificance)
        evidenceString.setTarget(ensemblGeneIdUri, 'http://identifiers.org/cttv.activity/unknown')
    evidenceString.setVariant('http://identifiers.org/dbsnp/' + rs, getCttvVariantType(record))
    evidenceString.setDate(clinvarRecord.get_date())
    evidenceString.setDbxrefUrl('http://identifiers.org/clinvar.record/' + clinvarRecord.get_acc())
    evidenceString.setUrl('http://www.ncbi.nlm.nih.gov/clinvar/' + clinvarRecord.get_acc())
    evidenceString.setAssociation(
        clinicalSignificance != 'non-pathogenic' and clinicalSignificance != 'probable-non-pathogenic'
        and clinicalSignificance != 'likely benign' and clinicalSignificance != 'benign')
    evidenceString.setGene2VariantEvidenceCodes(rcvToGeneEvidenceCodes)
    mostSevereSoTerm = consequenceType.getMostSevereSo()
    if mostSevereSoTerm.getAccession() is None:
        evidenceString.setGene2VariantFunctionalConsequence(
            'http://targetvalidation.org/sequence/' + mostSevereSoTerm.getName())
    else:
        evidenceString.setGene2VariantFunctionalConsequence(
            'http://purl.obolibrary.org/obo/' + mostSevereSoTerm.getAccession().replace(':', '_'))

    referenceList = list(set(traitRefsList[traitCounter] + observedRefsList + measureSetRefsList))
    if len(referenceList) > 0:
        evidenceString.setVariant2DiseaseLiterature(referenceList)
        # Arbitrarily select only one reference among all
        evidenceString.setUniqueReference(referenceList[0])
        evidenceString.set_top_level_literature(referenceList)
    EFOList.sort()
    # Just (arbitrarily) adding one of the potentially multiple EFO terms because of schema constraints
    evidenceString.setDisease(EFOList[0])
    evidenceString.addUniqueAssociationField('phenotype', EFOList[0])
    nMoreThanOneEfoTerm += (len(EFOList) > 1)
    traits.update(set(EFOList))
    ensemblGeneIdUris.add(ensemblGeneIdUri)
    return evidenceString, nMoreThanOneEfoTerm


def getCTTVSomaticEvidenceString(EFOList, clinicalSignificance, clinicalSignificance2Activity, clinvarRecord,
                                 ensemblGeneId, ensemblGeneIdUri, ensemblGeneIdUris, measureSetRefsList,
                                 nMoreThanOneEfoTerm, observedRefsList, traitCounter, traitRefsList, traits,
                                 unrecognizedClinicalSignificances, consequenceType):
    evidenceString = evidence_strings.CTTVSomaticEvidenceString()
    evidenceString.addUniqueAssociationField('gene', ensemblGeneId)
    evidenceString.addUniqueAssociationField('clinvarAccession', clinvarRecord.get_acc())
    evidenceString.addUniqueAssociationField('alleleOrigin', 'somatic')
    try:
        evidenceString.setTarget(ensemblGeneIdUri, clinicalSignificance2Activity[clinicalSignificance])
    except KeyError:
        unrecognizedClinicalSignificances.add(clinicalSignificance)
        evidenceString.setTarget(ensemblGeneIdUri, 'http://identifiers.org/cttv.activity/unknown')

    evidenceString.setDate(clinvarRecord.get_date())
    evidenceString.setDbxrefUrl('http://identifiers.org/clinvar.record/' + clinvarRecord.get_acc())
    evidenceString.setUrl('http://www.ncbi.nlm.nih.gov/clinvar/' + clinvarRecord.get_acc())
    evidenceString.setAssociation(
        clinicalSignificance != 'non-pathogenic' and clinicalSignificance != 'probable-non-pathogenic'
        and clinicalSignificance != 'likely benign' and clinicalSignificance != 'benign')

    evidenceString.set_known_mutations(consequenceType)

    referenceList = list(set(traitRefsList[traitCounter] + observedRefsList + measureSetRefsList))
    if len(referenceList) > 0:
        evidenceString.set_evidence_literature(referenceList)
        evidenceString.set_top_level_literature(referenceList)

    EFOList.sort()
    # Just (arbitrarily) adding one of the potentially multiple EFO terms because of schema constraints
    evidenceString.setDisease(EFOList[0])
    evidenceString.addUniqueAssociationField('phenotype', EFOList[0])
    nMoreThanOneEfoTerm += (len(EFOList) > 1)
    traits.update(set(EFOList))
    ensemblGeneIdUris.add(ensemblGeneIdUri)
    return evidenceString, nMoreThanOneEfoTerm


def addEvidenceString(clinvarRecord, evidenceString, evidenceStringList, nEvidenceStringsPerRecord):
    try:
        evidenceString.validate()
        evidenceStringList.append(evidenceString)
        nEvidenceStringsPerRecord += 1
    except jsonschema.exceptions.ValidationError as err:
        print('Error: evidence_string does not validate against schema.')
        print('ClinVar accession: ' + clinvarRecord.get_acc())
        print(err)
        print(json.dumps(evidenceString))
        sys.exit(1)
    except EFOTerm.EFOTerm.IsObsoleteException as err:
        print('Error: obsolete EFO term.')
        print('Term: ' + evidenceString.get_disease().getId())
        print(err)
        print(json.dumps(evidenceString))
        sys.exit(1)

    return nEvidenceStringsPerRecord


def writeStringListToFile(stringList, filename):
    fd = open(filename, 'w')
    fd.write('\n'.join(stringList))
    fd.close()


def appendNsv(nsvList, clinvarRecord):
    nsv = clinvarRecord.get_nsv()
    if nsv is not None:
        nsvList.append(nsv)
    return nsvList


def getCttvVariantType(record):
    if len(record['reference']) < 2 and len(record['alternate']) < 2:
        cttvVariantType = 'snp single'
    elif len(record['reference']) > 50 or len(record['alternate']) > 50:
        cttvVariantType = 'structural variant'
    else:
        cttvVariantType = 'snp single'  # Sam asked for this in his email 21/05/2015
        # cttvVariantType = 'snp multiple'

    return cttvVariantType


def mapEFO(trait2EFO, traitList):
    EFOList = []
    traitListToReturn = []
    traitString = traitList[0].lower()
    if traitString in trait2EFO:
        for EFOtrait in trait2EFO[traitString]:
            if EFOtrait not in EFOList:  # First element in traitList mus always be the "Preferred" trait name
                traitListToReturn.append(traitList[0])
                EFOList.append(EFOtrait)
    else:
        for trait in traitList[1:]:
            traitString = trait.lower()
            if traitString in trait2EFO:
                for EFOtrait in trait2EFO[traitString]:
                    if EFOtrait not in EFOList:  # First element in traitList mus always be the "Preferred" trait name
                        traitListToReturn.append(trait)
                        EFOList.append(EFOtrait)

    return traitListToReturn, EFOList


def loadEFOMapping(ignoreTermsFile=None, adaptTermsFile=None):
    ignoreTerms = getTermsFromFile(ignoreTermsFile)
    adaptTerms = getTermsFromFile(adaptTermsFile)

    print('Loading phenotypes to EFO mapping...')
    EFOMappingReadBook = xlrd.open_workbook(EFOMAPPINGFILE, formatting_info=True)
    EFOMappingReadSheet = EFOMappingReadBook.sheet_by_index(0)
    trait2EFO = {}
    unavailableEFO = {}
    nEFOmappings = 0
    for i in range(1, EFOMappingReadSheet.nrows):
        if EFOMappingReadSheet.cell_value(rowx=i, colx=1) != '':
            validEFO, urlsToAdapt = getUrls(EFOMappingReadSheet.cell_value(rowx=i, colx=1).split(', '), ignoreTerms,
                                            adaptTerms)
            clinvarTrait = EFOMappingReadSheet.cell_value(rowx=i, colx=0).lower()
            if len(validEFO) > 0:
                trait2EFO[clinvarTrait] = validEFO
                nEFOmappings += 1
            elif len(urlsToAdapt) > 0:
                trait2EFO[clinvarTrait] = []
                for url in urlsToAdapt:
                    if url not in unavailableEFO:
                        unavailableEFO[url] = 1
                    else:
                        unavailableEFO[url] += 1
                    trait2EFO[clinvarTrait].append(getUnmappedUrl(url))

    print(str(nEFOmappings) + ' EFO mappings loaded')
    print(str(len(unavailableEFO)) + ' urls without an actual valid EFO mapping')

    return trait2EFO, unavailableEFO


def getUnmappedUrl(url):
    parts = url.split('/')
    if parts[-1].startswith("Orphanet_"):
        newUrl = "http://purl.bioontology.org/ORDO/" + parts[-1]
    elif parts[-1].startswith("HP_"):
        newUrl = "http://purl.bioontology.org/obo/" + parts[-1]
    else:
        print("Error. Unhandled url type: " + url)
        sys.exit(1)

    return newUrl


def getUrls(urlList, ignoreTerms, adaptTerms):
    validEFO = []
    urlsToAdapt = []
    for term in urlList:
        if term not in ignoreTerms:
            if term in adaptTerms:
                urlsToAdapt.append(term)
            else:
                validEFO.append(term)

    return validEFO, urlsToAdapt


def getTermsFromFile(termsFile):
    if termsFile is not None:
        print('Loading list of terms...')
        fd = open(termsFile, 'r')
        termList = [line.rstrip() for line in fd]
        fd.close()
        print(str(len(termsFile)) + ' terms found at ' + termsFile)
    else:
        termList = []

    return termList


class ArgParser(object):
    """
    For parsing command line arguments
    """
    def __init__(self, argv):
        usage = """
        ************************************************************************************************************************************************************
        Task: generate CTTV evidence strings from ClinVar mongo
        ************************************************************************************************************************************************************



        usage: %prog --clinSig <clinicalSignificanceList> --out <fileout>"""
        parser = argparse.ArgumentParser(usage)

        parser.add_argument("--clinSig", dest="clinSig", help="""Optional. String containing a comma-sparated list with the clinical significances that will be allowed to generate evidence-strings. By default all clinical significances will be considered. Possible tags: 'unknown','untested','non-pathogenic','probable-non-pathogenic','probable-pathogenic','pathogenic','drug-response','drug response','histocompatibility','other','benign','protective','not provided','likely benign','confers sensitivity','uncertain significance','likely pathogenic','conflicting data from submitters','risk factor','association' """, default="pathogenic,likely pathogenic")
        parser.add_argument("--ignore", dest="ignoreTermsFile", help="""Optional. String containing full path to a txt file containing a list of term urls which will be ignored during batch processing """, default=None)
        parser.add_argument("--adapt", dest="adaptTermsFile", help="""Optional. String containing full path to a txt file containing a list of invalid EFO urls which will be adapted to a general valid url during batch processing """, default=None)
        parser.add_argument("--out", dest="out", help="""String containing the name of the file were results will be stored.""", required=True)

        args = parser.parse_args(args=argv[1:])

        self.clinSig = args.clinSig
        self.ignoreTermsFile = args.ignoreTermsFile
        self.adaptTermsFile = args.adaptTermsFile
        self.out = args.out


def main():
    parser = ArgParser(sys.argv)

    utilities.check_for_local_schema()

    # call core function
    if parser.clinSig is None:
        clinvarToEvidenceStrings(parser.out, ignoreTermsFile=parser.ignoreTermsFile,
                                 adaptTermsFile=parser.adaptTermsFile)
    else:
        clinvarToEvidenceStrings(parser.out, allowedClinicalSignificance=parser.clinSig.split(','),
                                 ignoreTermsFile=parser.ignoreTermsFile, adaptTermsFile=parser.adaptTermsFile)

    print('>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> Finished <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<')


if __name__ == '__main__':
    main()
