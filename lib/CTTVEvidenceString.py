__author__ = 'Javier Lopez: javild@gmail.com'


import EFOTerm


class CTTVEvidenceString(dict):
    def __init__(self, aDictionary):
        dict.__init__(self, aDictionary)

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

    def getDisease(self):
        return EFOTerm.EFOTerm(self['disease']['id'][0])

    def setEvidenceCodes(self, evidenceCodeList):
        self['evidence']['evidence_codes'] = evidenceCodeList

    def getEvidenceCodes(self):
        return self['evidence']['evidence_codes']
