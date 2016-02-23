from eva_cttv_pipeline import EFOTerm

__author__ = 'Javier Lopez: javild@gmail.com'


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
