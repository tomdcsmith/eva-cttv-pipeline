__author__ = 'Javier Lopez: javild@gmail.com'

import os


def getAvailableTerms(sparqlep, user, password):
    query = """SELECT DISTINCT ?uri sample(?label) AS ?label
    FROM <http://www.targetvalidation.org/>
    FROM <http://www.ebi.ac.uk/efo/>
    WHERE {
        ?uri rdfs:subClassOf* <http://www.targetvalidation.org/cttv_root> .
        ?uri rdfs:label ?label .
    }"""

    return executeQuery(sparqlep, user, password, query)


def getObsoleteTerms(sparqlep, user, password):
    # What EFO classes are now obsolete and what is the reason for obsoletion?
    query = """SELECT * FROM <http://www.ebi.ac.uk/efo/>
     WHERE
     {
              ?uri rdfs:subClassOf <http://www.geneontology.org/formats/oboInOwl#ObsoleteClass> .
              ?uri <http://www.ebi.ac.uk/efo/reason_for_obsolescence> ?reason
     }"""

    return executeQuery(sparqlep, user, password, query)


def executeQuery(sparqlep, user, password, query):
    cmd = """curl --post301 -L --data-urlencode "query=""" + query + """" -u """ + user + ":" + password + """ -H "Accept: text/tab-separated-values" """ + sparqlep
    fd = os.popen(cmd)
    print('Loading obsolete terms...')
    terms = {}
    fd.readline()
    for line in fd.readlines():
        parts = line.rstrip("\n").split('\t')
        if len(parts) > 1:
            terms[parts[0].rstrip('"').lstrip('"')] = parts[1].lstrip('"').rstrip('"')
    print('Done.')

    return terms


class EFOTerm():
    print('Querying EFO  web services for valid terms...')
    # CONNECTION PARAMETERS
    sparqlep = ""
    user = ""
    password = ""

    obsoleteTerms = getObsoleteTerms(sparqlep, user, password)
    cttvAvailableTerms = getAvailableTerms(sparqlep, user, password)

    def __init__(self, _id=None):
        self.id = _id

    def getId(self):
        return self.id

    def isObsolete(self):
        if self.id is not None and self.id in EFOTerm.obsoleteTerms:
            raise EFOTerm.IsObsoleteException(
                "Term " + self.id + " is obsolete. Cause/Description/Action: " + EFOTerm.obsoleteTerms[self.id])
        else:
            return False

    def isCttvAvailable(self):
        if self.id is not None and self.id in EFOTerm.cttvAvailableTerms:
            return True
        else:
            raise EFOTerm.NotCttvAvailableException("Term " + self.id + " is not currently available at EFO.")

    class IsObsoleteException(Exception):
        def __init__(self, value):
            self.value = value

        def __str__(self):
            return repr(self.value)

    class NotCttvAvailableException(Exception):
        def __init__(self, value):
            self.value = value

        def __str__(self):
            return repr(self.value)
