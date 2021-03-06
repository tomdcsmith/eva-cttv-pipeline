import os

from eva_cttv_pipeline.evidence_string_generation import config


def get_available_terms(sparqlep, user, password):
    query = """SELECT DISTINCT ?uri sample(?label) AS ?label
    FROM <http://www.targetvalidation.org/>
    FROM <http://www.ebi.ac.uk/efo/>
    WHERE {
        ?uri rdfs:subClassOf* <http://www.targetvalidation.org/cttv_root> .
        ?uri rdfs:label ?label .
    }"""

    return execute_query(sparqlep, user, password, query)


def get_obsolete_terms(sparqlep, user, password):
    # What EFO classes are now obsolete and what is the reason for obsoletion?
    query = """SELECT * FROM <http://www.ebi.ac.uk/efo/>
     WHERE
     {
              ?uri rdfs:subClassOf <http://www.geneontology.org/formats/oboInOwl#ObsoleteClass> .
              ?uri <http://www.ebi.ac.uk/efo/reason_for_obsolescence> ?reason
     }"""

    return execute_query(sparqlep, user, password, query)


def execute_query(sparqlep, user, password, query):
    cmd = "curl --post301 -L --data-urlencode \"query=" + query + " -u " + user + ":" + \
          password + " -H \"Accept: text/tab-separated-values\" " + sparqlep
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


class EFOTerm:

    """
    Experimental factor ontology term.
    Includes methods to check for obsoleteness and if it is available for Open Targets.
    Used in pipeline by mapping from a trait to the corresponding EFO term.
    """

    print('Querying EFO  web services for valid terms...')
    # CONNECTION PARAMETERS
    sparqlep = config.SPARQLEP
    user = config.SPARQL_USER
    password = config.SPARQL_PASSWORD

    obsolete_terms = get_obsolete_terms(sparqlep, user, password)
    cttv_available_terms = get_available_terms(sparqlep, user, password)

    def __init__(self, efoid=None):
        self.efoid = efoid

    def __eq__(self, other):
        return isinstance(other, self.__class__) and self.efoid == other.efoid

    def __ne__(self, other):
        return not self.__eq__(other)

    def __repr__(self):
        return "%s(%r)" % (self.__class__, self.__dict__)

    def is_obsolete(self):
        if self.efoid is not None and self.efoid in EFOTerm.obsolete_terms:
            raise EFOTerm.IsObsoleteException(
                "Term " + self.efoid + " is obsolete. Cause/Description/Action: " +
                EFOTerm.obsolete_terms[self.efoid])
        else:
            return False

    def is_cttv_available(self):
        if self.efoid is not None and self.efoid in EFOTerm.cttv_available_terms:
            return True
        else:
            raise EFOTerm.TermUnavailableException("Term " + self.efoid +
                                                    " is not currently available at EFO.")

    class IsObsoleteException(Exception):
        def __init__(self, value):
            Exception.__init__(self, value)

    class TermUnavailableException(Exception):
        def __init__(self, value):
            Exception.__init__(self, value)
