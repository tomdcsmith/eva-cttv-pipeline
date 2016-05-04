local_schema = "resources/schema_local"

#########################################
# clinvar_to_evidence_strings.py settings
#########################################
# cellbase
BATCH_SIZE = 200
HOST = 'www.ebi.ac.uk'

# output settings
EVIDENCESTRINGSFILENAME = 'evidence_strings.json'
EVIDENCERECORDSFILENAME = 'evidence_records.tsv'
UNMAPPEDTRAITSFILENAME = 'unmappedTraits.tsv'
UNAVAILABLEEFOFILENAME = 'unavailableefo.tsv'
NSVLISTFILE = 'nsvlist.txt'
TMPDIR = '/tmp/'


sparqlep = ""
sparql_user = ""
sparql_password = ""
