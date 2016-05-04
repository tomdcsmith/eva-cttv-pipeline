local_schema = "resources/schema_local"


sparqlep = ""
sparql_user = ""
sparql_password = ""


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

######


##############################
# evidence_strings.py settings
##############################
GEN_EV_STRING_JSON = "resources/CTTVGeneticsEvidenceString.json"
SOM_EV_STRING_JSON = "resources/CTTVSomaticEvidenceString.json"
GEN_SCHEMA_FILE = local_schema + "/src/genetics.json"
SOM_SCHEMA_FILE = local_schema + "/src/literature_curated.json"

#############
