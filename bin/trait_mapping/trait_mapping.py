import argparse
from collections import Counter
import csv
from enum import Enum
from functools import lru_cache, total_ordering
import gzip
import json
import re
import sys
import urllib

import requests
import progressbar


##
# Classes
##


class OntologyUri:
    db_to_uri_dict = {
        "orphanet": "http://www.orpha.net/ORDO/Orphanet_{}",
        "omim": "http://identifiers.org/omim/{}",
        "efo": "http://www.ebi.ac.uk/efo/EFO_{}",
        "mesh": "http://identifiers.org/mesh/{}",
        "medgen": "http://identifiers.org/medgen/{}",
        "human phenotype ontology": "http://purl.obolibrary.org/obo/HP_{}",
        "hp": "http://purl.obolibrary.org/obo/HP_{}"
    }

    def __init__(self, id_, db):
        self.id_ = id_
        self.db = db
        if self.db.lower() == "human phenotype ontology":
            self.uri = self.db_to_uri_dict[self.db.lower()].format(self.id_[3:])
        else:
            self.uri = self.db_to_uri_dict[self.db.lower()].format(self.id_)

    def __str__(self):
        return self.uri


@total_ordering
class OxOMapping:
    """
    Individual mapping for an ontology ID mapped to one other ontology ID. An OxO result can consist
    of multiple mappings.
    """
    def __init__(self, label, curie, distance, query_id):
        self.label = label
        self.db, self.id_ = curie.split(":")
        self.uri = OntologyUri(self.id_, self.db)
        self.distance = distance
        self.query_id = query_id
        self.in_efo = False
        self.is_current = False
        self.ontology_label = ""

    def __eq__(self, other):
        if not isinstance(other, type(self)):
            return False
        return (self.label == other.label, self.db == other.db, self.id_ == other.id_,
                self.distance == other.distance, self.in_efo == other.in_efo,
                self.is_current == other.is_current, self.ontology_label == other.ontology_label)

    def __lt__(self, other):
        return ((self.distance, self.in_efo, self.is_current) <
                (other.distance, other.in_efo, other.is_current))


class OxOResult:
    """
    A single result from querying OxO for one ID. A result can contain multiple mappings. A response
    from OxO can contain multiple results- one per queried ID.
    """
    def __init__(self, query_id, label, curie):
        self.query_id = query_id
        self.label = label
        self.db, self.id_ = curie.split(":")
        self.uri = OntologyUri(self.id_, self.db)
        self.oxo_mapping_list = []


@total_ordering
class ZoomaConfidence(Enum):
    """Enum to represent the confidence of a mapping in Zooma."""
    LOW = 1
    MEDIUM = 2
    GOOD = 3
    HIGH = 4

    def __eq__(self, other):
        if not isinstance(other, type(self)):
            return False
        return self.value == other.value

    def __lt__(self, other):
        return self.value < other.value

    def __str__(self):
        return self.name


@total_ordering
class ZoomaEntry:
    """Representation of one ontology term in a mapping in Zooma."""
    def __init__(self, uri, confidence, source):
        self.uri = uri
        self.confidence = ZoomaConfidence[confidence.upper()]
        self.source = source
        self.ontology_label = ""
        self.in_efo = False
        self.is_current = False

    def __eq__(self, other):
        if not isinstance(other, type(self)):
            return False
        if (self.uri != other.uri or self.confidence != other.confidence
                or self.ontology_label != other.ontology_label or self.in_efo != other.in_efo
                or self.is_current != other.is_current):
            return False
        return True

    def __lt__(self, other):
        return ((self.confidence, self.in_efo, self.is_current) <
                (other.confidence, other.in_efo, other.is_current))


class ZoomaMapping:
    """
    A mapping in Zooma from one term, which can contain multiple ontology IDs mapped to. One
    term can be mapped to multiple mappings.
    """
    def __init__(self, uri_list, zooma_label, confidence, source):
        self.zooma_label = zooma_label
        self.confidence = confidence
        self.source = source
        self.zooma_entry_list = []
        for uri in uri_list:
            self.zooma_entry_list.append(ZoomaEntry(uri, confidence, source))


class OntologyEntry:
    """
    A representation of an ontology term mapped to using this pipeline. Includes a uri and a label
    which are the two details needed in a finished mapping, in addition to the ClinVar trait name.
    """
    def __init__(self, uri, label):
        self.uri = uri
        self.label = label

    def __eq__(self, other):
        if not isinstance(other, type(self)):
            return False
        if self.uri != other.uri or self.label != other.label:
            return False
        return True

    def __hash__(self):
        return hash((self.uri, self.label))


class Trait:
    """
    Object to hold data for one trait name. Including the number of ClinVar record's traits it
    appears in, any Zooma and OxO mappings, and any mappings which are ready to be output.
    """
    def __init__(self, name, frequency):
        self.name = name
        self.frequency = frequency
        self.zooma_mapping_list = []
        self.oxo_xref_list = []
        self.finished_mapping_set = set()

    @property
    def is_finished(self):
        """
        Return boolean which confirms whether the trait has finished mappings ready to be output
        """
        return len(self.finished_mapping_set) > 0

    def process_zooma_mappings(self):
        """
        Check whether any Zooma mappings can be output as a finished ontology mapping.
        Put any finished mappings in finished_mapping_set
        """
        for mapping in self.zooma_mapping_list:
            if mapping.confidence.lower() != "high":
                continue

            for entry in mapping.zooma_entry_list:
                if entry.in_efo and entry.is_current:
                    ontology_entry = OntologyEntry(entry.uri, entry.ontology_label)
                    self.finished_mapping_set.add(ontology_entry)

    def process_oxo_mappings(self):
        """
        Check whether any OxO mappings can be output as a finished ontology mapping.
        Put any finished mappings in finished_mapping_set
        """
        for result in self.oxo_xref_list:
            for mapping in result.oxo_mapping_list:
                if mapping.in_efo and mapping.is_current and mapping.distance == 1:
                    uri = str(mapping.uri)
                    ontology_label = mapping.ontology_label
                    ontology_entry = OntologyEntry(uri, ontology_label)
                    self.finished_mapping_set.add(ontology_entry)


##
# Main method
##


def main():
    parser = ArgParser(sys.argv)

    trait_names_list = parse_trait_names(parser.input_filepath)
    trait_names_counter = Counter(trait_names_list)

    with open(parser.output_mappings_filepath, "w", newline='') as mapping_file, open(parser.output_curation_filepath, "wt") as curation_file:
        mapping_writer = csv.writer(mapping_file, delimiter="\t")
        mapping_writer.writerow(["#clinvar_trait_name", "uri", "label"])
        curation_writer = csv.writer(curation_file, delimiter="\t")

        bar = progressbar.ProgressBar(max_value=len(trait_names_counter),
                                      widgets=[progressbar.AdaptiveETA(samples=1000)])

        for trait_name, freq in bar(trait_names_counter.items()):
            trait = Trait(trait_name, freq)
            trait = process_trait(trait, parser.filters, parser.zooma_host, parser.oxo_target_list,
                                  parser.oxo_distance)
            output_trait(trait, mapping_writer, curation_writer)


def output_trait_mapping(trait, mapping_writer):
    """
    Write any finished ontology mappings for a trait to a csv file writer.

    :param trait: A trait with finished ontology mappings in finished_mapping_set
    :param mapping_writer: A csv.writer to write the finished mappings
    """
    for ontology_entry in trait.finished_mapping_set:
        mapping_writer.writerow([trait.name, ontology_entry.uri, ontology_entry.label])


def output_for_curation(trait, curation_writer):
    """
    Write any non-finished Zooma or OxO mappings of a trait to a file for manual curation.
    Also outputs traits without any ontology mappings.

    :param trait: A trait with no finished ontology mappings in finished_mapping_set
    :param curation_writer: A csv.writer to write non-finished ontology mappings for manual curation
    """
    output_row = []
    output_row.extend([trait.name, trait.frequency])

    zooma_entry_list = [zooma_entry for zooma_mapping in trait.zooma_mapping_list for zooma_entry in zooma_mapping.zooma_entry_list]
    zooma_entry_list.sort(reverse=True)

    for zooma_entry in zooma_entry_list:
        cell = [zooma_entry.uri, zooma_entry.ontology_label, str(zooma_entry.in_efo),
                str(zooma_entry.is_current), str(zooma_entry.confidence), zooma_entry.source]
        output_row.append("|".join(cell))

    oxo_mapping_list = [oxo_mapping for oxo_result in trait.oxo_xref_list for oxo_mapping in oxo_result.oxo_mapping_list]
    oxo_mapping_list.sort(reverse=True)

    for oxo_mapping in oxo_mapping_list:
        cell = [str(oxo_mapping.uri), oxo_mapping.ontology_label, str(oxo_mapping.in_efo),
                str(oxo_mapping.is_current), str(oxo_mapping.distance), oxo_mapping.query_id]
        output_row.append("|".join(cell))

    curation_writer.writerow(output_row)


def output_trait(trait, mapping_writer, curation_writer):
    if trait.is_finished:
        output_trait_mapping(trait, mapping_writer)
    else:
        output_for_curation(trait, curation_writer)


def get_uris_for_oxo(zooma_mapping_list):
    uri_set = set()
    for mapping in zooma_mapping_list:
        # Only use high confidence Zooma mappings for querying OxO
        if mapping.confidence.lower() != "high":
            continue
        uri_set.update([entry.uri for entry in mapping.zooma_entry_list])
    return uri_set


def process_trait(trait, filters, zooma_host, oxo_target_list, oxo_distance):
    zooma_mappings = get_ontology_mappings(trait.name, filters, zooma_host)
    trait.zooma_mapping_list = zooma_mappings
    trait.process_zooma_mappings()
    if (trait.is_finished
            or len(trait.zooma_mapping_list) == 0
            or any([entry.is_current for mapping in trait.zooma_mapping_list for entry in mapping.zooma_entry_list])):
        return trait
    uris_for_oxo_set = get_uris_for_oxo(trait.zooma_mapping_list)
    if len(uris_for_oxo_set) == 0:
        return trait
    oxo_input_id_list = uris_to_oxo_format(uris_for_oxo_set)
    oxo_result_list = get_oxo_results(oxo_input_id_list, oxo_target_list, oxo_distance)
    trait.oxo_xref_list = oxo_result_list
    trait.process_oxo_mappings()

    return trait


##
# Parsing trait names
##


def clinvar_jsons(filepath):
    with gzip.open(filepath, "rt") as f:
        for line in f:
            line = line.rstrip()
            yield json.loads(line)


def get_trait_names(clinvar_json):
    # This if-else block is due to the change in the format of the CellBase JSON that holds the
    # ClinVar data. Originally "clinvarSet" was the top level, but this level was removed and
    # referenceClinVarAssertion is now the top level.
    if "clinvarSet" in clinvar_json:
        trait_set = clinvar_json["clinvarSet"]["referenceClinVarAssertion"]["traitSet"]
    else:
        trait_set = clinvar_json["referenceClinVarAssertion"]["traitSet"]
    trait_list = []
    for trait in trait_set['trait']:
        trait_list.append([])
        for name in trait['name']:
            # First trait name in the list will always be the "Preferred" one
            if name['elementValue']['type'] == 'Preferred':
                trait_list[-1] = [name['elementValue']['value']] + trait_list[-1]
            elif name['elementValue']['type'] in ["EFO URL", "EFO id", "EFO name"]:
                continue  # if the trait name not originally from clinvar
            else:
                trait_list[-1].append(name['elementValue']['value'])

    trait_names_to_return = []
    for trait in trait_list:
        if len(trait) == 0:
            continue
        trait_names_to_return.append(trait[0].lower())

    return trait_names_to_return


def parse_trait_names(filepath):
    trait_name_list = []
    for clinvar_json in clinvar_jsons(filepath):
        new_trait_names = get_trait_names(clinvar_json)
        trait_name_list.extend(new_trait_names)
    return trait_name_list


##
# Zooma functions
##


def request_retry_helper(function, retry_count, url):
    for retry_num in range(retry_count):
        return_value = function(url)
        if return_value is not None:
            return return_value
        print("attempt {}: failed running function {} with url {}".format(retry_num, function, url))
    print("error on last attempt, skipping")
    return None


def zooma_query_helper(url):
    try:
        json_response_1 = requests.get(url).json()
        return json_response_1
    except json.decoder.JSONDecodeError as e:
        return None


def ols_query_helper(url):
    try:
        json_response = requests.get(url).json()
        for term in json_response["_embedded"]["terms"]:
            if term["is_defining_ontology"]:
                return term["label"]
    except:
        return None


def get_ontology_mappings(trait_name, filters, zooma_host):
    '''
    First get the URI, label from a selected source, confidence and source:
    http://snarf.ebi.ac.uk:8580/spot/zooma/v2/api/services/annotate?propertyValue=intellectual+disability
    Then the ontology label to replace the label from a source:
    http://www.ebi.ac.uk/ols/api/terms?iri=http%3A%2F%2Fwww.ebi.ac.uk%2Fefo%2FEFO_0003847
    '''
    url = build_zooma_query(trait_name, filters, zooma_host)
    zooma_response = request_retry_helper(zooma_query_helper, 4, url)

    if zooma_response is None:
        return None

    mappings = get_mappings_for_trait(zooma_response)

    for mapping in mappings:
        for zooma_entry in mapping.zooma_entry_list:
            label = get_ontology_label_from_ols(zooma_entry.uri)
            # If no label is returned (shouldn't really happen) keep the existing one
            if label is not None:
                zooma_entry.ontology_label = label
            else:
                print(
                    "Couldn't retrieve ontology label from OLS for trait '{}'".format(
                        trait_name))

            uri_is_current_and_in_efo = is_current_and_in_efo(zooma_entry.uri)
            if not uri_is_current_and_in_efo:
                uri_is_in_efo = is_in_efo(zooma_entry.uri)
                zooma_entry.in_efo = uri_is_in_efo
            else:
                zooma_entry.in_efo = uri_is_current_and_in_efo
                zooma_entry.is_current = uri_is_current_and_in_efo

    return mappings


def build_zooma_query(trait_name, filters, zooma_host):
    url = "{}/spot/zooma/v2/api/services/annotate?propertyValue={}".format(zooma_host, trait_name)
    url_filters = [
                    "required:[{}]".format(filters["required"]),
                    "ontologies:[{}]".format(filters["ontologies"]),
                    "preferred:[{}]".format(filters["preferred"])
                  ]
    url += "&filter={}".format(",".join(url_filters))
    return url


def get_mappings_for_trait(zooma_response):
    mappings = []
    for result in zooma_response:
        # uri_list = ",".join(result["semanticTags"])
        uris = result["semanticTags"]
        zooma_label = result["annotatedProperty"]["propertyValue"]
        confidence = result["confidence"]
        source_name = result["derivedFrom"]["provenance"]["source"]["name"]
        mappings.append(ZoomaMapping(uris, zooma_label, confidence, source_name))
    return mappings


@lru_cache(maxsize=8192)
def get_ontology_label_from_ols(ontology_uri):
    url = build_ols_query(ontology_uri)
    label = request_retry_helper(ols_query_helper, 4, url)
    return label


def build_ols_query(ontology_uri):
    url = "http://www.ebi.ac.uk/ols/api/terms?iri={}".format(ontology_uri)
    return url


##
# OxO functions
##


URI_DB_TO_DB_DICT = {
    "ordo": "Orphanet",
    "omim": "OMIM",
    "efo": "EFO",
    "mesh": "MeSH",
    "obo": "HP"
}

NON_NUMERIC_RE = re.compile(r'[^\d]+')


def uri_to_oxo_format(uri):
    if not any(x in uri.lower() for x in URI_DB_TO_DB_DICT.keys()):
        return None
    uri = uri.rstrip("/")
    uri_list = uri.split("/")
    id_ = NON_NUMERIC_RE.sub("", uri_list[-1])
    db = URI_DB_TO_DB_DICT[uri_list[-2].lower()]
    return "{}:{}".format(db, id_)


def uris_to_oxo_format(uri_list):
    oxo_id_list = []
    for uri in uri_list:
        oxo_id = uri_to_oxo_format(uri)
        if oxo_id is not None:
            oxo_id_list.append(oxo_id)
    return oxo_id_list


def build_oxo_payload(id_list, target_list, distance):
    payload = {}
    payload["ids"] = id_list
    payload["mappingTarget"] = target_list
    payload["distance"] = distance
    return payload


def oxo_query_helper(url, payload):
    try:
        json_response = requests.post(url, data=payload).json()
        return json_response
    except json.decoder.JSONDecodeError as e:
        return None


def oxo_request_retry_helper(retry_count, url, id_list, target_list, distance):
    payload = build_oxo_payload(id_list, target_list, distance)
    for retry_num in range(retry_count):
        return_value = oxo_query_helper(url, payload)
        if return_value is not None:
            return return_value
        print("attempt {}: failed running function oxo_query_helper with url {}".format(retry_num, url))
    print("error on last attempt, skipping")
    return None


def get_oxo_results_from_response(oxo_response):
    oxo_result_list = []
    results = oxo_response["_embedded"]["searchResults"]
    for result in results:
        if len(result["mappingResponseList"]) == 0:
            continue
        query_id = result["queryId"]
        label = result["label"]
        curie = result["curie"]
        oxo_result = OxOResult(query_id, label, curie)
        for mapping_response in result["mappingResponseList"]:
            mapping_label = mapping_response["label"]
            mapping_curie = mapping_response["curie"]
            mapping_distance = mapping_response["distance"]
            oxo_mapping = OxOMapping(mapping_label, mapping_curie, mapping_distance, query_id)

            uri = str(oxo_mapping.uri)

            ontology_label = get_ontology_label_from_ols(uri)
            if ontology_label is not None:
                oxo_mapping.ontology_label = ontology_label

            uri_is_current_and_in_efo = is_current_and_in_efo(uri)
            if not uri_is_current_and_in_efo:
                uri_is_in_efo = is_in_efo(uri)
                oxo_mapping.in_efo = uri_is_in_efo
            else:
                oxo_mapping.in_efo = uri_is_current_and_in_efo
                oxo_mapping.is_current = uri_is_current_and_in_efo

            oxo_result.oxo_mapping_list.append(oxo_mapping)

        oxo_result_list.append(oxo_result)

    return oxo_result_list


def get_oxo_results(id_list, target_list, distance):
    url = "http://www.ebi.ac.uk/spot/oxo/api/search?size=5000"
    oxo_response = oxo_request_retry_helper(4, url, id_list, target_list, distance)

    if oxo_response is None:
        return None

    oxo_results = get_oxo_results_from_response(oxo_response)
    return oxo_results


# http://www.ebi.ac.uk/spot/oxo/api/search
# ?size=1000
# POST
# payload:
# {"ids":["EFO:0001360","DOID:162","OMIM:180200","MESH:D009202","UBERON_0002107","HP_0005978"],"mappingTarget":[],"distance":"2"}

# Targets:
# "Orphanet"
# "efo"
# "hp"


##
# OLS functions
##


def double_encode_uri(uri):
    return urllib.parse.quote(urllib.parse.quote(uri, safe=""), safe="")


def ols_efo_query(uri):
    double_encoded_uri = double_encode_uri(uri)
    return requests.get(
        "http://www.ebi.ac.uk/ols/api/ontologies/efo/terms/{}".format(double_encoded_uri))


@lru_cache(maxsize=8192)
def is_current_and_in_efo(uri):
    response = ols_efo_query(uri)
    if response.status_code != 200:
        return False
    response_json = response.json()
    return not response_json["is_obsolete"]


def is_in_efo(uri):
    response = ols_efo_query(uri)
    return response.status_code == 200


##
# Argparser
##


class ArgParser:

    def __init__(self, argv):
        description = """
                Script for running terms through Zooma, retrieving mapped uri, label from OLS,
                confidence of the mapping, and source of the mapping.
                """
        parser = argparse.ArgumentParser(description=description)

        parser.add_argument("-i", dest="input_filepath", required=True, help="ClinVar json file. One record per line.")
        parser.add_argument("-o", dest="output_mappings_filepath", required=True, help="path to output file for mappings")
        parser.add_argument("-c", dest="output_curation_filepath", required=True, help="path to output file for curation")
        parser.add_argument("-n", dest="ontologies", default="efo,ordo,hp", help="ontologies to use in query")
        parser.add_argument("-r", dest="required", default="cttv,eva-clinvar,gwas", help="data sources to use in query.")
        parser.add_argument("-p", dest="preferred", default="eva-clinvar,cttv,gwas", help="preference for data sources, with preferred data source first.")
        parser.add_argument("-z", dest="zooma_host", default="https://www.ebi.ac.uk", help="the host to use for querying zooma")
        parser.add_argument("-t", dest="oxo_target_list", default="Orphanet,efo,hp", help="target ontologies to use with OxO")
        parser.add_argument("-d", dest="oxo_distance", default=3, help="distance to use to query OxO.")

        args = parser.parse_args(args=argv[1:])

        self.input_filepath = args.input_filepath
        self.output_mappings_filepath = args.output_mappings_filepath
        self.output_curation_filepath = args.output_curation_filepath

        self.filters = {"ontologies": args.ontologies, "required": args.required, "preferred": args.preferred}

        self.zooma_host = args.zooma_host
        self.oxo_target_list = [target.strip() for target in args.oxo_target_list.split(",")]
        self.oxo_distance = args.oxo_distance


if __name__ == '__main__':
    main()
