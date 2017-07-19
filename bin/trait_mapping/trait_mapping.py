import argparse
import gzip
import json
import sys
import urllib

import requests
from collections import Counter


##
# Classes
##


class ZoomaMapping:
    def __init__(self, uri_list, zooma_label, confidence, source):
        self.uri_list = uri_list
        self.zooma_label = zooma_label
        self.label_list = ["" for _ in range(len(uri_list))]
        self.uri_in_efo_list = [False for _ in range(len(uri_list))]
        self.uri_current_list = [False for _ in range(len(uri_list))]
        self.confidence = confidence
        self.source = source


class OntologyEntry:
    def __init__(self, uri, label):
        self.uri = uri
        self.label = label


class Trait:
    def __init__(self, name, frequency):
        self.name = name
        self.frequency = frequency
        self.zooma_mapping_list = []
        self.oxo_xref_list = []
        self.finished_mapping_list = []

    @property
    def is_finished(self):
        return len(self.finished_mapping_list) > 0

    def process_zooma_mappings(self):
        for mapping in self.zooma_mapping_list:
            for uri, uri_is_in_efo, uri_is_current, label in zip(mapping.uri_list,
                                                                 mapping.uri_in_efo_list,
                                                                 mapping.uri_current_list,
                                                                 mapping.label_list):
                if uri_is_in_efo and uri_is_current:
                    ontology_entry = OntologyEntry(uri, label)
                    self.finished_mapping_list.append(ontology_entry)


##
# Main method
##


def main():
    parser = ArgParser(sys.argv)

    trait_names_list = parse_trait_names(parser.input_filepath)
    trait_names_counter = Counter(trait_names_list)

    for trait_name, freq in trait_names_counter.items():
        trait = Trait(trait_name, freq)
        process_trait_name(trait, parser.filters, parser.zooma_host)


def process_trait_name(trait, filters, zooma_host):
    zooma_mappings = get_ontology_mappings(trait.name, filters, zooma_host)
    trait.zooma_mapping_list = zooma_mappings
    trait.process_zooma_mappings()
    if trait.is_finished:
        return




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
        for idx, uri in enumerate(mapping.uri_list):
            label = get_ontology_label_from_ols(uri)
            # If no label is returned (shouldn't really happen) keep the existing one
            if label is not None:
                mapping.label_list[idx] = label
            else:
                print(
                    "Couldn't retrieve ontology label from OLS for trait '{}', will use the one from Zooma".format(
                        trait_name))

            uri_is_current_and_in_efo = is_current_and_in_efo(uri)
            if not uri_is_current_and_in_efo:
                uri_is_in_efo = is_in_efo(uri)
                mapping.uri_in_efo_list[idx] = uri_is_in_efo
            else:
                mapping.uri_in_efo_list[idx] = uri_is_current_and_in_efo
                mapping.uri_current_list[idx] = uri_is_current_and_in_efo

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


def get_ontology_label_from_ols(uri_mapping):
    url = build_ols_query(uri_mapping)
    json_response_1 = request_retry_helper(ols_query_helper, 4, url)
    return json_response_1


def build_ols_query(ontology_uri):
    url = "http://www.ebi.ac.uk/ols/api/terms?iri={}".format(ontology_uri)
    return url



##
# OLS functions
##


def double_encode_uri(uri):
    return urllib.parse.quote(urllib.parse.quote(uri, safe=""), safe="")


def ols_efo_query(uri):
    double_encoded_uri = double_encode_uri(uri)
    return requests.get(
        "http://www.ebi.ac.uk/ols/api/ontologies/efo/terms/{}".format(double_encoded_uri))


def is_current_and_in_efo(uri):
    response = ols_efo_query(uri)
    if response.status_code == 400:
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
        parser.add_argument("-o", dest="output_filepath", required=True, help="path to output file")
        parser.add_argument("-n", dest="ontologies", default="efo,ordo,hp", help="ontologies to use in query")
        parser.add_argument("-r", dest="required", default="cttv,eva-clinvar,gwas", help="data sources to use in query.")
        parser.add_argument("-p", dest="preferred", default="eva-clinvar,cttv,gwas", help="preference for data sources, with preferred data source first.")
        parser.add_argument("-z", dest="zooma_host", default="https://www.ebi.ac.uk", help="the host to use for querying zooma")

        args = parser.parse_args(args=argv[1:])

        self.input_filepath = args.input_filepath
        self.output_filepath = args.output_filepath

        self.filters = {"ontologies": args.ontologies, "required": args.required, "preferred": args.preferred}

        self.zooma_host = args.zooma_host


if __name__ == '__main__':
    main()
