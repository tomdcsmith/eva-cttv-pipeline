import argparse
import gzip
import itertools
import json
import requests
import sys
from time import gmtime, strftime

import progressbar

from parse_trait_names import clinvar_jsons, get_traits_from_json
from extract_pathogenic_and_likely_pathogenic_variants import has_allowed_clinical_significance


DATE = strftime("%d/%m/%y %H:%M", gmtime())


def main():
    parser = ArgParser(sys.argv)

    with open(parser.outfile_path, "wt") as outfile:
        outfile.write("STUDY\tBIOENTITY\tPROPERTY_TYPE\tPROPERTY_VALUE\tSEMANTIC_TAG\tANNOTATOR\tANNOTATION_DATE\n")
        line_count = file_len(parser.infile_path)
        bar = progressbar.ProgressBar(max_value=line_count, widgets=[progressbar.AdaptiveETA(samples=1000)])
        is_zooma_mapping_dict = {}
        for clinvar_json in bar(clinvar_jsons(parser.infile_path)):
            if not has_allowed_clinical_significance(clinvar_json):
                continue
            process_clinvar_json(clinvar_json, outfile, parser.zooma_host, parser.filters,
                                 is_zooma_mapping_dict)


def process_clinvar_json(clinvar_json, outfile, zooma_host, filters, is_zooma_mapping_dict):
    clinvar_acc = get_clinvar_acc(clinvar_json)
    variant_id_list = get_variant_ids(clinvar_json)

    trait_dict = {}
    trait_dict = get_traits_from_json(clinvar_json, trait_dict)

    for variant_id, trait_dict_item in itertools.product(variant_id_list,
                                                         trait_dict.items()):
        trait_dict_key, trait_dict_value = trait_dict_item
        trait_name = trait_dict_key
        if trait_name == "not provided":
            continue

        if trait_name in is_zooma_mapping_dict:
            is_zooma_mapping = is_zooma_mapping_dict[trait_name]
        else:
            zooma_uri_set = get_zooma_uris(trait_name, zooma_host, filters)
            is_zooma_mapping = len(zooma_uri_set) > 0
            is_zooma_mapping_dict[trait_name] = is_zooma_mapping

        for xref in trait_dict_value.xref_set:
            if xref.status.lower() != "current" or xref.db.lower() not in OntologyUri.db_to_uri_dict:
                continue
            ontology_uri = OntologyUri(xref.id_, xref.db)

            if is_zooma_mapping:
                continue

            write_zooma_record(clinvar_acc, variant_id, trait_name, ontology_uri, DATE, outfile)


def write_zooma_record(clinvar_acc, variant_id, trait_name, ontology_uri, date, outfile):
    zooma_output_list = [clinvar_acc,
                         variant_id,
                         "disease",
                         trait_name,
                         str(ontology_uri),
                         "clinvar-xrefs",
                         date]
    outfile.write("\t".join(zooma_output_list) + "\n")


def get_variant_ids(clinvar_json):
    rs_id_list = get_rs_ids(clinvar_json)
    sv_id_list = get_sv_ids(clinvar_json)
    return rs_id_list + sv_id_list


def get_rs_ids(clinvar_json):
    rs_id_list = []
    for measure in clinvar_json["clinvarSet"]['referenceClinVarAssertion']["measureSet"]["measure"]:
        if "xref" in measure:
            for xref in measure["xref"]:
                if xref["db"].lower() == "dbsnp":
                    rs_id_list.append("rs{}".format(xref["id"]))
    return rs_id_list


def get_sv_ids(clinvar_json):
    sv_id_list = []
    for measure in clinvar_json["clinvarSet"]['referenceClinVarAssertion']["measureSet"]["measure"]:
        if "xref" in measure:
            for xref in measure["xref"]:
                if xref["db"].lower() == "dbvar" and xref["id"].lower()[:3] in ("nsv", "esv"):
                    sv_id_list.append(xref["id"])
    return sv_id_list


def get_clinvar_acc(clinvar_json):
    return clinvar_json["clinvarSet"]["referenceClinVarAssertion"]["clinVarAccession"]["acc"]


def get_zooma_uris(trait_name, zooma_host, filters):
    url = build_zooma_query(trait_name, filters, zooma_host)
    json_response = request_retry_helper(zooma_query_helper, 4, url)

    if json_response is None:
        return None

    uri_set = set()
    for result in json_response:
        if result["confidence"].lower() == "high":
            uri_set.update(result["semanticTags"])

    return uri_set


def build_zooma_query(trait_name, filters, zooma_host):
    url = "{}/spot/zooma/v2/api/services/annotate?propertyValue={}".format(zooma_host, trait_name)
    url_filters = [
                    "required:[{}]".format(filters["required"]),
                    "ontologies:[{}]".format(filters["ontologies"]),
                    "preferred:[{}]".format(filters["preferred"])
                  ]
    url += "&filter={}".format(",".join(url_filters))
    return url


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
        json_response = requests.get(url).json()
        return json_response
    except json.decoder.JSONDecodeError as e:
        return None


def open_file(file_path, mode):
    if file_path.endswith(".gz"):
        return gzip.open(file_path, mode)
    else:
        return open(file_path, mode)


def file_len(fname):
    with open_file(fname, "rt") as f:
        for i, l in enumerate(f):
            pass
    return i + 1


class OntologyUri:
    db_to_uri_dict = {
        "orphanet": "http://www.orpha.net/ORDO/Orphanet_{}",
        "omim": "http://identifiers.org/omim/{}",
        "efo": "http://www.ebi.ac.uk/efo/{}",
        "mesh": "http://identifiers.org/mesh/{}",
        "medgen": "http://identifiers.org/medgen/{}",
        "human phenotype ontology": "http://purl.obolibrary.org/obo/HP_{}"
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


class ArgParser:
    def __init__(self, argv):
        description = """
                Script for extracting the trait names of ClinVar records from a file with a list
                of CellBase, ClinVar JSONs, and the number of traits with this trait name.
                """
        parser = argparse.ArgumentParser(description=description)

        parser.add_argument("-i", dest="infile_path", required=True, help="Path to a file containing one CellBase ClinVar JSON per line.'")
        parser.add_argument("-o", dest="outfile_path", required=True, help="Path to file to output trait names in the zooma-accepted format")

        parser.add_argument("-n", dest="ontologies", default="efo,ordo,hp",
                            help="ontologies to use in query")
        parser.add_argument("-r", dest="required", default="cttv,eva-clinvar,gwas",
                            help="data sources to use in query.")
        parser.add_argument("-p", dest="preferred", default="eva-clinvar,cttv,gwas",
                            help="preference for data sources, with preferred data source first.")
        parser.add_argument("-z", dest="zooma_host", default="https://www.ebi.ac.uk",
                            help="the host to use for querying zooma")

        args = parser.parse_args(args=argv[1:])

        self.infile_path = args.infile_path
        self.outfile_path = args.outfile_path

        self.filters = {"ontologies": args.ontologies, "required": args.required,
                        "preferred": args.preferred}

        self.zooma_host = args.zooma_host


if __name__ == '__main__':
    main()
