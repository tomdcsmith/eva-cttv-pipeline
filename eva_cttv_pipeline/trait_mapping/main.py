import argparse
from collections import Counter
import csv
import progressbar
import sys

from eva_cttv_pipeline.trait_mapping.output import output_trait
from eva_cttv_pipeline.trait_mapping.oxo import get_oxo_results
from eva_cttv_pipeline.trait_mapping.oxo import uris_to_oxo_format
from eva_cttv_pipeline.trait_mapping.trait import Trait
from eva_cttv_pipeline.trait_mapping.trait_names_parsing import parse_trait_names
from eva_cttv_pipeline.trait_mapping.zooma import get_ontology_mappings


def get_uris_for_oxo(zooma_mapping_list: list) -> set:
    """
    For a list of Zooma mappings return a list of uris for the mappings in that list with a high
    confidence but not in EFO.

    :param zooma_mapping_list: List with elements of class ZoomaMapping
    :return: set of uris from high confidence Zooma mappings, for which to query OxO
    """
    uri_set = set()
    for mapping in zooma_mapping_list:
        # Only use high confidence Zooma mappings for querying OxO
        if mapping.confidence.lower() != "high":
            continue
        uri_set.update([entry.uri for entry in mapping.zooma_entry_list])
    return uri_set


def process_trait(trait: Trait, filters: dict, zooma_host: str, oxo_target_list: list,
                  oxo_distance: int) -> Trait:
    """
    Process a single trait. Find any mappings in Zooma. If there are no high confidence Zooma
    mappings that are in EFO then query OxO with any high confidence mappings not in EFO.

    :param trait: A trait of class Trait.
    :param filters: A dictionary of filters to use for querying Zooma.
    :param zooma_host: A string with the hostname to use for querying Zooma
    :param oxo_target_list: A list of strings, each being an OxO ID for an ontology. Used to specify
                            which ontologies should be queried using OxO.
    :param oxo_distance: int specifying the maximum number of steps to use to query OxO. i.e. OxO's
                         "distance" parameter.
    :return: The original trait after querying Zooma and possibly OxO, with any results found.
    """
    zooma_mappings = get_ontology_mappings(trait.name, filters, zooma_host)
    trait.zooma_mapping_list = zooma_mappings
    trait.process_zooma_mappings()
    if (trait.is_finished
            or len(trait.zooma_mapping_list) == 0
            or any([entry.is_current
                    for mapping in trait.zooma_mapping_list
                    for entry in mapping.zooma_entry_list])):
        return trait
    uris_for_oxo_set = get_uris_for_oxo(trait.zooma_mapping_list)
    if len(uris_for_oxo_set) == 0:
        return trait
    oxo_input_id_list = uris_to_oxo_format(uris_for_oxo_set)
    oxo_result_list = get_oxo_results(oxo_input_id_list, oxo_target_list, oxo_distance)
    trait.oxo_xref_list = oxo_result_list
    trait.process_oxo_mappings()

    return trait


def main():
    parser = ArgParser(sys.argv)

    trait_names_list = parse_trait_names(parser.input_filepath)
    trait_names_counter = Counter(trait_names_list)

    with open(parser.output_mappings_filepath, "w", newline='') as mapping_file, \
            open(parser.output_curation_filepath, "wt") as curation_file:
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


class ArgParser:

    def __init__(self, argv):
        description = """
                Script for running terms through Zooma, retrieving mapped uri, label from OLS,
                confidence of the mapping, and source of the mapping.
                """
        parser = argparse.ArgumentParser(description=description)

        parser.add_argument("-i", dest="input_filepath", required=True,
                            help="ClinVar json file. One record per line.")
        parser.add_argument("-o", dest="output_mappings_filepath", required=True,
                            help="path to output file for mappings")
        parser.add_argument("-c", dest="output_curation_filepath", required=True,
                            help="path to output file for curation")
        parser.add_argument("-n", dest="ontologies", default="efo,ordo,hp",
                            help="ontologies to use in query")
        parser.add_argument("-r", dest="required", default="cttv,eva-clinvar,gwas",
                            help="data sources to use in query.")
        parser.add_argument("-p", dest="preferred", default="eva-clinvar,cttv,gwas",
                            help="preference for data sources, with preferred data source first.")
        parser.add_argument("-z", dest="zooma_host", default="https://www.ebi.ac.uk",
                            help="the host to use for querying zooma")
        parser.add_argument("-t", dest="oxo_target_list", default="Orphanet,efo,hp",
                            help="target ontologies to use with OxO")
        parser.add_argument("-d", dest="oxo_distance", default=3,
                            help="distance to use to query OxO.")

        args = parser.parse_args(args=argv[1:])

        self.input_filepath = args.input_filepath
        self.output_mappings_filepath = args.output_mappings_filepath
        self.output_curation_filepath = args.output_curation_filepath

        self.filters = {"ontologies": args.ontologies,
                        "required": args.required,
                        "preferred": args.preferred}

        self.zooma_host = args.zooma_host
        self.oxo_target_list = [target.strip() for target in args.oxo_target_list.split(",")]
        self.oxo_distance = args.oxo_distance


if __name__ == '__main__':
    main()
