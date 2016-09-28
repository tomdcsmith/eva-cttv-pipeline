import argparse
import requests
import sys


def main():
    parser = ArgParser(sys.argv)

    trait_name_tuples = get_trait_name_tuples_from_file(parser.input_filepath)
    parse_and_write_output(trait_name_tuples, parser.output_filepath, parser.filters)


def parse_and_write_output(trait_name_tuples, output_filepath, filters):
    with open(output_filepath, "wt") as output_file:
        output_file.write("#{}\n".format(filters))
        for trait_name_tuple in trait_name_tuples:
            trait_name = trait_name_tuple[0]
            zooma_response = query_zooma(trait_name, filters)
            if zooma_response == 0:
                continue
            output_list = zooma_response_to_list(zooma_response, trait_name_tuple)

            output_file.write("\t".join(output_list) + "\n")


def get_trait_name_tuples_from_file(filepath):
    trait_name_tuples = []
    with open(filepath) as f:
        for line in f:
            line = line.rstrip()
            line_list = line.split("\t")
            trait_name = line_list[0]
            freq = line_list[1]
            trait_name_tuples.append((trait_name, freq))
    return trait_name_tuples


def form_zooma_query(trait_name, filters):
    url_filters = []
    url = "http://snarf.ebi.ac.uk:8580/spot/zooma/v2/api/services/annotate?propertyValue={}".format(trait_name)
    url_filters.append("required:[{}]".format(filters["required"]))
    url_filters.append("ontologies:[{}]".format(filters["ontologies"]))
    url_filters.append("preferred:[{}]".format(filters["preferred"]))
    url += "&filter={}".format(",".join(url_filters))

    return url


def query_zooma(trait_name, filters):
    url = form_zooma_query(trait_name, filters)
    r = requests.get(url)
    try:
        json_response = r.json()
    except Exception as e:
        # print("\t".join([trait_name, e, r]))
        return 0
    return json_response


def result_to_tuple(result):
    label = result["derivedFrom"]["annotatedProperty"]["propertyValue"]
    uri = ",".join(result["semanticTags"])
    confidence = result["confidence"]
    source_name = result["derivedFrom"]["provenance"]["source"]["name"]

    result_tuple = (label, uri, confidence, source_name)

    return result_tuple


def zooma_response_to_list(zooma_response, trait_name_tuple):
    trait_name = trait_name_tuple[0]
    freq = trait_name_tuple[1]
    output_list = [trait_name, freq]
    for result in zooma_response:
        output_list.extend(result_to_tuple(result))

    return output_list


class ArgParser:

    def __init__(self, argv):
        parser = argparse.ArgumentParser()

        parser.add_argument("-i", dest="input_filepath", required=True, help="path to input file, with trait names in first column, number of variants the trait name appears in in the second column. delimeted using tab")
        parser.add_argument("-o", dest="output_filepath", required=True, help="path to output file (not just the directory). outputs a file with a header (line starting with \"#\") which shows the filters used. then the first column is trait name, then number of variants for the trait, then zooma label, uri(s), confidence, source. these zooma columns repeat when there are multiple mappings.")
        parser.add_argument("-n", dest="ontologies", default="efo,hp,ordo", help="ontologies to use in query")
        parser.add_argument("-r", dest="required", default="bmb-wp7,cttv,eva-clinvar,sysmicro,atlas,uniprot,gwas,ebisc", help="data sources to use in query.")
        parser.add_argument("-p", dest="preferred", default="efo,hp,ordo,cttv,eva-clinvar,gwas,atlas", help="preference for data sources, with preferred data source first.")

        args = parser.parse_args(args=argv[1:])

        self.input_filepath = args.input_filepath
        self.output_filepath = args.output_filepath

        self.filters = {"ontologies": args.ontologies, "required": args.required, "preferred": args.preferred}


if __name__ == "__main__":
    main()
