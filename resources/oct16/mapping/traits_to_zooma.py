import argparse
import requests
import sys


def main():
    parser = ArgParser(sys.argv)

    trait_name_tuples = get_trait_name_tuples_from_file(parser.input_filepath)
    parse_and_write_output(trait_name_tuples, parser.output_filepath)


def parse_and_write_output(trait_name_tuples, output_filepath):
    with open(output_filepath, "wt") as output_file:
        for trait_name_tuple in trait_name_tuples:
            trait_name = trait_name_tuple[0]
            zooma_response = query_zooma(trait_name)
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


# def form_zooma_query(property_value, ontologies, required, preferred):
def form_zooma_query(trait_name):
    url = "http://snarf.ebi.ac.uk:8580/spot/zooma/v2/api/services/annotate?propertyValue={}&filter=ontologies:[efo,hp,ordo],required:[none]".format(trait_name)
    return url


# def query_zooma(property_value, ontologies=None, required=None, preferred=None):
def query_zooma(trait_name):
    url = form_zooma_query(trait_name)
    r = requests.get(url)
    return r.json()


def result_to_tuple(result):
    label = result["annotatedProperty"]["propertyValue"]
    if len(result["semanticTags"]) > 1:
        print(result)
        sys.exit()
    uri = result["semanticTags"][0]
    confidence = result["confidence"]

    result_tuple = (label, uri, confidence)

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

        parser.add_argument("-i", dest="input_filepath")
        parser.add_argument("-o", dest="output_filepath")

        args = parser.parse_args(args=argv[1:])

        self.input_filepath = args.input_filepath
        self.output_filepath = args.output_filepath


if __name__ == "__main__":
    main()
