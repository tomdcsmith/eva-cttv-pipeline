import argparse
import requests
import sys


class Trait:
    
    def __init__(self, name, frequency):
        self.name = name
        self.frequency = frequency
    
    def __str__(self):
        return "{}\t{}".format(self.name, self.frequency)


class OntologyMapping:
    
    def __init__(self, uri, label, confidence, source):
        self.uri = uri
        self.label = label
        self.confidence = confidence
        self.source = source
    
    def __str__(self):
        return "{}\t{}\t{}\t{}".format(self.label, self.uri, self.confidence, self.source)



def main():
    parser = ArgParser(sys.argv)

    # Read trait names from input file
    traits = read_traits(parser.input_filepath)
    
    with open(parser.output_filepath, "wt") as output_file:
        output_file.write("#{}\n".format(parser.filters))
        for trait in traits:
            output_file.write(str(trait))
            for mapping in get_ontology_mappings(trait, parser.filters):
                output_file.write("\t" + str(mapping))
            output_file.write("\n")
            


def read_traits(filepath):
    traits = []
    with open(filepath) as f:
        for line in f:
            line = line.rstrip()
            line_list = line.split("\t")
            name = line_list[0]
            frequency = line_list[1]
            traits.append(Trait(name, frequency))
    return traits


def write_ontology_mappings(ontology_mappings, filters, output_filepath):
    with open(output_filepath, "wt") as output_file:
        output_file.write("#{}\n".format(filters))
        for mapping in ontology_mappings:
            output_file.write("\t".join(mapping) + "\n")


def get_ontology_mappings(trait, filters):
    '''
    First get the URI, label from a selected source, confidence and source:
    http://snarf.ebi.ac.uk:8580/spot/zooma/v2/api/services/annotate?propertyValue=intellectual+disability
    Then the ontology label to replace the label from a source:
    http://www.ebi.ac.uk/ols/api/terms?iri=http%3A%2F%2Fwww.ebi.ac.uk%2Fefo%2FEFO_0003847
    '''
    url = build_zooma_query(trait.name, filters)
    json_response_1 = requests.get(url).json()
    
    mappings = get_mappings_for_trait(json_response_1, trait)
    
    for mapping in mappings:
        try:
            label = get_ontology_label_from_ols(mapping.uri)
            # If no label is returned (shouldn't really happen) keep the existing one
            if label:
                mapping.label = label
        except:
            print("Couldn't retrieve ontology label from OLS for trait '{}', will use the one from Zooma".format(trait.name))
    
    return mappings


def build_zooma_query(trait_name, filters):
    url_filters = []
    url_filters.append("required={}".format(filters["required"]))
    url_filters.append("ontologies={}".format(filters["ontologies"]))
    url_filters.append("preferred={}".format(filters["preferred"]))
    
    url = "http://snarf.ebi.ac.uk:8580/spot/zooma/v2/api/services/annotate?propertyValue={}&".format(trait_name)
    url += "&".join(url_filters)
    return url


def get_mappings_for_trait(zooma_response, trait):
    mappings = []
    for result in zooma_response:
        uri = ",".join(result["semanticTags"])
        label = result["annotatedProperty"]["propertyValue"]
        confidence = result["confidence"]
        source_name = result["derivedFrom"]["provenance"]["source"]["name"]
        mappings.append(OntologyMapping(uri, label, confidence, source_name))
    return mappings


def get_ontology_label_from_ols(uri_mapping):
    url = build_ols_query(uri_mapping)
    json_response = requests.get(url).json()
    for term in json_response["_embedded"]["terms"]:
        if term["is_defining_ontology"]:
            return term["label"]
    
    return None


def build_ols_query(ontology_uri):
    url = "http://www.ebi.ac.uk/ols/api/terms?iri={}".format(ontology_uri)
    return url



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
