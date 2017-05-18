import argparse
import json
import gzip
import sys
from collections import defaultdict


def main():
    parser = ArgParser(sys.argv)

    trait_names = defaultdict(int)

    for clinvar_json in clinvar_jsons(parser.infile_path):
        var_trait_names = get_trait_names(clinvar_json)
        for name in var_trait_names:
            trait_names[name] += 1

    with gzip.open(parser.outfile_path, "wt") as outfile:
        for name, count in trait_names.items():
            outfile.write(name + "\t" + str(count) + "\n")


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
            if name['elementValue']['type'].lower() == 'preferred':
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


def clinvar_jsons(filepath):
    with gzip.open(filepath, "rt") as f:
        for line in f:
            line = line.rstrip()
            yield json.loads(line)


class ArgParser:
    def __init__(self, argv):
        description = """
                Script for extracting the trait names of ClinVar records from a file with a list
                of CellBase, ClinVar JSONs, and the number of traits with this trait name.
                """
        parser = argparse.ArgumentParser(description=description)

        parser.add_argument("-i", dest="infile_path", required=True, help="Path to a file containing one CellBase ClinVar JSON per line. This should usually contain just the ClinVar records that have clinical significances of 'pathogenic' and 'likely pathogenic'")
        parser.add_argument("-o", dest="outfile_path", required=True, help="Path to file to output trait names. This is a tab separated file with the trait name in the first column, and the number of traits (not records, since some records have multiple traits, although that count would not differ much) that have this trait name.")

        args = parser.parse_args(args=argv[1:])

        self.infile_path = args.infile_path
        self.outfile_path = args.outfile_path


if __name__ == "__main__":
    main()
