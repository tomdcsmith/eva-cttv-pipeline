import argparse
import csv
import sys

from clinvar_jsons_shared_lib import clinvar_jsons, get_traits_from_json, has_allowed_clinical_significance


def main():
    parser = ArgParser(sys.argv)

    trait_dict = {}
    for clinvar_json in clinvar_jsons(parser.infile_path):
        if not is_allowed_clinical_significance(clinvar_json):
            continue
        trait_dict = get_traits_from_json(clinvar_json, trait_dict)

    with open(parser.outfile_path, "w", newline="") as outfile:
        writer = csv.writer(outfile, delimiter="\t")
        for trait in trait_dict.values():
            writer.writerow([trait.name, trait.xref_string, trait.count])


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
