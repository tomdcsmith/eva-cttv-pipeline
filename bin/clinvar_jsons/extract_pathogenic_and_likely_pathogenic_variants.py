import argparse
import json
import gzip
import sys

from clinvar_jsons_shared_lib import clinvar_jsons, get_traits_from_json, has_allowed_clinical_significance


def main():
    parser = ArgParser(sys.argv)

    with gzip.open(parser.outfile_path, "wt") as outfile:
        for clinvar_json in clinvar_jsons(parser.infile_path):
            if has_allowed_clinical_significance(clinvar_json):
                outfile.write(json.dumps(clinvar_json) + "\n")


class ArgParser:
    def __init__(self, argv):
        description = """
        Script for extracting the ClinVar records that have clinical significance of 'pathogenic'
        and 'likely pathogenic' from a file with a list of CellBase ClinVar JSONs. The output is
        in the same format as the input.
        """
        parser = argparse.ArgumentParser(description=description)

        parser.add_argument("-i", dest="infile_path", required=True, help="a path to a file containing one CellBase ClinVar JSON per line")
        parser.add_argument("-o", dest="outfile_path", required=True, help="a path to the output file which will contain the JSONs from the input file, but only those that are specified as having a pathogenic or likely pathogenic clincar significance")

        args = parser.parse_args(args=argv[1:])

        self.infile_path = args.infile_path
        self.outfile_path = args.outfile_path


if __name__ == "__main__":
    main()
