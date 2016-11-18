import argparse
import json
import gzip
import sys
from collections import defaultdict


def main():
    parser = ArgParser(sys.argv)

    with gzip.open(parser.outfile_path, "wt") as outfile:
        for clinvar_json in clinvar_jsons(parser.infile_path):
            if is_path_or_likely_path(clinvar_json):
                outfile.write(json.dumps(clinvar_json) + "\n")


def clinvar_jsons(filepath):
    with gzip.open(filepath, "rt") as f:
        for line in f:
            line = line.rstrip()
            yield json.loads(line)


def is_path_or_likely_path(clinvar_json):
    clin_sigs = set()
    for clinvar_assertion in clinvar_json["clinvarSet"]["clinVarAssertion"]:
        if "description" in clinvar_assertion["clinicalSignificance"]:
            for description in clinvar_assertion["clinicalSignificance"]["description"]:
                clin_sigs.add(description)
    return len(clin_sigs.intersection({"Pathogenic", "Likely pathogenic"})) > 0


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
