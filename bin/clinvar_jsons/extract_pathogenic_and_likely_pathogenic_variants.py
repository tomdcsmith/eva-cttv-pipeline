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
    for clinvar_assertion in clinvar_json["clinVarAssertion"]:
        if "description" in clinvar_assertion["clinicalSignificance"]:
            for description in clinvar_assertion["clinicalSignificance"]["description"]:
                clin_sigs.add(description)
    return len(clin_sigs.intersection({"Pathogenic", "Likely pathogenic"})) > 0


class ArgParser:
    def __init__(self, argv):
        parser = argparse.ArgumentParser()

        parser.add_argument("-i", dest="infile_path", required=True)
        parser.add_argument("-o", dest="outfile_path", required=True)

        args = parser.parse_args(args=argv[1:])

        self.infile_path = args.infile_path
        self.outfile_path = args.outfile_path


if __name__ == "__main__":
    main()
