import argparse
import sys


from .parse_trait_names import get_trait_names, clinvar_jsons
from .extract_pathogenic_and_likely_pathogenic_variants import has_allowed_clinical_significance


def main():
    parser = ArgParser(sys.argv)

    for clinvar_json in clinvar_jsons(parser.infile_path):
        if not has_allowed_clinical_significance(clinvar_json):
            continue
        clinva_acc = get_clinvar_acc(clinvar_json)

        var_trait_names = get_trait_names(clinvar_json)


def get_clinvar_acc(clinvar_json):
    return clinvar_json["clinvarSet"]["referenceClinVarAssertion"]["clinVarAccession"]["acc"]


class ArgParser:
    def __init__(self, argv):
        description = """
                Script for extracting the trait names of ClinVar records from a file with a list
                of CellBase, ClinVar JSONs, and the number of traits with this trait name.
                """
        parser = argparse.ArgumentParser(description=description)

        parser.add_argument("-i", dest="infile_path", required=True, help="Path to a file containing one CellBase ClinVar JSON per line. This should usually contain just the ClinVar records that have allowed clinical significances'")
        parser.add_argument("-o", dest="outfile_path", required=True, help="Path to file to output trait names in the zooma-accepted format")

        args = parser.parse_args(args=argv[1:])

        self.infile_path = args.infile_path
        self.outfile_path = args.outfile_path


if __name__ == '__main__':
    main()
