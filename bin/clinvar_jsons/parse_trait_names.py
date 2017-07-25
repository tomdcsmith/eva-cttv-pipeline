import argparse
import csv
import json
import gzip
import sys
from collections import defaultdict, namedtuple


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


def is_allowed_clinical_significance(clinvar_json):
    allowed_clinical_significance_list = ["pathogenic", "likely pathogenic", "protective",
                                          "association", "risk_factor", "affects", "drug response"]
    if "description" in clinvar_json["clinvarSet"]["referenceClinVarAssertion"]["clinicalSignificance"]:
        if clinvar_json["clinvarSet"]["referenceClinVarAssertion"]["clinicalSignificance"]["description"].lower() in allowed_clinical_significance_list:
            return True
    return False


def get_traits_from_json(clinvar_json, trait_dict):

    # This if-else block is due to the change in the format of the CellBase JSON that holds the
    # ClinVar data. Originally "clinvarSet" was the top level, but this level was removed and
    # referenceClinVarAssertion is now the top level.
    if "clinvarSet" in clinvar_json:
        trait_set = clinvar_json["clinvarSet"]["referenceClinVarAssertion"]["traitSet"]
    else:
        trait_set = clinvar_json["referenceClinVarAssertion"]["traitSet"]
    for trait_doc in trait_set['trait']:
        preferred_trait_name = None
        non_preferred_names = []
        for name in trait_doc['name']:
            # First trait name in the list will always be the "Preferred" one
            if name['elementValue']['type'].lower() == 'preferred':
                preferred_trait_name = name['elementValue']['value']
                break
            elif name['elementValue']['type'] in ["EFO URL", "EFO id", "EFO name"]:
                continue  # if the trait name not originally from clinvar
            else:
                non_preferred_names.append(name['elementValue']['value'])
        if preferred_trait_name is None:
            preferred_trait_name = non_preferred_names[0]

        if preferred_trait_name in trait_dict:
            trait = trait_dict[preferred_trait_name]
            trait.count += 1
        else:
            trait = Trait(preferred_trait_name)

        if "xref" in trait_doc:
            for xref in trait_doc["xref"]:
                trait.xref_set.add(TraitXref(xref["db"], xref["id"], xref["status"]))

        trait_dict[preferred_trait_name] = trait

    return trait_dict


def clinvar_jsons(filepath):
    with gzip.open(filepath, "rt") as f:
        for line in f:
            line = line.rstrip()
            yield json.loads(line)


class Trait:
    def __init__(self, name):
        self.name = name
        self.xref_set = set()
        self.count = 1

    # def __eq__(self, other):
    #     if isinstance(other, Trait):
    #         return self.name == other.name and set(self.xref_set) == set(other.xref_set)
    #     else:
    #         return False
    #
    # def __hash__(self):
    #     return hash((self.name, self.xref_set))
    #
    # def __ne__(self, other):
    #     return not (self == other)

    @property
    def xref_string(self):
        return "|".join(["{}/{}".format(xref.db, xref.id_)
                         for xref in self.xref_set if xref.status.lower() == "current"])

    def __str__(self):
        return "{}\t{}\t{}".format(self.name, self.xref_string, self.count)


TraitXref = namedtuple("TraitXref", ["db", "id_", "status"])


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
