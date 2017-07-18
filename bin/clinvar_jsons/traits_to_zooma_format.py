import argparse
import itertools
import sys
from time import gmtime, strftime

from parse_trait_names import clinvar_jsons, get_traits_from_json
from extract_pathogenic_and_likely_pathogenic_variants import has_allowed_clinical_significance


DATE = strftime("%d/%m/%y %H:%M", gmtime())


def main():
    parser = ArgParser(sys.argv)

    with open(parser.outfile_path, "wt") as outfile:
        outfile.write("STUDY\tBIOENTITY\tPROPERTY_TYPE\tPROPERTY_VALUE\tSEMANTIC_TAG\tANNOTATOR\tANNOTATION_DATE\n")
        for clinvar_json in clinvar_jsons(parser.infile_path):
            if not has_allowed_clinical_significance(clinvar_json):
                continue
            process_clinvar_json(clinvar_json, outfile)


def process_clinvar_json(clinvar_json, outfile):
    clinvar_acc = get_clinvar_acc(clinvar_json)
    variant_id_list = get_variant_ids(clinvar_json)

    trait_dict = {}
    trait_dict = get_traits_from_json(clinvar_json, trait_dict)

    for variant_id, trait_dict_item in itertools.product(variant_id_list,
                                                          trait_dict.items()):
        trait_dict_key, trait_dict_value = trait_dict_item
        trait_name = trait_dict_key
        if trait_name == "not provided":
            continue
        for xref in trait_dict_value.xref_set:
            if xref.status.lower() != "current" or xref.db.lower() not in OntologyUri.db_to_uri_dict:
                continue
            ontology_uri = OntologyUri(xref.id_, xref.db)

            write_zooma_record(clinvar_acc, variant_id, trait_name, ontology_uri, DATE, outfile)


def write_zooma_record(clinvar_acc, variant_id, trait_name, ontology_uri, date, outfile):
    zooma_output_list = [clinvar_acc,
                         variant_id,
                         "disease",
                         trait_name,
                         str(ontology_uri),
                         "clinvar-xrefs",
                         date]
    outfile.write("\t".join(zooma_output_list) + "\n")


def get_variant_ids(clinvar_json):
    rs_id_list = get_rs_ids(clinvar_json)
    sv_id_list = get_sv_ids(clinvar_json)
    return rs_id_list + sv_id_list


def get_rs_ids(clinvar_json):
    rs_id_list = []
    for measure in clinvar_json["clinvarSet"]['referenceClinVarAssertion']["measureSet"]["measure"]:
        if "xref" in measure:
            for xref in measure["xref"]:
                if xref["db"].lower() == "dbsnp":
                    rs_id_list.append("rs{}".format(xref["id"]))
    return rs_id_list


def get_sv_ids(clinvar_json):
    sv_id_list = []
    for measure in clinvar_json["clinvarSet"]['referenceClinVarAssertion']["measureSet"]["measure"]:
        if "xref" in measure:
            for xref in measure["xref"]:
                if xref["db"].lower() == "dbvar" and xref["id"].lower()[:3] in ("nsv", "esv"):
                    sv_id_list.append(xref["id"])
    return sv_id_list


def get_clinvar_acc(clinvar_json):
    return clinvar_json["clinvarSet"]["referenceClinVarAssertion"]["clinVarAccession"]["acc"]


class OntologyUri:
    db_to_uri_dict = {
        "orphanet": "http://www.orpha.net/ORDO/Orphanet_{}",
        "omim": "http://identifiers.org/omim/{}",
        "efo": "http://www.ebi.ac.uk/efo/{}",
        "mesh": "http://identifiers.org/mesh/{}",
        "medgen": "http://identifiers.org/medgen/{}",
        "human phenotype ontology": "http://purl.obolibrary.org/obo/HP_{}"
    }

    def __init__(self, id_, db):
        self.id_ = id_
        self.db = db
        if self.db.lower() == "human phenotype ontology":
            self.uri = self.db_to_uri_dict[self.db.lower()].format(self.id_[3:])
        else:
            self.uri = self.db_to_uri_dict[self.db.lower()].format(self.id_)

    def __str__(self):
        return self.uri


class ArgParser:
    def __init__(self, argv):
        description = """
                Script for extracting the trait names of ClinVar records from a file with a list
                of CellBase, ClinVar JSONs, and the number of traits with this trait name.
                """
        parser = argparse.ArgumentParser(description=description)

        parser.add_argument("-i", dest="infile_path", required=True, help="Path to a file containing one CellBase ClinVar JSON per line.'")
        parser.add_argument("-o", dest="outfile_path", required=True, help="Path to file to output trait names in the zooma-accepted format")

        args = parser.parse_args(args=argv[1:])

        self.infile_path = args.infile_path
        self.outfile_path = args.outfile_path


if __name__ == '__main__':
    main()
