import gzip
import json


def clinvar_jsons(filepath: str) -> dict:
    """
    Yields a dict object, parsed from a file with one json per line of ClinVar records

    :param filepath: String giving the path to a gzipped file containing jsons of ClinVar records,
                     one per line.
    :return: Yields a dictionary for each json.
    """
    with gzip.open(filepath, "rt") as f:
        for line in f:
            line = line.rstrip()
            yield json.loads(line)


def get_trait_names(clinvar_json: dict) -> list:
    """
    Given a dictionary for a ClinVar record parse the strings for the trait name of each trait in
    the record, prioritising trait names which are specified as being "Preferred". Return a list of
    there strings.

    :param clinvar_json: A dict object containing a ClinVar record
    :return: List of strings, the elements being one trait name for each trait in the input record.
    """
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
            if name['elementValue']['type'] == 'Preferred':
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


def parse_trait_names(filepath: str) -> list:
    """
    For a file containing ClinVar records in the format of one json per line, return a list of the
    trait names for the records in the file.

    :param filepath: String giving the path to a gzipped file containing jsons of ClinVar records,
                     one per line.
    :return: List (important that it is a list and not a set, it is used to calculate the
             frequencies) containing the trait names of the traits of the records in the file
             containing ClinVar records.
    """
    trait_name_list = []
    for clinvar_json in clinvar_jsons(filepath):
        new_trait_names = get_trait_names(clinvar_json)
        trait_name_list.extend(new_trait_names)
    return trait_name_list
