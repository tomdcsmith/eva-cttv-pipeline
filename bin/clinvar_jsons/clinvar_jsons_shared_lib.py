import gzip
import json
from collections import namedtuple


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


def clinvar_jsons(filepath):
    with gzip.open(filepath, "rt") as f:
        for line in f:
            line = line.rstrip()
            yield json.loads(line)


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


def has_allowed_clinical_significance(clinvar_json):
    if "description" in clinvar_json["clinvarSet"]["referenceClinVarAssertion"]["clinicalSignificance"]:
        if clinvar_json["clinvarSet"]["referenceClinVarAssertion"]["clinicalSignificance"]["description"].lower() \
                in ["pathogenic", "likely pathogenic", "protective", "association", "risk_factor", "affects", "drug response"]:
            return True
    return False


