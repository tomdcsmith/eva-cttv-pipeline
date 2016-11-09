def map_efo(trait_2_efo_dict, name_list):
    """
    Function accepts a dict with mappings from clinvar trait names to tuples containing ontology
    ids and labels, and a list with clinvar trait names for one trait. Returns the clinvar trait
    name (in lowercase) that is earliest in the list and has a mapping in the dict, along with the
    ontology id and label.

    :param trait_2_efo_dict: dict with clinvar trait names as keys, and tuples as values. The first
    element of each tuple is an ontology id (uri), the second element is the ontology label for
    that id.
    :param name_list: List of clinvar trait names for the trait, with the preferred trait name (if
    one exists) in the first element.
    :return: The clinvar name, ontology id, and ontology label, for the clinvar name in the
    name_list with the lowest element that is in the trait_2_efo_dict. If none of the names in the
    name_list are keys in the trait_2_efo_dict then None is returned in each position.
    """
    trait_string = name_list[0].lower()  # Try first element, which should be preferred name if it exists
    if trait_string in trait_2_efo_dict:
        ontology_id = trait_2_efo_dict[trait_string][0]
        ontology_label = trait_2_efo_dict[trait_string][1]
        return trait_string, ontology_id, ontology_label
    else:
        for trait in name_list[1:]:  # Otherwise cycle through the list, the first name that is in
            trait_string = trait.lower()  # the dict is returned along with the ontology id and
            if trait_string in trait_2_efo_dict:  # label from the dict
                ontology_id = trait_2_efo_dict[trait_string][0]
                ontology_label = trait_2_efo_dict[trait_string][1]
                return trait_string, ontology_id, ontology_label

    return None, None, None  # If none of the names in the list are keys in the dict then None is returned


class Trait:
    def __init__(self, clinvar_trait_name_list, trait_counter, trait_2_efo_dict):
        self.trait_counter = trait_counter  # number of trait for record
        self.clinvar_name, self.ontology_id, self.ontology_label = map_efo(trait_2_efo_dict,
                                                                           clinvar_trait_name_list)
