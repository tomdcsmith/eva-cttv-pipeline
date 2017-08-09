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
    :return: The clinvar name, and mapping to an ontology id and term, for the clinvar name in the
    name_list with the lowest element that is in the trait_2_efo_dict. If none of the names in the
    name_list are keys in the trait_2_efo_dict then None is returned in each position.
    """
    for trait in name_list:
        trait_string = trait.lower()
        if trait_string in trait_2_efo_dict:
            return trait_string, trait_2_efo_dict[trait_string]

    return None, None  # If none of the names in the list are keys in the dict then None is returned


class Trait:
    def __init__(self, clinvar_name, ontology_id, ontology_label, trait_counter):
        self.trait_counter = trait_counter  # number of trait for record
        self.clinvar_name = clinvar_name
        self.ontology_id = ontology_id
        self.ontology_label = ontology_label

    def __str__(self):
        return "clinvar name: {} ontology id: {} ontology label: {} trait_counter: {}".format(self.clinvar_name,
                                                                                              self.ontology_id,
                                                                                              self.ontology_label,
                                                                                              self.trait_counter)

    def __eq__(self, other):
        if isinstance(other, self.__class__):
            return self.__dict__ == other.__dict__
        else:
            return False
