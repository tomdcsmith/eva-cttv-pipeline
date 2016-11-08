def map_efo(trait_2_efo_dict, name_list):
    trait_string = name_list[0].lower()
    if trait_string in trait_2_efo_dict:
        for efo_trait in trait_2_efo_dict[trait_string]:
            # First element in trait_list mus always be the "Preferred" trait name
            return trait_string, efo_trait
    else:
        for trait in name_list[1:]:
            trait_string = trait.lower()
            if trait_string in trait_2_efo_dict:
                for efo_trait in trait_2_efo_dict[trait_string]:
                    # First element in trait_list mus always be the "Preferred" trait name
                    return trait_string, efo_trait

    return None, None


class Trait:
    def __init__(self, clinvar_trait_name_list, trait_counter, trait_2_efo_dict):
        self.trait_counter = trait_counter  # number of trait for record
        self.clinvar_name, self.ontology_id = map_efo(trait_2_efo_dict,
                                                      clinvar_trait_name_list)
