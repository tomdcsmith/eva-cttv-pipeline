def map_efo(trait_2_efo_dict, name_list):
    trait_string = name_list[0].lower()
    if trait_string in trait_2_efo_dict:
        ontology_id = trait_2_efo_dict[trait_string][0]
        ontology_label = trait_2_efo_dict[trait_string][1]
        return trait_string, ontology_id, ontology_label
    else:
        for trait in name_list[1:]:
            trait_string = trait.lower()
            if trait_string in trait_2_efo_dict:
                ontology_id = trait_2_efo_dict[trait_string][0]
                ontology_label = trait_2_efo_dict[trait_string][1]
                return trait_string, ontology_id, ontology_label

    return None, None, None


class Trait:
    def __init__(self, clinvar_trait_name_list, trait_counter, trait_2_efo_dict):
        self.trait_counter = trait_counter  # number of trait for record
        self.clinvar_name, self.ontology_id, self.ontology_label = map_efo(trait_2_efo_dict,
                                                                           clinvar_trait_name_list)
