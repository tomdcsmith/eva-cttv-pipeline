import csv

from eva_cttv_pipeline.trait_mapping.trait import Trait


def output_trait_mapping(trait: Trait, mapping_writer: csv.writer):
    """
    Write any finished ontology mappings for a trait to a csv file writer.

    :param trait: A trait with finished ontology mappings in finished_mapping_set
    :param mapping_writer: A csv.writer to write the finished mappings
    """
    for ontology_entry in trait.finished_mapping_set:
        mapping_writer.writerow([trait.name, ontology_entry.uri, ontology_entry.label])


def get_zooma_mappings_for_curation(trait: Trait) -> list:
    """Sorted in reverse so the highest ranked zooma entries are shown first"""
    mapping_list = []
    for zooma_result in trait.zooma_result_list:
        for zooma_mapping in zooma_result.mapping_list:
            if zooma_mapping.in_efo and zooma_mapping.is_current:
                mapping_list.append(zooma_mapping)
    mapping_list.sort(reverse=True)
    return mapping_list


# TODO best way to unite this and the above very similar method?
def get_oxo_mappings_for_curation(trait: Trait) -> list:
    """Sorted in reverse so the highest ranked oxo mappings are shown first"""
    oxo_mapping_list = []
    for oxo_result in trait.oxo_xref_list:
        for oxo_mapping in oxo_result.oxo_mapping_list:
            if oxo_mapping.in_efo and oxo_mapping.is_current:
                oxo_mapping_list.append(oxo_mapping)
    oxo_mapping_list.sort(reverse=True)
    return oxo_mapping_list


def output_for_curation(trait: Trait, curation_writer: csv.writer):
    """
    Write any non-finished Zooma or OxO mappings of a trait to a file for manual curation.
    Also outputs traits without any ontology mappings.

    :param trait: A Trait with no finished ontology mappings in finished_mapping_set
    :param curation_writer: A csv.writer to write non-finished ontology mappings for manual curation
    """
    output_row = [trait.name, trait.frequency]

    zooma_mapping_list = get_zooma_mappings_for_curation(trait)

    for zooma_mapping in zooma_mapping_list:
        cell = [zooma_mapping.uri, zooma_mapping.ontology_label, str(zooma_mapping.confidence),
                zooma_mapping.source]
        output_row.append("|".join(cell))

    oxo_mapping_list = get_oxo_mappings_for_curation(trait)

    for oxo_mapping in oxo_mapping_list:
        cell = [str(oxo_mapping.uri), oxo_mapping.ontology_label, str(oxo_mapping.distance),
                oxo_mapping.query_id]
        output_row.append("|".join(cell))

    curation_writer.writerow(output_row)


def output_trait(trait: Trait, mapping_writer: csv.writer, curation_writer: csv.writer):
    """
    Output finished ontology mappings of a trait, or non-finished mappings (if any) for curation.

    :param trait: A trait which has been used to query Zooma and possibly OxO.
    :param mapping_writer: A csv.writer to write the finished mappings
    :param curation_writer: A csv.writer to write non-finished ontology mappings for manual curation
    """
    if trait.is_finished:
        output_trait_mapping(trait, mapping_writer)
    else:
        output_for_curation(trait, curation_writer)
