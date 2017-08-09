class OntologyEntry:
    """
    A representation of an ontology term mapped to using this pipeline. Includes a uri and a label
    which are the two details needed in a finished mapping, in addition to the ClinVar trait name.
    """
    def __init__(self, uri, label):
        self.uri = uri
        self.label = label

    def __eq__(self, other):
        if not isinstance(other, type(self)):
            return False
        if self.uri != other.uri or self.label != other.label:
            return False
        return True

    def __hash__(self):
        return hash((self.uri, self.label))


class Trait:
    """
    Object to hold data for one trait name. Including the number of ClinVar record's traits it
    appears in, any Zooma and OxO mappings, and any mappings which are ready to be output.
    """
    def __init__(self, name, frequency):
        self.name = name
        self.frequency = frequency
        self.zooma_mapping_list = []
        self.oxo_xref_list = []
        self.finished_mapping_set = set()

    @property
    def is_finished(self):
        """
        Return boolean which confirms whether the trait has finished mappings ready to be output
        """
        return len(self.finished_mapping_set) > 0

    def process_zooma_mappings(self):
        """
        Check whether any Zooma mappings can be output as a finished ontology mapping.
        Put any finished mappings in finished_mapping_set
        """
        for mapping in self.zooma_mapping_list:
            if mapping.confidence.lower() != "high":
                continue

            for entry in mapping.zooma_entry_list:
                if entry.in_efo and entry.is_current:
                    ontology_entry = OntologyEntry(entry.uri, entry.ontology_label)
                    self.finished_mapping_set.add(ontology_entry)

    def process_oxo_mappings(self):
        """
        Check whether any OxO mappings can be output as a finished ontology mapping.
        Put any finished mappings in finished_mapping_set
        """
        for result in self.oxo_xref_list:
            for mapping in result.oxo_mapping_list:
                if mapping.in_efo and mapping.is_current and mapping.distance == 1:
                    uri = str(mapping.uri)
                    ontology_label = mapping.ontology_label
                    ontology_entry = OntologyEntry(uri, ontology_label)
                    self.finished_mapping_set.add(ontology_entry)
