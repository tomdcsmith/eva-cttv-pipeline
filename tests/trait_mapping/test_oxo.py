import unittest

import requests_mock

import eva_cttv_pipeline.trait_mapping.oxo as oxo


class TestUriToOxoFormat(unittest.TestCase):
    def test_ordo(self):
        self.assertEqual(oxo.uri_to_oxo_format("http://www.orpha.net/ORDO/Orphanet_140162"),
                         "Orphanet:140162")

    def test_omim(self):
        self.assertEqual(oxo.uri_to_oxo_format("http://purl.obolibrary.org/obo/OMIM_314580"),
                         "OMIM:314580")
        self.assertEqual(oxo.uri_to_oxo_format("http://identifiers.org/omim/314580"),
                         "OMIM:314580")

    def test_efo(self):
        self.assertEqual(oxo.uri_to_oxo_format("http://www.ebi.ac.uk/efo/EFO_0000313"),
                         "EFO:0000313")

    def test_mesh(self):
        self.assertEqual(oxo.uri_to_oxo_format("http://purl.obolibrary.org/obo/MESH_D002277"),
                         "MeSH:D002277")
        self.assertEqual(oxo.uri_to_oxo_format("http://identifiers.org/mesh/D002277"),
                         "MeSH:D002277")

    def test_hp(self):
        self.assertEqual(oxo.uri_to_oxo_format("http://purl.obolibrary.org/obo/HP_0030731"),
                         "HP:0030731")

    def test_nonexistent(self):
        self.assertEqual(oxo.uri_to_oxo_format("not_a_real_uri"),
                         None)


class TestBuildOxoPayload(unittest.TestCase):
    def test_build_payload(self):
        id_list = ["OMIM:314580", "MeSH:D002277"]
        target_list = ["Orphanet", "efo", "hp"]
        distance = 3
        self.assertEqual(oxo.build_oxo_payload(id_list, target_list, distance),
                         {"ids": id_list, "mappingTarget": target_list, "distance": distance})


class TestGetOxoResultsFromResponse(unittest.TestCase):
    def test_get_oxo_results_from_response(self):
        url = "http://www.ebi.ac.uk/ols/api/terms?iri=http://www.orpha.net/ORDO/Orphanet_660"
        with requests_mock.mock() as m:
            m.get(url,
                  json = {'_embedded': {'terms': [{'has_children': False, '_links': {'has_inheritance': {'href': 'http://www.ebi.ac.uk/ols/api/ontologies/ordo/terms/http%253A%252F%252Fwww.orpha.net%252FORDO%252FOrphanet_660/http%253A%252F%252Fwww.orpha.net%252FORDO%252FOrphanet_C016'}, 'parents': {'href': 'http://www.ebi.ac.uk/ols/api/ontologies/ordo/terms/http%253A%252F%252Fwww.orpha.net%252FORDO%252FOrphanet_660/parents'}, 'self': {'href': 'http://www.ebi.ac.uk/ols/api/ontologies/ordo/terms/http%253A%252F%252Fwww.orpha.net%252FORDO%252FOrphanet_660'}, 'hierarchicalParents': {'href': 'http://www.ebi.ac.uk/ols/api/ontologies/ordo/terms/http%253A%252F%252Fwww.orpha.net%252FORDO%252FOrphanet_660/hierarchicalParents'}, 'part_of': {'href': 'http://www.ebi.ac.uk/ols/api/ontologies/ordo/terms/http%253A%252F%252Fwww.orpha.net%252FORDO%252FOrphanet_660/http%253A%252F%252Fpurl.obolibrary.org%252Fobo%252FBFO_0000050'}, 'hierarchicalAncestors': {'href': 'http://www.ebi.ac.uk/ols/api/ontologies/ordo/terms/http%253A%252F%252Fwww.orpha.net%252FORDO%252FOrphanet_660/hierarchicalAncestors'}, 'jstree': {'href': 'http://www.ebi.ac.uk/ols/api/ontologies/ordo/terms/http%253A%252F%252Fwww.orpha.net%252FORDO%252FOrphanet_660/jstree'}, 'has_AgeOfOnset': {'href': 'http://www.ebi.ac.uk/ols/api/ontologies/ordo/terms/http%253A%252F%252Fwww.orpha.net%252FORDO%252FOrphanet_660/http%253A%252F%252Fwww.orpha.net%252FORDO%252FOrphanet_C017'}, 'graph': {'href': 'http://www.ebi.ac.uk/ols/api/ontologies/ordo/terms/http%253A%252F%252Fwww.orpha.net%252FORDO%252FOrphanet_660/graph'}, 'ancestors': {'href': 'http://www.ebi.ac.uk/ols/api/ontologies/ordo/terms/http%253A%252F%252Fwww.orpha.net%252FORDO%252FOrphanet_660/ancestors'}}, 'iri': 'http://www.orpha.net/ORDO/Orphanet_660', 'ontology_prefix': 'ORDO', 'obo_synonym': None, 'is_defining_ontology': True, 'label': 'Omphalocele', 'is_root': False, 'obo_xref': [{'description': 'NTBT (narrower term maps to a broader term)', 'url': None, 'database': 'ICD-10', 'id': 'Q79.2'}, {'description': 'E (exact mapping (the terms and the concepts are equivalent))', 'url': None, 'database': 'MedDRA', 'id': '10030308'}, {'description': 'BTNT (broader term maps to a narrower term)', 'url': 'http://omim.org/entry/164750', 'database': 'OMIM', 'id': '164750'}, {'description': 'E (exact mapping (the terms and the concepts are equivalent))', 'url': None, 'database': 'UMLS', 'id': 'C0795690'}, {'description': 'BTNT (broader term maps to a narrower term)', 'url': 'http://omim.org/entry/310980', 'database': 'OMIM', 'id': '310980'}], 'description': ['Omphalocele is an embryopathy classified in the group of abdominal celosomias and is characterized by a large hernia of the abdominal wall, centered on the umbilical cord, in which the protruding viscera are protected by a sac.'], 'annotation': {'hasDbXref': ['UMLS:C0795690', 'OMIM:310980', 'ICD-10:Q79.2', 'MedDRA:10030308', 'OMIM:164750'], 'definition_citation': ['orphanet']}, 'ontology_name': 'ordo', 'in_subset': None, 'ontology_iri': 'http://www.orpha.net/ontology/orphanet.owl', 'term_replaced_by': None, 'is_obsolete': False, 'short_form': 'Orphanet_660', 'obo_definition_citation': None, 'synonyms': None, 'obo_id': 'Orphanet:660'}]}, 'page': {'totalPages': 1, 'size': 20, 'number': 0, 'totalElements': 1}, '_links': {'self': {'href': 'http://www.ebi.ac.uk/ols/api/terms'}}})

            oxo_response = {'_embedded': {'searchResults': [{'curie': 'HP:0001537', '_links': {'self': {'href': 'http://www.ebi.ac.uk/spot/oxo/api/terms/HP:0001537'}, 'mappings': {'href': 'http://www.ebi.ac.uk/spot/oxo/api/mappings?fromId=HP:0001537'}}, 'label': 'Umbilical hernia', 'querySource': None, 'queryId': 'HP:0001537', 'mappingResponseList': [{'curie': 'Orphanet:660', 'targetPrefix': 'Orphanet', 'distance': 3, 'sourcePrefixes': ['ONTONEO', 'Orphanet', 'DOID', 'UMLS'], 'label': 'Omphalocele'}, {'curie': 'Orphanet:3164', 'targetPrefix': 'Orphanet', 'distance': 3, 'sourcePrefixes': ['ONTONEO', 'EFO', 'Orphanet', 'DOID'], 'label': 'Omphalocele syndrome, Shprintzen-Goldberg type'}]}]}, 'page': {'totalPages': 1, 'size': 1000, 'number': 0, 'totalElements': 1}, '_links': {'self': {'href': 'http://www.ebi.ac.uk/spot/oxo/api/search'}}}
            expected_oxo_result = oxo.OxOResult("HP:0001537", "Umbilical hernia", "HP:0001537")
            expected_oxo_result.oxo_mapping_list = [oxo.OxOMapping("Omphalocele", "Orphanet:660", 3, "HP:0001537"),
                                                    oxo.OxOMapping("Omphalocele syndrome, Shprintzen-Goldberg type", "Orphanet:3164", 3, "HP:0001537")]
            expected_oxo_results = [expected_oxo_result]

            self.assertEqual(oxo.get_oxo_results_from_response(oxo_response),
                             expected_oxo_results)







