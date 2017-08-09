import unittest

import eva_cttv_pipeline.trait_mapping.oxo as oxo


class TestUriToOxoFormat(unittest.TestCase):
    def test_ordo(self):
        self.assertEqual(oxo.uri_to_oxo_format("http://www.orpha.net/ORDO/Orphanet_140162"),
                         "Orphanet:140162")

    def test_omim(self):
        self.assertEqual(oxo.uri_to_oxo_format("http://purl.obolibrary.org/obo/OMIM_314580"),
                         "OMIM:314580")

    def test_efo(self):
        self.assertEqual(oxo.uri_to_oxo_format("http://www.ebi.ac.uk/efo/EFO_0000313"),
                         "EFO:0000313")

    def test_mesh(self):
        self.assertEqual(oxo.uri_to_oxo_format("http://purl.obolibrary.org/obo/MESH_D002277"),
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
        oxo_response = {'_embedded': {'searchResults': [{'curie': 'HP:0001537', '_links': {'self': {'href': 'http://www.ebi.ac.uk/spot/oxo/api/terms/HP:0001537'}, 'mappings': {'href': 'http://www.ebi.ac.uk/spot/oxo/api/mappings?fromId=HP:0001537'}}, 'label': 'Umbilical hernia', 'querySource': None, 'queryId': 'HP:0001537', 'mappingResponseList': [{'curie': 'Orphanet:660', 'targetPrefix': 'Orphanet', 'distance': 3, 'sourcePrefixes': ['ONTONEO', 'Orphanet', 'DOID', 'UMLS'], 'label': 'Omphalocele'}, {'curie': 'Orphanet:3164', 'targetPrefix': 'Orphanet', 'distance': 3, 'sourcePrefixes': ['ONTONEO', 'EFO', 'Orphanet', 'DOID'], 'label': 'Omphalocele syndrome, Shprintzen-Goldberg type'}]}]}, 'page': {'totalPages': 1, 'size': 1000, 'number': 0, 'totalElements': 1}, '_links': {'self': {'href': 'http://www.ebi.ac.uk/spot/oxo/api/search'}}}
        expected_oxo_result = oxo.OxOResult("HP:0001537", "Umbilical hernia", "HP:0001537")
        expected_oxo_result.oxo_mapping_list = [oxo.OxOMapping("Omphalocele", "Orphanet:660", 3, "HP:0001537"),
                                                oxo.OxOMapping("Omphalocele syndrome, Shprintzen-Goldberg type", "Orphanet:3164", 3, "HP:0001537")]
        expected_oxo_results = [expected_oxo_result]

        self.assertEqual(oxo.get_oxo_results_from_response(oxo_response),
                         expected_oxo_results)







