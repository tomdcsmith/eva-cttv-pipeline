import unittest

import clinvar_to_evidence_strings


class TestGetUnmappedUrl(unittest.TestCase):
    def test_orphanet(self):
        urls = ["www.test.com/Orphanet_1234", "bio.co.uk/Orphanet_431"]
        for url in urls:
            self.assertEquals(MAIN.get_unmapped_url(url), "http://purl.bioontology.org/ORDO/" + url.split("/")[-1])
