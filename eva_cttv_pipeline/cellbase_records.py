import codecs
import json
import urllib.request
import urllib.error

from eva_cttv_pipeline import config, utilities


class CellbaseRecords:

    """Assists in the requesting and iteration of clinvar cellbase records."""

    def __init__(self, json_file):
        """

        :param json_file: Path to a file containing a list of json strings of the Clinvar records
        from Cellbase, one per line. This can be used to potentially save time since requests to
        Cellbase are subsequently not needed.
        """
        self.json_file = json_file

    def __iter__(self):
        for record in self.__each_line_in_file():
            yield record

    def __each_line_in_file(self):
        with utilities.open_file(self.json_file, "rt") as f:
            for line in f:
                yield json.loads(line.rstrip())
