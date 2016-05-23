import codecs
import json
import urllib.request

from eva_cttv_pipeline import config


class CellbaseRecords:
    def __init__(self, limit=config.BATCH_SIZE):
        self.limit = limit

    def __iter__(self):
        for curr_result_list in self.__get_curr_result_lists():
            for record in curr_result_list:
                yield record

    def __get_curr_response(self, skip):
        reader = codecs.getreader("utf-8")
        answer = urllib.request.urlopen('http://' + config.HOST + '/cellbase/webservices/rest/v3/hsapiens/feature/clinical/all?source=clinvar&skip=' + str(skip) + '&limit=' + str(self.limit))
        curr_response = json.load(reader(answer))['response'][0]
        return curr_response

    def __get_curr_result_list(self, curr_response):
        curr_result_list = curr_response['result']
        return curr_result_list

    def __get_curr_result_lists(self):
        skip = 0
        while True:
            curr_response = self.__get_curr_response(skip)
            curr_result_list = self.__get_curr_result_list(curr_response)
            if len(curr_result_list) == 0:
                break
            skip += config.BATCH_SIZE
            yield curr_result_list
