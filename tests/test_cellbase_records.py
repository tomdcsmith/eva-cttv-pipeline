import unittest

from eva_cttv_pipeline import cellbase_records, config


class CellbaseRecordsTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.cb_records = cellbase_records.CellbaseRecords()

    def test_get_curr_result_list(self):
        curr_response = self.cb_records._CellbaseRecords__get_curr_response(0)
        curr_result_list = self.cb_records._CellbaseRecords__get_curr_result_list(curr_response)
        self.assertEqual(len(curr_result_list), config.BATCH_SIZE)

        curr_response = self.cb_records._CellbaseRecords__get_curr_response(0)
        curr_response_test = self.cb_records._CellbaseRecords__get_curr_response(curr_response['numTotalResults'])
        curr_result_list_test = self.cb_records._CellbaseRecords__get_curr_result_list(curr_response_test)
        self.assertEqual(len(curr_result_list_test), 0)

        curr_response = self.cb_records._CellbaseRecords__get_curr_response(999999)
        curr_result_list = self.cb_records._CellbaseRecords__get_curr_result_list(curr_response)
        self.assertEqual(len(curr_result_list), 0)

        curr_response = self.cb_records._CellbaseRecords__get_curr_response(0)
        len_to_expect = 20
        curr_response = self.cb_records._CellbaseRecords__get_curr_response(curr_response['numTotalResults'] - len_to_expect)
        curr_result_list = self.cb_records._CellbaseRecords__get_curr_result_list(curr_response)
        self.assertEqual(len(curr_result_list), len_to_expect)

    # def test_get_curr_result_lists(self):
    #     curr_response = clinvar_to_evidence_strings.get_curr_response(0)
    #     num_total_results = curr_response['numTotalResults']
    #     list_counter = 0
    #     for list in clinvar_to_evidence_strings.get_curr_result_lists():
    #         list_counter += 1
    #     pred_num_lists = math.ceil(num_total_results / config.BATCH_SIZE)
    #     self.assertEqual(pred_num_lists, list_counter)