import unittest
from d3_rtdetective_adopt import exact_record

def fixture():
    bridge = {
        "chrom": "chr1",
        "bed_chromStart": "99",
        "bed_chromEnd": "102",
    }
    status = {
        "chrom": "chr1", "target_start": "99", "target_end": "102",
        "status": "unique_exact_ascending_coordinate_match",
        "n_exact_records": "1",
        "n_nonforward_records_near_target": "0",
    }
    record = {
        "chrom": "chr1", "target_start": "99", "target_end": "102",
        "coordinate_order": "ascending",
        "raw_start": "99", "raw_end": "102",
        "author_score_raw": "1.234567890",
        "source_line_number": "20",
    }
    return bridge, status, record

class AdoptionTests(unittest.TestCase):
    def test_accept_and_preserve_raw_score(self):
        bridge, status, record = fixture()
        result, _ = exact_record(bridge, [status], [record])
        self.assertEqual(result["author_score_raw"], "1.234567890")

    def test_missing_status(self):
        bridge, _, record = fixture()
        result, reason = exact_record(bridge, [], [record])
        self.assertIsNone(result)
        self.assertEqual(reason, "no_cached_target_status")

    def test_duplicate_status(self):
        bridge, status, record = fixture()
        result, _ = exact_record(bridge, [status, status], [record])
        self.assertIsNone(result)

    def test_nonforward_ambiguity(self):
        bridge, status, record = fixture()
        status["n_nonforward_records_near_target"] = "1"
        result, _ = exact_record(bridge, [status], [record])
        self.assertIsNone(result)

    def test_duplicate_exact_records(self):
        bridge, status, record = fixture()
        result, _ = exact_record(bridge, [status], [record, record])
        self.assertIsNone(result)

    def test_descending_not_normalized(self):
        bridge, status, record = fixture()
        record["coordinate_order"] = "descending"
        result, _ = exact_record(bridge, [status], [record])
        self.assertIsNone(result)

    def test_interval_disagreement(self):
        bridge, status, record = fixture()
        record["raw_start"] = "100"
        result, _ = exact_record(bridge, [status], [record])
        self.assertIsNone(result)

    def test_nonfinite_score(self):
        bridge, status, record = fixture()
        record["author_score_raw"] = "NaN"
        result, _ = exact_record(bridge, [status], [record])
        self.assertIsNone(result)

    def test_wrong_target(self):
        bridge, status, record = fixture()
        status["target_start"] = "199"
        result, _ = exact_record(bridge, [status], [record])
        self.assertIsNone(result)

    def test_invalid_source_line(self):
        bridge, status, record = fixture()
        record["source_line_number"] = "0"
        result, _ = exact_record(bridge, [status], [record])
        self.assertIsNone(result)

if __name__ == "__main__":
    unittest.main()
