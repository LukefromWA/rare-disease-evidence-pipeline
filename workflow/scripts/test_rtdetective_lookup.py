import unittest
from d3_rtdetective_lookup import scan_lines
HEADER=b"#chrom chromStart chromEnd score\n"
TARGET=("chr1",99,102)

def scan(*rows): return scan_lines([HEADER]+[r.encode()+b"\n" for r in rows],[TARGET])
class Tests(unittest.TestCase):
    def test_exact_preserves_score_and_line(self):
        statuses,exact,_=scan("chr1 99 102 1.234567890")
        self.assertEqual(statuses[0]["status"],"unique_exact_ascending_coordinate_match")
        self.assertEqual(exact[0]["author_score_raw"],"1.234567890")
        self.assertEqual(exact[0]["source_line_number"],2)
    def test_duplicate_not_aggregated(self):
        statuses,exact,_=scan(*(["chr1 99 102 1.0"]*2))
        self.assertEqual(statuses[0]["n_exact_records"],2); self.assertEqual(len(exact),2)
    def test_descending_not_adopted(self):
        statuses,exact,_=scan("chr1 102 99 1.0")
        self.assertEqual(exact,[]); self.assertEqual(statuses[0]["n_nonforward_records_near_target"],1)
    def test_forward_plus_descending_is_blocked(self):
        statuses,_,_=scan("chr1 99 102 1.0","chr1 102 99 2.0")
        self.assertNotEqual(statuses[0]["status"],"unique_exact_ascending_coordinate_match")
    def test_wrong_chromosome(self):
        statuses,exact,_=scan("chr2 99 102 1.0")
        self.assertEqual(exact,[]); self.assertEqual(statuses[0]["n_exact_records"],0)
    def test_off_by_one(self):
        statuses,exact,_=scan("chr1 100 103 1.0")
        self.assertEqual(exact,[])
    def test_nan(self):
        with self.assertRaises(ValueError): scan("chr1 99 102 NaN")
    def test_bad_row(self):
        with self.assertRaises(ValueError): scan("chr1 bad 102 1.0")
    def test_missing_header(self):
        with self.assertRaises(ValueError): scan_lines([b"chr1 1 4 1.0"],[TARGET])
    def test_adjacent_descending_not_overlap(self):
        statuses,_,_=scan("chr1 105 102 1.0")
        self.assertEqual(statuses[0]["n_nonforward_records_near_target"],0)
    def test_empty_targets_still_scans(self):
        statuses,exact,count=scan_lines([HEADER,b"chr2 1 4 2\n"],[])
        self.assertEqual((statuses,exact,count),([],[],1))
    def test_no_exact_hit_reported(self):
        statuses,_,_=scan("chr1 10 13 1.0")
        self.assertEqual(statuses[0]["status"],"no_exact_ascending_coordinate_match")
if __name__=="__main__": unittest.main()
