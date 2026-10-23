import unittest
from d3_claim_ledger import make_claim, normalize, stable_id, hypotheses_from_bridge

def fixture():
    seed = {
        "seed_key": "synthetic",
        "pmid": "TEST",
        "compound": "synthetic_compound",
        "evidence_class": "synthetic_experiment",
        "source_gene": "SYN",
        "population_or_model": "synthetic",
        "outcome": "synthetic_endpoint",
        "claim_text": "Synthetic extraction claim.",
        "anchor": {"kind": "abstract", "quotes": ["signal increased"]},
        "statistics": {},
        "limitations": "Synthetic test only.",
        "hypothesis_relation": "not_clinical_evidence",
    }
    metadata = {"TEST": {"abstract": "The signal increased in this model."}}
    return seed, metadata

class LedgerTests(unittest.TestCase):
    def test_abstract_anchor(self):
        seed, metadata = fixture()
        claim = make_claim(seed, metadata, [])
        self.assertEqual(claim["source_id"], "PMIDTEST")

    def test_missing_abstract_anchor_fails(self):
        seed, metadata = fixture()
        seed["anchor"]["quotes"] = ["unsupported statement"]
        with self.assertRaises(ValueError):
            make_claim(seed, metadata, [])

    def test_missing_publication_fails(self):
        seed, _ = fixture()
        with self.assertRaises(ValueError):
            make_claim(seed, {}, [])

    def test_normalization(self):
        self.assertEqual(normalize("−1.0   to  2.0"), "-1.0 to 2.0")

    def test_stable_identifier(self):
        self.assertEqual(
            stable_id("X", {"a": 1, "b": 2}),
            stable_id("X", {"b": 2, "a": 1}),
        )

    def test_pending_review_not_automatic_adjudication(self):
        seed, metadata = fixture()
        claim = make_claim(seed, metadata, [])
        self.assertEqual(
            claim["review_status"],
            "source_anchors_checked_adjudication_pending",
        )

    def test_table_anchor(self):
        seed, metadata = fixture()
        seed["anchor"] = {
            "kind": "table_cells", "table_id": "T1",
            "cells": [{"row": 1, "column": 1, "expected": "0.5"}],
        }
        rows = [{
            "pmid": "TEST", "table_id": "T1", "row_number": "1",
            "cells_json": '["effect", "0.5"]',
        }]
        self.assertIn(
            "table=T1", make_claim(seed, metadata, rows)["source_locator"]
        )

    def test_wrong_table_value_fails(self):
        seed, metadata = fixture()
        seed["anchor"] = {
            "kind": "table_cells", "table_id": "T1",
            "cells": [{"row": 1, "column": 0, "expected": "0.5"}],
        }
        rows = [{
            "pmid": "TEST", "table_id": "T1", "row_number": "1",
            "cells_json": '["0.6"]',
        }]
        with self.assertRaises(ValueError):
            make_claim(seed, metadata, rows)


class BridgeLinkageTests(unittest.TestCase):
    def synthetic_row(self):
        return {
            "hypothesis_id": "synthetic_hypothesis",
            "gene": "SYN",
            "variant": "synthetic_allele",
            "q_j_raw": "0.2",
            "selecting_views": "synthetic_view",
            "phase_status": "unknown",
            "haplotype_scope": "synthetic",
            "NMD_status": "not_assessed",
            "upstream_route_status": "eligible_for_context_verification",
            "bridge_status": "missing_context",
        }

    def test_hypothesis_retained_without_prediction_or_context(self):
        result = hypotheses_from_bridge([self.synthetic_row()])
        self.assertIn("synthetic_hypothesis", result)
        self.assertEqual(
            result["synthetic_hypothesis"]["q_j_raw"], "0.2"
        )

    def test_nonapplicable_hypothesis_not_linked(self):
        row = self.synthetic_row()
        row["upstream_route_status"] = "not_applicable"
        self.assertEqual(hypotheses_from_bridge([row]), {})

    def test_duplicate_hypothesis_rejected(self):
        row = self.synthetic_row()
        with self.assertRaises(ValueError):
            hypotheses_from_bridge([row, row.copy()])


class TableColumnRegressionTests(unittest.TestCase):
    def synthetic_table(self):
        seed, metadata = fixture()
        seed["anchor"] = {
            "kind": "table_cells",
            "table_id": "T2",
            "cells": [
                {"row": 10, "column": 2, "expected": "0.75"},
                {"row": 10, "column": 6, "expected": "0.5"},
            ],
        }
        rows = [{
            "pmid": "TEST",
            "table_id": "T2",
            "row_number": "10",
            "cells_json": '["A/B", "", "0.75c", "", "", "", "0.5"]',
        }]
        return seed, metadata, rows

    def test_empty_cells_preserve_column_positions_and_footnote(self):
        seed, metadata, rows = self.synthetic_table()
        claim = make_claim(seed, metadata, rows)
        self.assertIn("0.75c", claim["source_excerpt"])

    def test_adjacent_empty_column_is_rejected(self):
        seed, metadata, rows = self.synthetic_table()
        seed["anchor"]["cells"][0]["column"] = 3
        with self.assertRaises(ValueError):
            make_claim(seed, metadata, rows)

if __name__ == "__main__":
    unittest.main()
