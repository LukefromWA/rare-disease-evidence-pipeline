import unittest
from d3_merge_mechanistic_ledger import (
    combine_sources, combine_claims, mechanism_relation
)

def source(snapshot="first.xml", doi="10.0000/synthetic"):
    return {
        "source_id": "PMIDTEST",
        "doi": doi,
        "title": "Synthetic source",
        "snapshot": snapshot,
        "snapshot_sha256": "a" * 64,
    }

class MergeTests(unittest.TestCase):
    def test_source_deduplication_preserves_snapshots(self):
        sources, snapshots = combine_sources(
            [source()], [source(snapshot="second.xml")]
        )
        self.assertEqual(len(sources), 1)
        self.assertEqual(len(snapshots), 2)

    def test_conflicting_doi_fails(self):
        with self.assertRaises(ValueError):
            combine_sources(
                [source()], [source(doi="10.0000/different")]
            )

    def test_claim_deduplication(self):
        row = {"claim_id": "synthetic_claim", "claim_text": "Synthetic"}
        self.assertEqual(len(combine_claims([row], [row.copy()])), 1)

    def test_conflicting_claim_fails(self):
        with self.assertRaises(ValueError):
            combine_claims(
                [{"claim_id": "x", "claim_text": "First"}],
                [{"claim_id": "x", "claim_text": "Different"}],
            )

    def test_g418_is_comparator_only(self):
        self.assertEqual(
            mechanism_relation("G418"),
            "biochemical_comparator_only_not_candidate_nomination",
        )

    def test_mechanism_does_not_claim_functional_rescue(self):
        self.assertEqual(
            mechanism_relation("ataluren"),
            "translation_mechanism_rationale_not_case_functional_validation",
        )

if __name__ == "__main__":
    unittest.main()
