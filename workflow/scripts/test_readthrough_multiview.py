import unittest
from d3_readthrough_multiview import select_routes

def scenarios():
    return [
        {
            "hypothesis_id": "A", "gene": "SYN_A",
            "variant": "synthetic_A", "q_j": "0.2",
            "consequence": "stop_gained",
        },
        {
            "hypothesis_id": "B", "gene": "SYN_B",
            "variant": "synthetic_B", "q_j": "0.8",
            "consequence": "missense_variant",
        },
    ]

def view(hypothesis_id, name, rank):
    return {
        "hypothesis_id": hypothesis_id,
        "gene": "SYN_" + hypothesis_id,
        "variant": "synthetic_" + hypothesis_id,
        "model_view": name,
        "view_status": "scored",
        "score_type": "synthetic_view_specific_quantity",
        "score": "0.123",
        "rank": str(rank),
        "source_field": "synthetic_source_field",
        "source_sha256": "a" * 64,
    }

def views():
    return [
        view("A", "stan", 2), view("B", "stan", 1),
        view("A", "a10_informed", 1),
        view("B", "a10_informed", 2),
    ]

class MultiviewTests(unittest.TestCase):
    def test_disagreement_preserved(self):
        routed, witnesses, _ = select_routes(
            scenarios(), views(), ["stan", "a10_informed"]
        )
        self.assertEqual(
            {row["hypothesis_id"] for row in routed
             if row["selected_by_policy"]},
            {"A", "B"},
        )
        self.assertEqual(len(witnesses), 2)

    def test_routing_separate_from_selection(self):
        routed, _, _ = select_routes(
            scenarios(), views(), ["stan", "a10_informed"]
        )
        statuses = {
            row["hypothesis_id"]: row["route_status"] for row in routed
        }
        self.assertEqual(statuses["A"], "eligible_for_context_verification")
        self.assertEqual(statuses["B"], "not_applicable")

    def test_original_q_preserved(self):
        routed, _, _ = select_routes(
            scenarios(), views(), ["stan", "a10_informed"]
        )
        self.assertEqual(
            {row["hypothesis_id"]: row["q_j_raw"] for row in routed},
            {"A": "0.2", "B": "0.8"},
        )

    def test_ties_retained(self):
        records = [view("A", "stan", 1), view("B", "stan", 1)]
        _, witnesses, _ = select_routes(scenarios(), records, ["stan"])
        self.assertEqual(len(witnesses), 2)

    def test_no_promotion_when_rank1_absent(self):
        records = [view("A", "stan", 2), view("B", "stan", 3)]
        routed, witnesses, coverage = select_routes(
            scenarios(), records, ["stan"]
        )
        self.assertEqual(witnesses, [])
        self.assertFalse(any(row["selected_by_policy"] for row in routed))
        self.assertEqual(coverage[0]["coverage_status"], "no_rank1_in_handoff")

    def test_missing_view_record_fails(self):
        with self.assertRaises(ValueError):
            select_routes(scenarios(), views()[:-1], ["stan", "a10_informed"])

    def test_identity_mismatch_fails(self):
        records = views()
        records[0]["gene"] = "WRONG"
        with self.assertRaises(ValueError):
            select_routes(scenarios(), records, ["stan"])

    def test_duplicate_view_record_fails(self):
        records = views()
        records.append(records[0].copy())
        with self.assertRaises(ValueError):
            select_routes(scenarios(), records, ["stan"])

    def test_invalid_rank_fails(self):
        records = views()
        records[0]["rank"] = "NaN"
        with self.assertRaises(ValueError):
            select_routes(scenarios(), records, ["stan"])

    def test_disabled_still_preserves_witnesses(self):
        routed, witnesses, _ = select_routes(
            scenarios(), views(), ["stan", "a10_informed"], mode="off"
        )
        self.assertEqual(len(witnesses), 2)
        self.assertTrue(all(row["route_status"] == "disabled" for row in routed))
        self.assertFalse(any(row["selected_for_assessment"] for row in routed))

    def test_witness_score_is_not_converted(self):
        _, witnesses, _ = select_routes(
            scenarios(), views(), ["stan"]
        )
        self.assertEqual(witnesses[0]["source_score"], "0.123")
        self.assertEqual(
            witnesses[0]["source_score_type"],
            "synthetic_view_specific_quantity",
        )

    def test_unscored_record_is_not_selected(self):
        records = [view("A", "stan", 2), view("B", "stan", 1)]
        records[1]["view_status"] = "not_scored"
        records[1]["rank"] = ""
        records[1]["score"] = ""
        _, witnesses, coverage = select_routes(
            scenarios(), records, ["stan"]
        )
        self.assertEqual(witnesses, [])
        self.assertEqual(coverage[0]["not_scored"], 1)

if __name__ == "__main__":
    unittest.main()
