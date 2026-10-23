import unittest
from d3_readthrough_context_bridge import check_record

def fixture():
    request = {
        "hypothesis_id": "synthetic_hypothesis",
        "gene": "SYNTHETIC_GENE",
        "variant": "1:101:A:G",
        "q_j_raw": "0.2",
        "selecting_views": "synthetic_view",
        "selected_for_assessment": "True",
        "route_status": "eligible_for_context_verification",
    }
    context = {
        "gene": "SYNTHETIC_GENE",
        "candidate_id": "1:101:A:G",
        "transcript_version": "SYNTHETIC_TRANSCRIPT.1",
        "HGVSc": "synthetic_coding_change",
        "HGVSp": "synthetic_stop_change",
        "alternate_codon_RNA": "UGA",
        "alternate_context_bracketed_RNA": "AAA[UGA]CCC",
        "sequence_verified": "True",
        "alternate_is_stop": "True",
        "context_boundary_truncated": "False",
        "reference_assembly": "GRCh38",
        "phase_status": "unknown",
        "haplotype_scope": "synthetic",
        "NMD_status": "unknown",
    }
    mapping = {
        field: context[field] for field in [
            "gene", "candidate_id", "transcript_version",
            "HGVSc", "HGVSp", "alternate_codon_RNA",
            "alternate_context_bracketed_RNA",
        ]
    }
    mapping.update({
        "assembly": "GRCh38", "chrom": "chr1", "strand": "1",
        "bed_chromStart": "99", "bed_chromEnd": "102",
        "genomic_coordinate_mapping_verified": "True",
        "genomic_context_matches_frozen_CDS": "True",
        "codon_contiguous": "True",
        "track_lookup_status": "contiguous_stop_codon_lookup_key",
    })
    return request, context, mapping

class BridgeTests(unittest.TestCase):
    def test_valid_cached_binding(self):
        request, context, mapping = fixture()
        result = check_record(request, [context], [mapping])
        self.assertEqual(
            result["bridge_status"], "verified_cached_context_and_mapping"
        )

    def test_missing_context(self):
        request, _, mapping = fixture()
        self.assertEqual(
            check_record(request, [], [mapping])["bridge_status"],
            "missing_context",
        )

    def test_duplicate_context(self):
        request, context, mapping = fixture()
        self.assertEqual(
            check_record(request, [context, context], [mapping])["bridge_status"],
            "ambiguous_context",
        )

    def test_missing_mapping(self):
        request, context, _ = fixture()
        self.assertEqual(
            check_record(request, [context], [])["bridge_status"],
            "missing_mapping",
        )

    def test_transcript_mismatch(self):
        request, context, mapping = fixture()
        mapping["transcript_version"] = "OTHER.2"
        self.assertEqual(
            check_record(request, [context], [mapping])["bridge_status"],
            "identity_or_transcript_mismatch",
        )

    def test_unverified_sequence(self):
        request, context, mapping = fixture()
        context["sequence_verified"] = "False"
        self.assertEqual(
            check_record(request, [context], [mapping])["bridge_status"],
            "verification_failed",
        )

    def test_negative_strand_not_reinterpreted(self):
        request, context, mapping = fixture()
        mapping["strand"] = "-1"
        self.assertEqual(
            check_record(request, [context], [mapping])["bridge_status"],
            "unsupported_resource_scope",
        )

    def test_variant_outside_codon(self):
        request, context, mapping = fixture()
        mapping["bed_chromStart"] = "199"
        mapping["bed_chromEnd"] = "202"
        self.assertEqual(
            check_record(request, [context], [mapping])["bridge_status"],
            "identity_or_transcript_mismatch",
        )

    def test_context_boundary(self):
        request, context, mapping = fixture()
        context["context_boundary_truncated"] = "True"
        self.assertEqual(
            check_record(request, [context], [mapping])["bridge_status"],
            "verification_failed",
        )

    def test_nonapplicable_stays_nonapplicable(self):
        request, context, mapping = fixture()
        request["route_status"] = "not_applicable"
        self.assertEqual(
            check_record(request, [context], [mapping])["bridge_status"],
            "not_applicable_upstream",
        )

if __name__ == "__main__":
    unittest.main()
