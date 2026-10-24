import csv
import hashlib
import json
from pathlib import Path
import tempfile
import unittest
from audit_full_extended import inspect

def sha(path): return hashlib.sha256(path.read_bytes()).hexdigest()
def write_table(path,rows):
    path.parent.mkdir(parents=True,exist_ok=True)
    fields=list(rows[0])
    with path.open("w",newline="") as stream:
        w=csv.DictWriter(stream,fieldnames=fields,delimiter="\t"); w.writeheader(); w.writerows(rows)
def repin(directory,extra=None):
    data={"output_sha256":{p.name:sha(p) for p in directory.iterdir() if p.is_file() and p.name!="manifest.json"}}
    if extra: data.update(extra)
    (directory/"manifest.json").write_text(json.dumps(data))
def fixture(case):
    for relative in ("d3/codon_context","d3/codon_context/genomic_mapping","d3/readthrough_route","d3/readthrough_route_multiview","d3/readthrough_context_bridge","d3/rtdetective_adoption","d3/claim_ledger","d3/claim_ledger_mechanistic"):
        directory=case/relative; directory.mkdir(parents=True,exist_ok=True); (directory/"marker.txt").write_text("synthetic\n")
    identity=dict(hypothesis_id="SYNTHETIC",gene="SYNTHETIC",variant="1:101:A:T")
    write_table(case/"d3/readthrough_route_multiview/routing.tsv",[{**identity,"selected_for_assessment":"True","route_status":"eligible_for_context_verification"}])
    write_table(case/"d3/readthrough_context_bridge/verification.tsv",[{**identity,"bridge_status":"verified_cached_context_and_mapping"}])
    tracks=["DAP_tga.txt","Clitocine_tga.txt","SRI_tga.txt"]
    write_table(case/"d3/rtdetective_adoption/predictions.tsv",[{**identity,"track_name":t} for t in tracks])
    write_table(case/"d3/rtdetective_adoption/adoption_status.tsv",[{**identity,"track_name":t,"adoption_status":"adopted_verified_cached_lookup"} for t in tracks])
    write_table(case/"d3/claim_ledger_mechanistic/claims.tsv",[{"claim_id":"C1"}])
    write_table(case/"d3/claim_ledger_mechanistic/review_queue.tsv",[{"claim_id":"C1","status":"pending"}])
    mapping=case/"d3/codon_context/genomic_mapping/codon_genomic_mapping.tsv"; mapping.write_text("synthetic_mapping\n")
    for t in tracks:
        directory=case/"d3/rtdetective_lookup"/t.removesuffix(".txt"); directory.mkdir(parents=True)
        write_table(directory/"lookup_status.tsv",[{"status":"unique_exact_ascending_coordinate_match"}])
        repin(directory,{"published_md5_verified":True,"descending_intervals_accepted_as_forward_matches":False,"duplicate_matches_aggregated":False,"mapping_table_sha256":sha(mapping)})
    for relative in ("d3/codon_context","d3/codon_context/genomic_mapping","d3/readthrough_route","d3/readthrough_route_multiview","d3/readthrough_context_bridge","d3/rtdetective_adoption","d3/claim_ledger","d3/claim_ledger_mechanistic"):
        repin(case/relative,{"automatic_scientific_adjudication":False})
    for relative in ("a10/ranked_candidates.tsv.gz","c4_stan_union/current_Stan_final_hypothesis_set.tsv","d7/final_report.tsv"):
        p=case/relative; p.parent.mkdir(parents=True,exist_ok=True); p.write_text("synthetic nonempty output\n")

class AuditTests(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory(); self.case=Path(self.tmp.name); fixture(self.case)
    def tearDown(self): self.tmp.cleanup()
    def test_valid_completion_retains_limitations(self):
        result=inspect(self.case)
        self.assertEqual(result["prediction_rows"],3)
        self.assertFalse(result["scientific_review_complete"])
        self.assertFalse(result["core_D7_report_recomputed_with_extended_evidence"])
    def test_tampered_output(self):
        (self.case/"d3/codon_context/marker.txt").write_text("changed")
        with self.assertRaisesRegex(ValueError,"hash mismatch"): inspect(self.case)
    def test_unverified_request(self):
        directory=self.case/"d3/readthrough_context_bridge"
        write_table(directory/"verification.tsv",[{"hypothesis_id":"SYNTHETIC","gene":"SYNTHETIC","variant":"1:101:A:T","bridge_status":"missing_mapping"}]); repin(directory)
        with self.assertRaisesRegex(ValueError,"bridge identities"): inspect(self.case)
    def test_missing_track_status(self):
        directory=self.case/"d3/rtdetective_adoption"
        write_table(directory/"adoption_status.tsv",[{"hypothesis_id":"SYNTHETIC","gene":"SYNTHETIC","variant":"1:101:A:T","track_name":"DAP_tga.txt","adoption_status":"adopted_verified_cached_lookup"}]); repin(directory)
        with self.assertRaisesRegex(ValueError,"track coverage"): inspect(self.case)
    def test_unsafe_direction_policy(self):
        path=self.case/"d3/rtdetective_lookup/DAP_tga/manifest.json"
        manifest=json.loads(path.read_text()); manifest["descending_intervals_accepted_as_forward_matches"]=True; path.write_text(json.dumps(manifest))
        with self.assertRaisesRegex(ValueError,"direction policy"): inspect(self.case)
    def test_missing_core_output(self):
        (self.case/"d7/final_report.tsv").unlink()
        with self.assertRaisesRegex(ValueError,"Missing core"): inspect(self.case)
if __name__=="__main__": unittest.main()
