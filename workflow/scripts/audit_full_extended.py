#!/usr/bin/env python3
"""Verify extended output bundles and record completion without interpreting scores."""
import argparse
from collections import Counter
import csv
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import tempfile

def sha(path): return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def require(condition,message):
    if not condition: raise ValueError(message)
def table(path):
    with path.open(newline="",encoding="utf-8-sig") as stream:
        reader=csv.DictReader(stream,delimiter="\t")
        require(reader.fieldnames and len(set(reader.fieldnames))==len(reader.fieldnames),f"Invalid headers: {path}")
        rows=list(reader)
    require(all(None not in row and None not in row.values() for row in rows),f"Malformed table: {path}")
    return rows

def inspect(case):
    bundles={}; hashes={}
    for relative in ("d3/codon_context","d3/codon_context/genomic_mapping","d3/readthrough_route",
        "d3/readthrough_route_multiview","d3/readthrough_context_bridge","d3/rtdetective_adoption",
        "d3/claim_ledger","d3/claim_ledger_mechanistic"):
        directory=case/relative; path=directory/"manifest.json"
        manifest=json.loads(path.read_text()); hashes[str(path)]=sha(path)
        outputs=manifest.get("output_sha256")
        require(isinstance(outputs,dict) and bool(outputs),f"Missing output hashes: {path}")
        for name,expected in outputs.items():
            rel=Path(name)
            require(not rel.is_absolute() and ".." not in rel.parts,"Unsafe manifest path")
            output=directory/rel
            require(sha(output)==expected,f"Output hash mismatch: {output}")
            hashes[str(output)]=expected
        bundles[relative]=manifest
    routing=table(case/"d3/readthrough_route_multiview/routing.tsv")
    requests={(r["hypothesis_id"],r["gene"],r["variant"]) for r in routing
        if r.get("selected_for_assessment","").lower()=="true" and
        r.get("route_status")=="eligible_for_context_verification"}
    bridge=table(case/"d3/readthrough_context_bridge/verification.tsv")
    verified={(r["hypothesis_id"],r["gene"],r["variant"]) for r in bridge
        if r.get("bridge_status")=="verified_cached_context_and_mapping"}
    require(requests==verified,"Eligible context requests do not exactly match verified bridge identities")
    predictions=table(case/"d3/rtdetective_adoption/predictions.tsv")
    statuses=table(case/"d3/rtdetective_adoption/adoption_status.tsv")
    require(all((r["hypothesis_id"],r["gene"],r["variant"]) in verified for r in predictions),"Prediction without verified bridge identity")
    keys=[(r["hypothesis_id"],r["track_name"]) for r in statuses]
    require(len(keys)==len(set(keys)),"Duplicate adoption statuses")
    expected_status_keys={(identity[0],track) for identity in verified for track in ("DAP_tga.txt","Clitocine_tga.txt","SRI_tga.txt")}
    require(set(keys)==expected_status_keys,"Incomplete or unexpected adoption track coverage")
    adopted={(r["hypothesis_id"],r["track_name"]) for r in statuses if r["adoption_status"]=="adopted_verified_cached_lookup"}
    prediction_keys=[(r["hypothesis_id"],r["track_name"]) for r in predictions]
    require(len(prediction_keys)==len(set(prediction_keys)) and set(prediction_keys)==adopted,"Adoption/prediction identity disagreement")
    lookup_statuses={}
    for name in ("DAP_tga","Clitocine_tga","SRI_tga"):
        directory=case/"d3/rtdetective_lookup"/name
        manifest_path=directory/"manifest.json"; manifest=json.loads(manifest_path.read_text())
        require(manifest.get("published_md5_verified") is True,"Lookup lacks published MD5 verification")
        require(manifest.get("descending_intervals_accepted_as_forward_matches") is False,"Unsafe lookup direction policy")
        require(manifest.get("duplicate_matches_aggregated") is False,"Unsafe duplicate policy")
        require(manifest.get("mapping_table_sha256")==sha(case/"d3/codon_context/genomic_mapping/codon_genomic_mapping.tsv"),"Lookup mapping fingerprint mismatch")
        hashes[str(manifest_path)]=sha(manifest_path)
        for filename,expected in manifest["output_sha256"].items():
            rel=Path(filename); require(not rel.is_absolute() and ".." not in rel.parts,"Unsafe lookup output path")
            require(sha(directory/rel)==expected,"Lookup output hash mismatch")
            hashes[str(directory/rel)]=expected
        lookup_statuses[name]=dict(Counter(r["status"] for r in table(directory/"lookup_status.tsv")))
    for relative in ("a10/ranked_candidates.tsv.gz","c4_stan_union/current_Stan_final_hypothesis_set.tsv","d7/final_report.tsv"):
        path=case/relative; require(path.is_file() and path.stat().st_size>0,f"Missing core output: {path}")
        hashes[str(path)]=sha(path)
    ledger=bundles["d3/claim_ledger_mechanistic"]
    require(ledger.get("automatic_scientific_adjudication") is False,"Unexpected automatic adjudication")
    claims=table(case/"d3/claim_ledger_mechanistic/claims.tsv")
    queue=table(case/"d3/claim_ledger_mechanistic/review_queue.tsv")
    require({r["claim_id"] for r in claims}=={r["claim_id"] for r in queue},"Incomplete claim review queue")
    return {"node":"fresh_full_extended_completion_audit_v1",
        "created_utc":datetime.now(timezone.utc).isoformat(),"case_directory":str(case),
        "eligible_context_requests":len(requests),"verified_context_requests":len(verified),
        "adoption_status_counts":dict(Counter(r["adoption_status"] for r in statuses)),
        "prediction_rows":len(predictions),"lookup_status_counts":lookup_statuses,
        "mechanistic_claim_rows":len(claims),"claim_review_rows":len(queue),
        "output_and_core_sha256":hashes,
        "completion_means":"requested output bundles exist and implemented integrity/identity checks pass",
        "all_track_predictions_required_for_completion":False,
        "core_D7_report_recomputed_with_extended_evidence":False,
        "prediction_is_independent_drug_discovery":False,"clinical_suitability_assessed":False,
        "scientific_review_complete":False,"network_requests":False}

def main():
    p=argparse.ArgumentParser(description=__doc__); p.add_argument("--case-dir",type=Path,required=True); p.add_argument("--out",type=Path,required=True)
    args=p.parse_args(); result=inspect(args.case_dir.resolve())
    result["implementation_sha256"]=sha(Path(__file__))
    args.out.parent.mkdir(parents=True,exist_ok=True)
    with tempfile.NamedTemporaryFile(mode="w",dir=args.out.parent,delete=False) as stream:
        json.dump(result,stream,indent=2); stream.write("\n"); temporary=Path(stream.name)
    temporary.replace(args.out)
    print(json.dumps({key:result[key] for key in ("eligible_context_requests","prediction_rows","mechanistic_claim_rows","adoption_status_counts")}))
if __name__=="__main__": main()
