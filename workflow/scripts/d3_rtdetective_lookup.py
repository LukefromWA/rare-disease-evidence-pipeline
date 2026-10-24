#!/usr/bin/env python3
"""Stream pinned four-column RTDetective tracks without reversing coordinates."""
import argparse
import csv
from datetime import datetime, timezone
import hashlib
import json
import math
from pathlib import Path
import tempfile

TRACKS = {"DAP_tga.txt", "Clitocine_tga.txt", "SRI_tga.txt"}
STATUS_FIELDS = ["chrom","target_start","target_end","status","n_exact_records",
                 "n_nonforward_records_near_target"]
EXACT_FIELDS = ["chrom","target_start","target_end","coordinate_order","raw_start",
                "raw_end","author_score_raw","source_line_number"]

def require(condition,message):
    if not condition: raise ValueError(message)

def sha(path): return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def targets_from_mapping(path):
    with Path(path).open(newline="",encoding="utf-8-sig") as stream:
        reader=csv.DictReader(stream,delimiter="\t")
        require(reader.fieldnames and len(set(reader.fieldnames))==len(reader.fieldnames),"Invalid mapping headers")
        rows=list(reader)
    require(all(None not in r and None not in r.values() for r in rows),"Malformed mapping rows")
    targets=set()
    for row in rows:
        if row.get("alternate_codon_RNA") != "UGA": continue
        require(row.get("assembly")=="GRCh38" and row.get("strand")=="1","Unsupported mapping scope")
        for key in ("genomic_coordinate_mapping_verified","genomic_context_matches_frozen_CDS","codon_contiguous"):
            require(row.get(key,"").lower()=="true",f"Unverified mapping flag: {key}")
        require(row.get("track_lookup_status")=="contiguous_stop_codon_lookup_key","Invalid lookup status")
        start,end=int(row["bed_chromStart"]),int(row["bed_chromEnd"])
        require(start>=0 and end-start==3,"Expected three-base BED interval")
        targets.add((row["chrom"],start,end))
    return sorted(targets)

def scan_lines(lines, targets):
    targets=set(targets)
    by_chrom={}
    for key in targets: by_chrom.setdefault(key[0],[]).append(key)
    exact=[]; nonforward={key:0 for key in targets}; counts={key:0 for key in targets}
    data_rows=0; header_seen=False
    for line_number,raw in enumerate(lines,1):
        line=raw.decode("utf-8-sig") if isinstance(raw,bytes) else raw
        words=line.strip().split()
        if not words: continue
        if words[0].startswith("#"):
            if not header_seen:
                require(words==["#chrom","chromStart","chromEnd","score"],"Unexpected track header")
                header_seen=True
            continue
        require(header_seen,"Missing track header")
        require(len(words)==4,f"Malformed track row {line_number}")
        chrom,a,b,score=words
        try: a,b=int(a),int(b); value=float(score)
        except ValueError as exc: raise ValueError(f"Invalid numeric track row {line_number}") from exc
        require(a>=0 and b>=0 and math.isfinite(value),f"Invalid coordinate/score at row {line_number}")
        data_rows+=1
        key=(chrom,a,b)
        if a<b and key in targets:
            counts[key]+=1
            exact.append(dict(chrom=chrom,target_start=a,target_end=b,coordinate_order="ascending",
                raw_start=a,raw_end=b,author_score_raw=score,source_line_number=line_number))
        elif a>=b:
            lo,hi=min(a,b),max(a,b)
            for target in by_chrom.get(chrom,[]):
                _,start,end=target
                if max(lo,start)<min(hi,end) or (lo==hi and start<=lo<end):
                    nonforward[target]+=1
    require(header_seen,"Empty track or missing header")
    statuses=[]
    for chrom,start,end in sorted(targets):
        key=(chrom,start,end); n=counts[key]; backwards=nonforward[key]
        status=("unique_exact_ascending_coordinate_match" if n==1 and backwards==0 else
                "nonforward_records_overlap_target" if backwards else
                "duplicate_exact_ascending_coordinate_matches" if n>1 else
                "no_exact_ascending_coordinate_match")
        statuses.append(dict(chrom=chrom,target_start=start,target_end=end,status=status,
            n_exact_records=n,n_nonforward_records_near_target=backwards))
    return statuses,exact,data_rows

def write_table(path,rows,fields):
    with path.open("w",newline="",encoding="utf-8") as stream:
        writer=csv.DictWriter(stream,fieldnames=fields,delimiter="\t",lineterminator="\n")
        writer.writeheader(); writer.writerows(rows)

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mapping-table",type=Path,required=True)
    parser.add_argument("--mapping-manifest",type=Path,required=True)
    parser.add_argument("--track",type=Path,required=True)
    parser.add_argument("--inventory-dir",type=Path,required=True)
    parser.add_argument("--outdir",type=Path,required=True)
    args=parser.parse_args()
    mapping=args.mapping_table.resolve(); track=args.track.resolve(); inventory=args.inventory_dir.resolve()
    out=args.outdir.resolve()
    require(track.name in TRACKS,"Unsupported track filename")
    require(out!=track.parent and out!=inventory,"Output directory must be separate from resources")
    metadata_path=inventory/"article_metadata.json"; inventory_manifest_path=inventory/"manifest.json"
    mapping_manifest_path=args.mapping_manifest.resolve()
    before={str(p):sha(p) for p in (mapping,mapping_manifest_path,metadata_path,inventory_manifest_path,Path(__file__).resolve())}
    mapping_manifest=json.loads(mapping_manifest_path.read_text())
    require(mapping_manifest.get("output_sha256",{}).get(mapping.name)==before[str(mapping)],"Mapping table is not pinned by its manifest")
    metadata=json.loads(metadata_path.read_text()); inventory_manifest=json.loads(inventory_manifest_path.read_text())
    require(before[str(metadata_path)]==inventory_manifest["metadata_sha256"],"Inventory metadata checksum mismatch")
    require(metadata.get("id")==23708901 and metadata.get("version")==2,"Unexpected article/version")
    entries=[entry for entry in metadata["files"] if entry.get("name")==track.name]
    require(len(entries)==1,"Track not uniquely pinned in inventory")
    entry=entries[0]; expected_md5=(entry.get("computed_md5") or entry.get("supplied_md5") or "").lower()
    require(len(expected_md5)==32,"Missing published track MD5")
    targets=targets_from_mapping(mapping)
    stat_before=track.stat(); md5=hashlib.md5(); sha256=hashlib.sha256(); byte_count=0
    def stream_hashed():
        nonlocal byte_count
        with track.open("rb") as stream:
            for raw in stream:
                md5.update(raw); sha256.update(raw); byte_count+=len(raw)
                yield raw
    print(f"Scanning and hash-verifying {track.name}; {len(targets)} unique UGA target(s)",flush=True)
    statuses,exact,data_rows=scan_lines(stream_hashed(),targets)
    stat_after=track.stat()
    require((stat_before.st_dev,stat_before.st_ino,stat_before.st_size,stat_before.st_mtime_ns,stat_before.st_ctime_ns)==
            (stat_after.st_dev,stat_after.st_ino,stat_after.st_size,stat_after.st_mtime_ns,stat_after.st_ctime_ns),"Track changed during scan")
    require(byte_count==int(entry["size"]) and md5.hexdigest()==expected_md5,"Published track size/MD5 mismatch")
    for path,expected in before.items(): require(sha(path)==expected,f"Input changed during scan: {path}")
    out.mkdir(parents=True,exist_ok=True)
    with tempfile.TemporaryDirectory(prefix=".lookup_",dir=out) as temporary:
        stage=Path(temporary)
        write_table(stage/"lookup_status.tsv",statuses,STATUS_FIELDS)
        write_table(stage/"exact_matches.tsv",exact,EXACT_FIELDS)
        manifest={"node":"rtdetective_direction_preserving_coordinate_lookup",
            "created_utc":datetime.now(timezone.utc).isoformat(),
            "source_article_id":23708901,"source_article_version":2,
            "mapping_table_sha256":before[str(mapping)],"source_track":str(track),
            "source_track_bytes":byte_count,"source_track_md5":md5.hexdigest(),
            "source_track_sha256":sha256.hexdigest(),"published_md5_verified":True,
            "descending_intervals_accepted_as_forward_matches":False,"duplicate_matches_aggregated":False,
            "nonforward_near_target_definition":"normalized half-open overlap with target, or zero-length point inside target; never adopted as forward",
            "target_scope":"verified_positive_strand_GRCh38_UGA_three_base_intervals",
            "target_count":len(targets),"source_data_rows":data_rows,
            "exact_records":len(exact),"raw_score_strings_preserved":True,
            "source_scan_performed":True,"network_requests":False,
            "input_and_implementation_sha256":before,
            "output_sha256":{name:sha(stage/name) for name in ("lookup_status.tsv","exact_matches.tsv")}}
        (stage/"manifest.json").write_text(json.dumps(manifest,indent=2)+"\n")
        for name in ("lookup_status.tsv","exact_matches.tsv","manifest.json"):
            (stage/name).replace(out/name)
    print(json.dumps({"track":track.name,"targets":len(targets),"exact_records":len(exact),
                      "status_counts":{status:sum(r['status']==status for r in statuses) for status in sorted({r['status'] for r in statuses})}}))

if __name__=="__main__": main()
