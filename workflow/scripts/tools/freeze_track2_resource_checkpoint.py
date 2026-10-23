#!/usr/bin/env python3
"""Freeze a local, resource-scoped checkpoint. Never publish or create readiness flags."""
import argparse, csv, hashlib, json, os, shutil, sys
from pathlib import Path
from datetime import datetime, timezone

def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(1024*1024),b''): h.update(b)
    return h.hexdigest()

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--root',default='.')
    ap.add_argument('--workflow',default='written_final_scripts/workflow')
    ap.add_argument('--out',required=True)
    ap.add_argument('--include-public-data',action='store_true')
    a=ap.parse_args()
    root=Path(a.root).resolve(); wf=(root/a.workflow).resolve()
    dest=Path(a.out)
    if not dest.is_absolute(): dest=root/dest
    if dest.exists(): raise SystemExit('Refusing to overwrite checkpoint: '+str(dest))
    if not wf.is_dir(): raise SystemExit('Workflow directory missing')
    resource=(root/'resources/track2').resolve()
    candidates={}
    def add(p,role):
        if not p.is_file() or p.is_symlink(): return
        rp=p.resolve()
        try: rel=rp.relative_to(root)
        except ValueError: raise SystemExit('Refusing source outside project: '+str(p))
        if any(t in str(rel).lower() for t in ('results/private','data/private','.renviron','patient_hpo','config.private')):
            raise SystemExit('Refusing private-path candidate: '+str(rel))
        candidates[str(rel)]=(rp,role)
    for p in (wf/'scripts').glob('e*.R'): add(p,'resource_node_code')
    # Source libraries may implement other layers, but only code is copied. Review before sharing.
    for pat in ('*.R','*.r'):
        for p in (wf/'lib').glob(pat): add(p,'library_code_review_before_sharing')
    for name in ('Snakefile.track2_resources','config.track2_resources.yaml'):
        add(wf/name,'resource_workflow_config')
    for d in ('docs','contracts'):
        for p in (wf/d).glob('*track2*'):
            if p.suffix in ('.md','.txt','.json'): add(p,'handoff_documentation')
    for parent in ('native','raw','audits','derived','graph','manifests'):
        base=resource/parent
        if base.exists():
            for p in sorted(base.rglob('*')):
                if p.is_file() and not p.is_symlink():
                    if p.name.endswith('.partial'): continue
                    add(p,'public_resource_'+parent)
    scratch=resource/'scratch'
    if scratch.exists():
        for d in sorted(scratch.iterdir()):
            if d.is_dir() and d.name.startswith(('e10_','e11_')):
                for p in sorted(d.rglob('*')): add(p,'resource_scratch_diagnostic')
    if not candidates: raise SystemExit('No eligible files found')
    dest.mkdir(parents=True)
    (dest/'INCOMPLETE').write_text('Checkpoint is incomplete until exporter succeeds.\n')
    rows=[]
    for rel,(src,role) in sorted(candidates.items()):
        # Code/config/docs and metadata copied by default; large reference data opt in.
        metadata=src.name.endswith(('.json','.sha256')) or role in ('resource_node_code','library_code_review_before_sharing','resource_workflow_config','handoff_documentation')
        copy_data=a.include_public_data or metadata
        before=src.stat(); before_sha=sha(src)
        row=dict(relative_path=rel,role=role,bytes=before.st_size,sha256=before_sha,
                 original_mtime_ns=before.st_mtime_ns,copied=copy_data)
        if copy_data:
            dst=dest/'files'/rel; dst.parent.mkdir(parents=True,exist_ok=True)
            shutil.copyfile(src,dst)
            if sha(dst)!=before_sha or sha(src)!=before_sha: raise RuntimeError('File changed while freezing: '+rel)
        rows.append(row)
    sentinel_names=('.ready','.signed_propagation_enabled','.signedpropagationenabled')
    report=dict(checkpoint_version='resource_checkpoint_v1',created_utc=datetime.now(timezone.utc).isoformat(),
                scope='local resource code/metadata; not a publication clearance or readiness validation',
                include_public_data=a.include_public_data,artifacts=rows,
                observed_sentinels={n:(resource/n).exists() for n in sentinel_names},
                sentinels_copied=False,private_case_data_copied=False,
                warning='Review code/config/manifests for embedded sensitive values before sharing. Path exclusions are not a content/privacy audit. Metadata-only mode requires the exact external bytes to restore the checkpoint.')
    (dest/'checkpoint.json').write_text(json.dumps(report,indent=2))
    copied=[r for r in rows if r['copied']]
    with (dest/'checksums.sha256').open('w') as f:
        for r in copied: f.write(r['sha256']+'  files/'+r['relative_path']+'\n')
    with (dest/'inventory.tsv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
    (dest/'README.txt').write_text('Local reproducibility checkpoint, not a ready resource bundle.\nVerify copied bytes from this directory: sha256sum -c checksums.sha256\nNo readiness or propagation flags are copied. Restore files deliberately into a new location and revalidate manifests; existing absolute paths may need audited rebasing.\nNative/prepared data included only with --include-public-data. Observe each source license.\nNo case outputs are collected. Review copied code/config and metadata before publishing; do not treat automatic path exclusions as privacy clearance.\n')
    (dest/'INCOMPLETE').unlink()
    print('CHECKPOINT COMPLETE')
    print('Inventoried files:',len(rows),'| copied:',len(copied),'| copied bytes:',sum(r['bytes'] for r in copied))
    print('Destination:',dest)
    print('No source edits, downloads, model fitting, ready flags, or publication performed.')
if __name__=='__main__': main()
