#!/usr/bin/env python3
"""Pre-publication scanner/exporter for the public repository. Reports; never rewrites files.

Scans text files for identifying strings (user paths, hostnames, emails, key assignments) plus your own
terms file (kept OUT of git), lists files it cannot scan (binary/compressed), flags large files, and
optionally copies a clean tree to an export directory. Export is refused if any blocking finding exists.
"""
import argparse, re, shutil, sys
from pathlib import Path

SKIP_DIRS = {".git", "__pycache__", ".snakemake", "native", "raw", "scratch", "cache", "checkpoints", "private"}
PATTERNS = [
    ("windows_user_path", "block", re.compile(r"[A-Za-z]:[\\/]+Users[\\/]+[^\\/\s\"']+")),
    ("wsl_user_path", "block", re.compile(r"/mnt/[a-z]/Users/[^/\s\"']+")),
    ("home_path", "block", re.compile(r"/home/[^/\s\"']+")),
    ("hostname", "block", re.compile(r"\b(?:DESKTOP|LAPTOP)-[A-Z0-9]{5,}\b")),
    ("email", "block", re.compile(r"[A-Za-z0-9._%+-]+@[A-Za-z0-9.-]+\.[A-Za-z]{2,}")),
    ("credential_assignment", "block", re.compile(r"(?i)\b(api[_-]?key|token|secret|password)\b\s*[=:]\s*[^\s\"'<>$]{8,}")),
    ("private_dir_reference", "review", re.compile(r"results[\\/]private|handoffs[\\/]private")),
]
MAX_MIB = 50

def mask(s): return s[:3] + "***" if len(s) > 6 else "***"

def iter_files(paths):
    for base in paths:
        base = Path(base)
        if base.is_file():
            yield base; continue
        for p in sorted(base.rglob("*")):
            if p.is_file() and not (set(p.relative_to(base).parts[:-1]) & SKIP_DIRS):
                yield p

def scan(paths, terms):
    findings, unscanned, large = [], [], []
    scanned = []
    for p in iter_files(paths):
        if p.stat().st_size > MAX_MIB * 1048576: large.append(p)
        raw = p.read_bytes()
        if b"\x00" in raw[:8192] or p.suffix in {".gz", ".rds", ".zip", ".sqlite", ".png", ".jpg", ".pdf"}:
            unscanned.append(p); continue
        scanned.append(p)
        for n, line in enumerate(raw.decode("utf-8", "replace").splitlines(), 1):
            for name, sev, rx in PATTERNS:
                for m in rx.finditer(line): findings.append((p, n, name, sev, mask(m.group(0))))
            low = line.lower()
            for t in terms:
                if t.lower() in low: findings.append((p, n, "user_term", "block", mask(t)))
    return findings, unscanned, large, scanned

def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--paths", nargs="+", required=True, help="Files/dirs to scan (and export)")
    ap.add_argument("--terms", help="File of identifying terms, one per line (do NOT commit it)")
    ap.add_argument("--export-dir", help="Copy scanned files here only if no blocking findings")
    a = ap.parse_args(argv)
    terms = []
    if a.terms:
        terms = [x.strip() for x in Path(a.terms).read_text(encoding="utf-8").splitlines() if x.strip() and not x.startswith("#")]
    findings, unscanned, large, scanned = scan(a.paths, terms)
    blocking = [f for f in findings if f[3] == "block"]
    for p, n, name, sev, m in findings: print("%s  %s:%d  %s  [%s]" % (sev.upper(), p, n, name, m))
    for p in unscanned: print("UNSCANNED  %s (binary/compressed: inspect or exclude)" % p)
    for p in large: print("LARGE  %s (> %d MiB: keep out of Git)" % (p, MAX_MIB))
    print("\nscanned=%d unscanned=%d blocking=%d review=%d large=%d" % (len(scanned), len(unscanned), len(blocking), len(findings) - len(blocking), len(large)))
    if a.export_dir:
        if blocking or large: print("EXPORT REFUSED: resolve blocking/large findings first."); return 1
        dest = Path(a.export_dir)
        if dest.exists() and any(dest.iterdir()): print("EXPORT REFUSED: destination not empty."); return 1
        keep = set(scanned)
        for base in map(Path, a.paths):
            for p in ([base] if base.is_file() else list(iter_files([base]))):
                if p not in keep: continue
                rel = Path(p.name) if base.is_file() else Path(base.name) / p.relative_to(base)
                (dest / rel).parent.mkdir(parents=True, exist_ok=True); shutil.copy2(p, dest / rel)
        print("EXPORTED scanned text files to", dest, "- unscanned binary/compressed files were NOT copied; review them yourself.")
    return 1 if blocking else 0

if __name__ == "__main__":
    sys.exit(main())
