#!/usr/bin/env python3
"""Regenerate CHECKSUMS.sha256, and keep PROVENANCE.yaml sha256 fields honest.

Two facts about every data file must agree: the digest in CHECKSUMS.sha256 (what
`check_no_leakage` reads) and the digest declared in the dataset's own
PROVENANCE.yaml (what a reader of that record sees). They are written in two
places on purpose — the record is meant to be readable on its own — so something
has to check they have not drifted apart.

    python tools/checksums.py            # rewrite CHECKSUMS.sha256, report drift
    python tools/checksums.py --check    # report only; non-zero exit on drift
"""
from __future__ import annotations
import hashlib, pathlib, sys, yaml

ROOT = pathlib.Path(__file__).resolve().parent.parent
DATA_EXT = {".csv", ".tab", ".txt", ".xls", ".xlsx", ".rda", ".mod"}
MANIFEST = ROOT / "CHECKSUMS.sha256"
HEADER = [
    "# sha256 of every data file in the corpus. Regenerate with tools/checksums.py;",
    "# verify with: sha256sum --check --strict CHECKSUMS.sha256",
]


def sha256(path: pathlib.Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


check_only = "--check" in sys.argv[1:]

data_files = sorted(p for p in ROOT.rglob("*")
                    if p.suffix in DATA_EXT and ".git" not in p.parts
                    and "LICENSES" not in p.parts and p.name != "PROVENANCE.yaml")
digests = {p: sha256(p) for p in data_files}

lines = HEADER + [f"{d}  {p.relative_to(ROOT)}" for p, d in digests.items()]
content = "\n".join(lines) + "\n"

problems: list[str] = []
stale = not MANIFEST.exists() or MANIFEST.read_text() != content
if stale and check_only:
    # Only a failure in --check mode; in write mode rewriting it *is* the fix.
    problems.append("CHECKSUMS.sha256 is out of date (run without --check to rewrite)")
if stale and not check_only:
    MANIFEST.write_text(content)
    print(f"rewrote {MANIFEST.relative_to(ROOT)}")

# Cross-check every digest declared inside a PROVENANCE.yaml against the file.
for prov in sorted(ROOT.rglob("PROVENANCE.yaml")):
    if ".git" in prov.parts:
        continue
    rec = yaml.safe_load(prov.read_text()) or {}
    for entry in rec.get("files", []):
        if not isinstance(entry, dict):
            continue
        target = prov.parent / entry.get("path", "")
        rel = prov.relative_to(ROOT)
        if not target.exists():
            problems.append(f"{rel}: declares missing file {entry.get('path')!r}")
        elif entry.get("sha256") != digests.get(target.resolve()):
            problems.append(f"{rel}: sha256 mismatch for {entry.get('path')!r}")

for p in problems:
    print(f"  x {p}", file=sys.stderr)
if problems:
    sys.exit(1)
print(f"OK: {len(digests)} data files, checksums and provenance records agree.")
