#!/usr/bin/env python3
"""Every dataset directory has a complete PROVENANCE.yaml; CI gate.

Unresolved licenses are reported but do not fail the build — they are a known,
tracked state (see LICENSE.md). A *missing* or *incomplete* record does fail:
the failure mode this guards against is a dataset quietly appearing with no
stated origin, which is how aggregation repos accumulate legal debt.
"""
from __future__ import annotations
import sys, pathlib, yaml

ROOT = pathlib.Path(__file__).resolve().parent.parent
REQUIRED = {"id", "title", "license", "upstream", "content", "files"}
LICENSE_REQUIRED = {"spdx", "holder", "verified", "redistribution"}
DATA_EXT = {".csv", ".tab", ".txt", ".xls", ".xlsx", ".rda", ".mod"}
# .txt is both a data extension here (Monolix outputs, warfarin) and the
# extension REUSE requires for license texts, so LICENSES/ has to be excluded by
# path rather than by suffix. Nothing in it is data and nothing in it needs a
# provenance record.
SKIP_DIRS = {".git", "LICENSES"}

errors: list[str] = []
unresolved: list[str] = []

def dataset_root(f: pathlib.Path) -> pathlib.Path | None:
    """Nearest ancestor holding a PROVENANCE.yaml.

    A dataset is a directory tree, not a single directory: ggPMX's theophylline
    set spans ggpmx/, ggpmx/Monolix/ and ggpmx/Monolix/ChartsData/Saem/ and is
    one dataset under one license. Requiring a record per subdirectory would be
    noise, so a record covers everything beneath it.
    """
    for d in f.parents:
        if d == ROOT.parent:
            break
        if (d / "PROVENANCE.yaml").exists():
            return d
        if d == ROOT:
            break
    return None


data_files = [p for p in ROOT.rglob("*")
              if p.suffix in DATA_EXT and not SKIP_DIRS & set(p.parts)]

uncovered: dict[pathlib.Path, list[str]] = {}
dirs: set[pathlib.Path] = set()
for f in data_files:
    if (root := dataset_root(f)) is not None:
        dirs.add(root)
    else:
        uncovered.setdefault(f.parent.relative_to(ROOT), []).append(f.name)

for d, names in sorted(uncovered.items()):
    errors.append(f"{d}: no PROVENANCE.yaml covers {sorted(names)[:3]}"
                  + (" ..." if len(names) > 3 else ""))

for d in sorted(dirs):
    rel = d.relative_to(ROOT)
    prov = d / "PROVENANCE.yaml"
    if not prov.exists():
        errors.append(f"{rel}: missing PROVENANCE.yaml")
        continue
    try:
        rec = yaml.safe_load(prov.read_text()) or {}
    except yaml.YAMLError as exc:
        errors.append(f"{rel}: unparseable PROVENANCE.yaml ({exc})")
        continue

    if missing := REQUIRED - rec.keys():
        errors.append(f"{rel}: missing fields {sorted(missing)}")
    lic = rec.get("license") or {}
    if missing := LICENSE_REQUIRED - lic.keys():
        errors.append(f"{rel}: license missing {sorted(missing)}")

    spdx = str(lic.get("spdx", ""))
    if spdx.startswith("LicenseRef-Unresolved") or lic.get("redistribution") == "unresolved":
        unresolved.append(f"{rel}: {spdx}")

    declared = {f["path"] for f in rec.get("files", []) if isinstance(f, dict)}
    actual = {str(p.relative_to(d)) for p in d.rglob("*")
              if p.suffix in DATA_EXT and dataset_root(p) == d}
    if undeclared := actual - declared:
        errors.append(f"{rel}: files present but not declared: {sorted(undeclared)}")

if unresolved:
    print("Unresolved licenses (tracked, see LICENSE.md):")
    for u in unresolved:
        print(f"  ~ {u}")

if errors:
    print("\nProvenance errors:", file=sys.stderr)
    for e in errors:
        print(f"  x {e}", file=sys.stderr)
    sys.exit(1)

print(f"\nOK: {len(dirs)} dataset directories, all with complete provenance.")
