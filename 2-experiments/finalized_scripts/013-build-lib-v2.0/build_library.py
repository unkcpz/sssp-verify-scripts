#!/usr/bin/env python3
"""Build the flat web-serving SSSP library layout.

Output layout (version = v2.0 by default):

  sssp/<version>/
    eff_pbe_recommended.tar.gz
    eff_pbesol_recommended.tar.gz
    prec_pbe_recommended.tar.gz
    prec_pbesol_recommended.tar.gz
    metadata/<element>.json
    library/<element>/
      eff_pbe_recommended/<element>...upf
      eff_pbesol_recommended/<element>...upf
      prec_pbe_recommended/<element>...upf
      prec_pbesol_recommended/<element>...upf
      all_upf_files/<element>...upf

Recommended upf files are the finalized sets in
  010-extract-eff-lib/mix-sssp-{pbe,pbesol}-eff-lib-v2/library and
  011-extract-prec-lib/mix-sssp-{pbe,pbesol}-prec-lib-v2/library.

all_upf_files/<element>/ holds every distinct upf content for that element,
collected from libraries-pbe, libraries-pbesol and common/pseudos, deduplicated
by md5.  When distinct contents share a basename the canonical one keeps the
original name and the colliding copies are renamed with a source-folder suffix
(preferred order: not *dont-use*, not *original_upfs*, gbrv upf2 over upf1).
"""

import argparse
import csv
import hashlib
import json
import os
import re
import shutil
import sys
import tarfile
from collections import defaultdict
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[3]

LIB_ROOTS = ["libraries-pbe", "libraries-pbesol"]
COMMON_ROOT = "common/pseudos"

FIN_DIR = REPO_ROOT / "2-experiments" / "finalized_scripts"
RECOMMENDED_DIRS = {
    ("eff", "pbe"): FIN_DIR / "010-extract-eff-lib" / "mix-sssp-pbe-eff-lib-v2" / "library",
    ("eff", "pbesol"): FIN_DIR / "010-extract-eff-lib" / "mix-sssp-pbesol-eff-lib-v2" / "library",
    ("prec", "pbe"): FIN_DIR / "011-extract-prec-lib" / "mix-sssp-prec-pbe-lib-v2" / "library",
    ("prec", "pbesol"): FIN_DIR / "011-extract-prec-lib" / "mix-sssp-prec-pbesol-lib-v2" / "library",
}
CUTOFFS_JSON = {
    "eff": FIN_DIR / "010-extract-eff-lib" / "cutoffs.json",
    "prec": FIN_DIR / "011-extract-prec-lib" / "cutoffs.json",
}

CSV_A = REPO_ROOT / "filtered_pseudopotential_data.csv"
CSV_B = REPO_ROOT / "filtered_pseudopotential_data_N.csv"

# source library folder -> csv `name` value (for the metadata csv cross-ref)
LIB_TO_CSV_NAME = {
    "nc-dojo-v0.4.1-std": "DOJO-041-std",
    "nc-dojo-v0.4.1-str": "DOJO-041-str",
    "nc-dojo-v0.5.0-std": "DOJO-050-std",
    "nc-sg15-oncvpsp4": "SG15",
    "nc-spms-oncvpsp4": "SPMS",
    "paw-jth-v1.1-std": "JTH-1.1-std",
    "paw-jth-v1.1-str": "JTH-1.1-str",
    "PAW-JTH1.1-standard": "JTH-1.1-std",
    "PAW-JTH1.1-stringent": "JTH-1.1-str",
    "paw-psl-v0.x": "PSL-PAW-v0.x",
    "paw-psl-v1.0.0-high": "PSL-PAW-v1-high",
    "paw-psl-v1.0.0-low": "PSL-PAW-v1-low",
    "us-psl-v0.x": "PSL-US-v0.x",
    "us-psl-v1.0.0-high": "PSL-US-v1-high",
    "us-psl-v1.0.0-low": "PSL-US-v1-low",
    "us-gbrv-v1.x-upf2": "GBRV-1.X",
    "us-gbrv-v1.0-upf2": "GBRV-1.X",
    "us-gbrv-v1.5-upf2": "GBRV-1.X",
    "us-gbrv-v1.4-upf2": "GBRV-1.X",
    "us-gbrv-v1.2-upf2": "GBRV-1.X",
    "us-gbrv-v1.x-upf1": "GBRV-1.X",
    "us-gbrv-v1.0-upf1": "GBRV-1.X",
    "us-gbrv-v1.5-upf1": "GBRV-1.X",
    "us-gbrv-v1.4-upf1": "GBRV-1.X",
    "us-gbrv-v1.2-upf1": "GBRV-1.X",
    "us-gbrv-v1.x-upf1-dont-use-this-lib-use-upf2": "GBRV-1.X",
}

ELEM_RE = re.compile(r"^([A-Z][a-z]?)\b")
Z_RE = re.compile(r"\.z_(\d+)\.")
FUNC_RE = re.compile(r"\.(pbe|pbesol)\.")


def md5_bytes(b):
    return hashlib.md5(b).hexdigest()


def sha256_bytes(b):
    return hashlib.sha256(b).hexdigest()


def extract_element(name):
    m = ELEM_RE.match(name)
    return m.group(1) if m else None


def infer_functional(path):
    m = FUNC_RE.search(path.name)
    if m:
        return m.group(1)
    rel = str(path.relative_to(REPO_ROOT))
    if rel.startswith("libraries-pbesol"):
        return "pbesol"
    return "pbe"


def infer_pp_type(name):
    for t in (".nc.", ".us.", ".paw."):
        if t in name:
            return t.strip(".")
    return None


def infer_z_valence(name):
    m = Z_RE.search(name)
    return int(m.group(1)) if m else None


def lib_folder(path, root):
    rel = path.relative_to(root)
    parts = rel.parts
    lib = parts[0] if len(parts) > 1 else str(rel)
    return lib


def pref_key(path):
    parts = [p for p in path.parts]
    key = 0
    if any("dont-use" in p for p in parts):
        key += 1000
    if any(p == "original_upfs" for p in parts):
        key += 500
    if any(re.search(r"upf1($|\D)", p) and "upf2" not in p for p in parts):
        key += 100
    return key


def unique_name(name, used):
    if name not in used:
        return name
    base, ext = os.path.splitext(name)
    i = 1
    while f"{base}_{i}{ext}" in used:
        i += 1
    return f"{base}_{i}{ext}"


def load_csv(path):
    rows = {}
    if not path.exists():
        return rows
    with open(path, newline="") as f:
        for r in csv.DictReader(f):
            rows[(r["element"], r["name"])] = r
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--version", default="v2.0")
    ap.add_argument("--out", default=REPO_ROOT / "sssp")
    args = ap.parse_args()

    out = Path(args.out) / args.version
    lib_out = out / "library"
    meta_out = out / "metadata"
    for d in (lib_out, meta_out):
        d.mkdir(parents=True, exist_ok=True)

    cutoffs = {}
    for calc, p in CUTOFFS_JSON.items():
        cutoffs[calc] = json.loads(p.read_text())
        shutil.copyfile(p, lib_out / f"{calc}_cutoffs.json")
    csv_rows = {}
    csv_rows.update(load_csv(CSV_A))
    csv_rows.update(load_csv(CSV_B))
    csv_rows = {k: v for k, v in csv_rows.items() if k[1]}  # drop empty names

    recommend = {}  # (calc, func, element) -> Path in mix library
    for (calc, func), d in RECOMMENDED_DIRS.items():
        for f in d.glob("*.upf"):
            elem = extract_element(f.name)
            if elem:
                recommend[(calc, func, elem)] = f

    roots = [REPO_ROOT / r for r in LIB_ROOTS]
    roots.append(REPO_ROOT / COMMON_ROOT)

    collected = defaultdict(list)  # element -> [paths]
    for root in roots:
        if not root.exists():
            continue
        for upf in root.rglob("*.upf"):
            elem = extract_element(upf.name)
            if elem:
                collected[elem].append(upf)

    elements = sorted(set(collected) | {e for (_, _, e) in recommend})
    print(f"{len(elements)} elements")

    summary = {"elements": 0, "distinct_contents": 0, "renamed": 0, "upf_files": 0}
    tar_handles = {}
    for func in ("pbe", "pbesol"):
        for calc in ("eff", "prec"):
            gz = str(lib_out / f"{calc}_{func}_recommended.tar.gz")
            tar_handles[(calc, func)] = tarfile.open(gz, "w:gz")

    try:
        for elem in elements:
            paths = sorted(collected.get(elem, []))
            groups = defaultdict(list)  # md5 -> [paths]
            digests = {}
            for p in paths:
                b = p.read_bytes()
                d = md5_bytes(b)
                digests[d] = sha256_bytes(b)
                groups[d].append(p)

            # canonical path per content
            canon = {}
            for d in groups:
                canon[d] = min(groups[d], key=lambda p: (pref_key(p), str(p)))

            # assign output names, resolving same-basename collisions
            used = {}
            assigned = {}  # md5 -> final name
            for d in sorted(groups, key=lambda k: (pref_key(canon[k]), str(canon[k]))):
                canon_path = canon[d]
                base = canon_path.name
                final = unique_name(base, used)
                used[final] = d
                assigned[d] = final

            # recommended copies + tagging
            for (calc, func), mix in RECOMMENDED_DIRS.items():
                mix_file = recommend.get((calc, func, elem))
                if not mix_file:
                    continue
                b = mix_file.read_bytes()
                d = md5_bytes(b)
                digest = digests.get(d)
                if digest is None:
                    digests[d] = sha256_bytes(b)
                    digest = digests[d]
                    groups[d].append(mix_file)
                    canon[d] = mix_file
                    assigned[d] = unique_name(mix_file.name, used)
                    used[assigned[d]] = d
                dst_dir = lib_out / elem / f"{calc}_{func}_recommended"
                dst_dir.mkdir(parents=True, exist_ok=True)
                dst = dst_dir / mix_file.name
                if not dst.exists():
                    dst.write_bytes(b)
                ti = tarfile.TarInfo(name=mix_file.name)
                ti.size = dst.stat().st_size
                with dst.open("rb") as fh:
                    tar_handles[(calc, func)].addfile(ti, fileobj=fh)

            # copies to all_upf_files
            elem_dir = lib_out / elem / "all_upf_files"
            elem_dir.mkdir(parents=True, exist_ok=True)
            file_records = []
            renames = []
            for d, final in assigned.items():
                src = canon[d]
                rel = src.relative_to(REPO_ROOT)
                top = rel.parts[0]
                if top == "libraries-pbe":
                    root = REPO_ROOT / "libraries-pbe"
                elif top == "libraries-pbesol":
                    root = REPO_ROOT / "libraries-pbesol"
                else:
                    root = REPO_ROOT / COMMON_ROOT
                lib = lib_folder(src, root)
                size = src.stat().st_size
                rec = {
                    "name": final,
                    "source_name": src.name,
                    "element": elem,
                    "functional": func,
                    "pp_type": infer_pp_type(src.name),
                    "z_valence": infer_z_valence(src.name),
                    "library": lib,
                    "source": str(rel),
                    "size": size,
                    "md5": d,
                    "sha256": digests[d],
                    "recommended": [],
                    "cutoffs": {},
                }
                csv_name = LIB_TO_CSV_NAME.get(lib)
                if csv_name and (elem, csv_name) in csv_rows:
                    r = csv_rows[(elem, csv_name)]
                    for k in ("ecut_eff", "ecut_prec", "eos_score", "score"):
                        v = r.get(k)
                        if v != "" and v is not None:
                            try:
                                rec["csv_" + k] = float(v)
                            except ValueError:
                                pass
                file_records.append(rec)
                if final != src.name:
                    renames.append({"name": src.name, "renamed_to": final, "source": str(rel)})
                dst = elem_dir / final
                if not dst.exists():
                    dst.write_bytes(src.read_bytes())

            # tag recommended files (by md5) with calc+cutoffs
            for (calc, func), mix in RECOMMENDED_DIRS.items():
                mix_file = recommend.get((calc, func, elem))
                if not mix_file:
                    continue
                d = md5_bytes(mix_file.read_bytes())
                tag = f"{calc}_{func}"
                for rec in file_records:
                    if rec["md5"] == d:
                        rec["recommended"].append(tag)
                        c = cutoffs[calc].get(elem)
                        if c:
                            rec["cutoffs"][calc] = c

            by_md5 = {r["md5"]: r for r in file_records}
            rec_md5 = {}
            for (calc, func), mix in RECOMMENDED_DIRS.items():
                mix_file = recommend.get((calc, func, elem))
                if mix_file:
                    rec_md5[f"{calc}_{func}"] = md5_bytes(mix_file.read_bytes())

            meta = {
                "element": elem,
                "count": len(file_records),
                "recommended": {
                    tag: by_md5.get(d, {}).get("name") for tag, d in rec_md5.items()
                },
                "files": file_records,
                "renamed": renames,
            }
            meta_out.joinpath(f"{elem}.json").write_text(json.dumps(meta, indent=2))

            summary["elements"] += 1
            summary["distinct_contents"] += len(file_records)
            summary["renamed"] += len(renames)
            summary["upf_files"] += sum(1 for r in file_records if r["name"].endswith(".upf"))
            print(f"  {elem}: {len(file_records)} contents ({len(renames)} renamed)")
    finally:
        for t in tar_handles.values():
            t.close()

    print(json.dumps(summary))


if __name__ == "__main__":
    sys.exit(main())