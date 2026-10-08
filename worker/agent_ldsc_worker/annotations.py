"""Normalize BED intervals or human gene identifiers to canonical reference-build intervals.

Every input record ends up in an audit table with a status and reason. Records are dropped
only under a reason category that the request lists in allow_drop; otherwise the request
fails with MAPPING_LOSS. Ambiguous gene identifiers always fail.
"""

from __future__ import annotations

import re
import shutil
import subprocess
import tempfile
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd

from .intervals import merge_intervals, per_chrom_summary, total_bp
from .refs import Bundle
from .util import WorkerError, normalize_chrom, sha256_bytes

BUILD_ALIASES = {"grch37": "GRCh37", "hg19": "GRCh37", "grch38": "GRCh38", "hg38": "GRCh38"}
DROP_CATEGORIES = {"non_autosomal", "liftover_unmapped", "gene_unmapped"}
ENSG_RE = re.compile(r"^(ENSG\d{11})(\.\d+)?$")
MAX_EXAMPLES = 20


@dataclass
class NormalizedAnnotation:
    id: str
    intervals: pd.DataFrame          # merged chrom/start/end in reference build
    audit: pd.DataFrame
    summary: dict = field(default_factory=dict)

    def canonical_bed(self) -> bytes:
        lines = [f"chr{c}\t{s}\t{e}\n" for c, s, e in self.intervals[["chrom", "start", "end"]].itertuples(index=False)]
        return "".join(lines).encode()

    @property
    def content_hash(self) -> str:
        return sha256_bytes(self.canonical_bed())


def normalize(spec: dict, input_path: Path, bundle: Bundle, profile: dict) -> NormalizedAnnotation:
    if spec.get("species", "human") != "human":
        raise WorkerError("BUILD_UNSUPPORTED", f"{spec['id']}: only human input is supported in v1")
    if spec["type"] == "bed":
        ann = _normalize_bed(spec, input_path, bundle, profile)
    elif spec["type"] == "genes":
        ann = _normalize_genes(spec, input_path, bundle, profile)
    else:
        raise WorkerError("INPUT_INVALID", f"unknown annotation type {spec['type']}")
    _enforce_drop_policy(spec, ann)
    if len(ann.intervals) == 0:
        raise WorkerError("INPUT_INVALID", f"{spec['id']}: no intervals remain after normalization")
    return ann


# --------------------------------------------------------------------------- BED


def parse_bed(path: Path) -> pd.DataFrame:
    rows, bad = [], []
    with open(path) as f:
        for n, line in enumerate(f, start=1):
            line = line.rstrip("\r\n")
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            p = line.split("\t") if "\t" in line else line.split()
            if len(p) < 3:
                bad.append((n, "fewer than 3 columns"))
                continue
            try:
                s, e = int(p[1]), int(p[2])
            except ValueError:
                bad.append((n, "non-integer start/end"))
                continue
            if s < 0:
                bad.append((n, "negative start"))
            elif e <= s:
                bad.append((n, "end <= start"))
            else:
                rows.append((n, p[0], s, e, p[3] if len(p) > 3 else ""))
    if bad:
        raise WorkerError("INPUT_INVALID", f"{path.name}: {len(bad)} malformed BED lines",
                          [{"line": n, "problem": why} for n, why in bad[:MAX_EXAMPLES]])
    if not rows:
        raise WorkerError("INPUT_INVALID", f"{path.name}: no BED records")
    return pd.DataFrame(rows, columns=["input_line", "input_chrom", "input_start", "input_end", "input_name"])


def _normalize_bed(spec, path, bundle, profile) -> NormalizedAnnotation:
    build = BUILD_ALIASES.get(str(spec.get("build", "")).lower())
    if build is None:
        raise WorkerError("BUILD_UNSUPPORTED", f"{spec['id']}: build must be declared as GRCh37/hg19 or GRCh38/hg38")
    ref_build = bundle.meta["genome_build"]
    bed = parse_bed(path)
    bed["chrom_num"] = [normalize_chrom(c) for c in bed.input_chrom]
    bed["status"], bed["reason"] = "mapped", ""
    bed["chrom"], bed["start"], bed["end"] = 0, -1, -1
    non_auto = bed.chrom_num.isna()
    bed.loc[non_auto, ["status", "reason"]] = ["dropped", "non_autosomal"]

    if build == ref_build:
        ok = ~non_auto
        bed.loc[ok, "chrom"] = bed.loc[ok, "chrom_num"].astype(int)
        bed.loc[ok, "start"] = bed.loc[ok, "input_start"]
        bed.loc[ok, "end"] = bed.loc[ok, "input_end"]
        liftover = None
    elif build == "GRCh38" and ref_build == "GRCh37":
        liftover = _liftover(bed[~non_auto], bundle.path("chain_hg38_to_ref"), profile["liftover"])
        bed = bed.drop(columns=["chrom", "start", "end"]).merge(liftover, on="input_line", how="left")
        lifted = bed.lift_status.notna()
        bed.loc[lifted & (bed.lift_status != "mapped"), "status"] = "dropped"
        bed.loc[lifted & (bed.lift_status != "mapped"), "reason"] = bed.lift_status
        bed = bed.drop(columns="lift_status").fillna({"chrom": 0, "start": -1, "end": -1})
        liftover = {"chain": bundle.meta["layout"]["chain_hg38_to_ref"],
                    "chain_sha256": bundle.meta["files"][bundle.meta["layout"]["chain_hg38_to_ref"]]["sha256"],
                    **profile["liftover"]}
    else:
        raise WorkerError("BUILD_UNSUPPORTED", f"{spec['id']}: no supported mapping {build} -> {ref_build}")

    bed[["chrom", "start", "end"]] = bed[["chrom", "start", "end"]].astype(np.int64)
    kept = bed[bed.status == "mapped"]
    _check_bounds(kept, bundle, spec["id"])
    merged = merge_intervals(kept[["chrom", "start", "end"]])
    autosomal_in = bed[~non_auto]
    summary = {
        "type": "bed", "input_build": build, "reference_build": ref_build,
        "records_in": len(bed), "records_kept": len(kept),
        "records_dropped": _reason_counts(bed),
        "input_bp": int((bed.input_end - bed.input_start).sum()),
        "input_autosomal_bp": int((autosomal_in.input_end - autosomal_in.input_start).sum()),
        "kept_bp_before_merge": int((kept.end - kept.start).sum()),
        "kept_input_bp": int((kept.input_end - kept.input_start).sum()),
        "merged_intervals": len(merged), "merged_bp": total_bp(merged),
        "liftover": liftover, "per_chromosome": per_chrom_summary(merged),
    }
    audit = bed[["input_line", "input_chrom", "input_start", "input_end", "input_name",
                 "status", "reason", "chrom", "start", "end"]]
    return NormalizedAnnotation(spec["id"], merged, audit, summary)


def _liftover(bed: pd.DataFrame, chain: Path, params: dict) -> pd.DataFrame:
    exe = shutil.which("liftOver")
    if exe is None:
        raise WorkerError("BUILD_UNSUPPORTED", "liftOver binary not found in image")
    with tempfile.TemporaryDirectory() as td:
        td = Path(td)
        with open(td / "in.bed", "w") as f:
            for r in bed.itertuples(index=False):
                f.write(f"chr{int(r.chrom_num)}\t{r.input_start}\t{r.input_end}\tr{r.input_line}\n")
        args = [exe, f"-minMatch={params['min_match']}"]
        if params.get("multiple"):
            raise WorkerError("BUILD_UNSUPPORTED", "multiple liftover mappings are not supported")
        args += [str(td / "in.bed"), str(chain), str(td / "out.bed"), str(td / "unmapped.bed")]
        proc = subprocess.run(args, capture_output=True, text=True)
        if proc.returncode != 0:
            raise WorkerError("BUILD_UNSUPPORTED", "liftOver failed", proc.stderr[-2000:])
        out = pd.read_csv(td / "out.bed", sep="\t", header=None, names=["c", "start", "end", "name"],
                          dtype={"c": str, "name": str}) if (td / "out.bed").stat().st_size else \
            pd.DataFrame(columns=["c", "start", "end", "name"])
        reasons = {}
        reason = "liftover_unmapped"
        for line in open(td / "unmapped.bed"):
            if line.startswith("#"):
                reason = "liftover_" + line[1:].strip().lower().replace(" in new", "").replace(" ", "_")
            elif line.strip():
                reasons[int(line.split("\t")[3][1:])] = reason
    if out.name.duplicated().any():
        raise WorkerError("BUILD_UNSUPPORTED", "liftOver produced multiple mappings for a record")
    res = pd.DataFrame({"input_line": out.name.str[1:].astype(int),
                        "chrom": [normalize_chrom(c) for c in out.c],
                        "start": out.start.astype(np.int64), "end": out.end.astype(np.int64)})
    res["lift_status"] = np.where(res.chrom.isna(), "liftover_non_autosomal_target", "mapped")
    res["chrom"] = res.chrom.fillna(0).astype(int)
    lost = pd.DataFrame({"input_line": list(reasons), "lift_status": list(reasons.values())})
    lost["chrom"], lost["start"], lost["end"] = 0, -1, -1
    missing = set(bed.input_line) - set(res.input_line) - set(lost.input_line)
    if missing:
        raise WorkerError("BUILD_UNSUPPORTED", f"liftOver lost {len(missing)} records without a reason")
    return pd.concat([res, lost], ignore_index=True)


def _check_bounds(kept: pd.DataFrame, bundle: Bundle, aid: str) -> None:
    sizes = pd.read_csv(bundle.path("chrom_sizes"), sep="\t")
    length = {normalize_chrom(c): int(n) for c, n in zip(sizes.chrom, sizes.length)}
    over = kept[kept.end > kept.chrom.map(length)]
    if len(over):
        raise WorkerError("INPUT_INVALID", f"{aid}: {len(over)} intervals extend past chromosome ends "
                          "(wrong declared build?)", over.head(MAX_EXAMPLES).to_dict("records"))


# --------------------------------------------------------------------------- genes


def read_gene_list(path: Path, header: bool) -> list[tuple[int, str]]:
    items = []
    for n, line in enumerate(open(path), start=1):
        v = line.strip()
        if not v or v.startswith("#"):
            continue
        items.append((n, v))
    if header:
        items = items[1:]
    if not items:
        raise WorkerError("INPUT_INVALID", f"{path.name}: no gene identifiers")
    bad = [(n, v) for n, v in items if len(v.split()) != 1]
    if bad:
        raise WorkerError("INPUT_INVALID", f"{path.name}: expected one identifier per line",
                          [{"line": n, "value": v} for n, v in bad[:MAX_EXAMPLES]])
    return items


def _normalize_genes(spec, path, bundle, profile) -> NormalizedAnnotation:
    gp = profile["genes"]
    if spec.get("gene_model", "gene_span") != "gene_span":
        raise WorkerError("INPUT_INVALID", f"{spec['id']}: gene_model must be gene_span")
    if "flank_bp" not in spec:
        raise WorkerError("INPUT_INVALID", f"{spec['id']}: flank_bp is required for gene annotations")
    flank = int(spec["flank_bp"])
    if not 0 <= flank <= gp["max_flank_bp"]:
        raise WorkerError("INPUT_INVALID", f"{spec['id']}: flank_bp outside [0, {gp['max_flank_bp']}]")
    id_type = spec["id_type"]
    genes = pd.read_csv(bundle.path("genes"), sep="\t", dtype=str, keep_default_na=False)
    hgnc = pd.read_csv(bundle.path("hgnc"), sep="\t", dtype=str, keep_default_na=False)
    by_gid = {g: r for g, r in genes[~genes.gene_id.str.endswith("_PAR_Y")].groupby("gene_id")}
    resolve = _symbol_resolver(hgnc) if id_type == "symbol" else None
    if id_type not in ("symbol", "ensembl_gene"):
        raise WorkerError("INPUT_INVALID", f"{spec['id']}: id_type must be symbol or ensembl_gene")

    rows, seen = [], {}
    for line_no, value in read_gene_list(path, bool(spec.get("header", False))):
        rec = {"input_line": line_no, "input_value": value, "status": "mapped", "reason": "",
               "hgnc_id": "", "gene_id": "", "gene_name": "", "gene_type": "", "candidates": "",
               "chrom": 0, "start": -1, "end": -1}
        if id_type == "ensembl_gene":
            m = ENSG_RE.match(value)
            gid = m.group(1) if m else None
            if gid is None:
                rec.update(status="dropped", reason="gene_invalid_id")
        else:
            gid, hgnc_id, why, cands = resolve(value)
            rec["hgnc_id"] = hgnc_id or ""
            if why == "ambiguous":
                rec.update(status="ambiguous", reason="gene_ambiguous_symbol", candidates="|".join(cands))
            elif why:
                rec.update(status="dropped", reason=why)
        if rec["status"] == "mapped":
            if gid in seen:
                rec.update(status="duplicate", reason=f"duplicate_of_line_{seen[gid]}", gene_id=gid)
                rows.append(rec)
                continue
            seen[gid] = line_no
            rec["gene_id"] = gid
            g = by_gid.get(gid)
            if g is None:
                rec.update(status="dropped", reason="gene_not_in_gene_model")
            elif len(g) != 1:
                rec.update(status="dropped", reason="gene_multiple_records")
            else:
                g = g.iloc[0]
                rec.update(gene_name=g.gene_name, gene_type=g.gene_type)
                chrom = normalize_chrom(g.chrom)
                if chrom is None:
                    rec.update(status="dropped", reason="non_autosomal")
                elif g.remap_status not in gp["allowed_remap_status"] or \
                        g.remap_num_mappings not in gp["allowed_remap_num_mappings"]:
                    rec.update(status="dropped", reason=f"gene_remap_{g.remap_status}_{g.remap_num_mappings or 'na'}")
                else:
                    rec.update(chrom=chrom, start=int(g.start) - 1, end=int(g.end))
        rows.append(rec)
    audit = pd.DataFrame(rows)

    amb = audit[audit.status == "ambiguous"]
    if len(amb):
        raise WorkerError("GENE_AMBIGUOUS", f"{spec['id']}: {len(amb)} ambiguous gene symbols",
                          amb[["input_line", "input_value", "candidates"]].head(MAX_EXAMPLES).to_dict("records"))

    sizes = pd.read_csv(bundle.path("chrom_sizes"), sep="\t")
    length = {normalize_chrom(c): int(n) for c, n in zip(sizes.chrom, sizes.length)}
    ok = audit.status == "mapped"
    audit["span_start"] = np.where(ok, np.maximum(0, audit.start - flank), -1)
    audit["span_end"] = np.where(ok, [min(length.get(c, 0), e + flank) if o else -1
                                      for c, e, o in zip(audit.chrom, audit.end, ok)], -1)
    kept = audit[ok]
    merged = merge_intervals(kept.rename(columns={"start": "s0", "end": "e0", "span_start": "start",
                                                  "span_end": "end"})[["chrom", "start", "end"]])
    summary = {
        "type": "genes", "id_type": id_type, "gene_model": "gene_span", "flank_bp": flank,
        "reference_build": bundle.meta["genome_build"],
        "records_in": len(audit), "genes_mapped": int(ok.sum()),
        "duplicates": int((audit.status == "duplicate").sum()),
        "records_dropped": _reason_counts(audit),
        "gene_span_bp_before_flank": int((kept.end - kept.start).sum()),
        "merged_intervals": len(merged), "merged_bp": total_bp(merged),
        "gene_coordinates": bundle.meta["genes"]["coordinate_source"],
        "gene_identifiers": bundle.meta["genes"]["identifier_source"],
        "per_chromosome": per_chrom_summary(merged),
    }
    return NormalizedAnnotation(spec["id"], merged, audit, summary)


def _symbol_resolver(hgnc: pd.DataFrame):
    approved = hgnc[hgnc.status == "Approved"]
    by_symbol = dict(zip(approved.symbol, zip(approved.hgnc_id, approved.ensembl_gene_id)))
    alt: dict[str, set] = {}
    for col in ("alias_symbol", "prev_symbol"):
        for hid, vals in zip(approved.hgnc_id, approved[col]):
            for v in filter(None, vals.split("|")):
                alt.setdefault(v, set()).add(hid)
    ens = dict(zip(approved.hgnc_id, approved.ensembl_gene_id))

    def resolve(symbol: str):
        if symbol in by_symbol:
            hid, gid = by_symbol[symbol]
            return (gid or None), hid, ("" if gid else "gene_no_ensembl_id"), []
        cands = sorted(alt.get(symbol, ()))
        if len(cands) == 1:
            gid = ens[cands[0]]
            return (gid or None), cands[0], ("" if gid else "gene_no_ensembl_id"), []
        if len(cands) > 1:
            return None, None, "ambiguous", [f"{h}:{approved.symbol[approved.hgnc_id == h].iloc[0]}" for h in cands]
        return None, None, "gene_symbol_not_found", []

    return resolve


# --------------------------------------------------------------------------- policy


def _reason_category(reason: str) -> str:
    if reason == "non_autosomal":
        return "non_autosomal"
    if reason.startswith("liftover_"):
        return "liftover_unmapped"
    return "gene_unmapped"


def _reason_counts(audit: pd.DataFrame) -> dict:
    d = audit[audit.status == "dropped"]
    return {str(k): int(v) for k, v in d.reason.value_counts().sort_index().items()}


def _enforce_drop_policy(spec: dict, ann: NormalizedAnnotation) -> None:
    allowed = set(spec.get("allow_drop", []))
    unknown = allowed - DROP_CATEGORIES
    if unknown:
        raise WorkerError("INPUT_INVALID", f"{spec['id']}: unknown allow_drop categories {sorted(unknown)}")
    dropped = ann.audit[ann.audit.status == "dropped"]
    cats = dropped.reason.map(_reason_category)
    blocked = dropped[~cats.isin(allowed)]
    ann.summary["allow_drop"] = sorted(allowed)
    if len(blocked):
        raise WorkerError("MAPPING_LOSS",
                          f"{spec['id']}: {len(blocked)} records would be dropped; add the category to allow_drop to accept",
                          {"by_reason": _reason_counts(blocked),
                           "categories": sorted(set(blocked.reason.map(_reason_category))),
                           "examples": blocked.head(MAX_EXAMPLES).astype(str).to_dict("records")})
