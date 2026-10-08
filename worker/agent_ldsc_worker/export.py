"""Strict conversion of validated task outputs into an ldsR parquet directory.

Joins are by verified identity, never position alone: each task's SNP column must equal the
bundle's regression index for that chromosome, and each annotation file must have exactly
the BIM row count. Chromosomes are processed numerically 1..22.
"""

from __future__ import annotations

import json
import os
import shutil
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq

from .ldscores import Task, read_task_annot, read_task_ld
from .refs import Bundle
from .util import AUTOSOMES, WorkerError, sha256_file, sha256_strings, write_json_atomic

HOME = Path(os.environ.get("AGENT_LDSC_HOME", Path(__file__).resolve().parents[2]))
VALIDATE_R = HOME / "R" / "validate_ldsR.R"


def export_ldsr(work: Path, out_tmp: Path, annotation_ids: list[str], tasks_by_annotation: dict,
                records: dict, bundle: Bundle) -> dict:
    """Write ld/annot/annot_ref parquet into out_tmp. Returns per-annotation counts."""
    reg = bundle.regression_snps()
    ref = bundle.reference_snps(columns=["CHR", "SNP", "common"])
    ld_cols, ref_cols, counts = {}, {}, {}
    for a in annotation_ids:
        ld_parts, annot_parts, m, m50 = [], [], 0, 0
        for chrom in AUTOSOMES:
            task: Task = tasks_by_annotation[a][chrom]
            ld = read_task_ld(work, task)
            exp = reg.SNP[reg.CHR == chrom].to_numpy()
            if len(ld) != len(exp) or not (ld.SNP.to_numpy() == exp).all():
                raise WorkerError("EXPORT_INVALID", f"{a} chr{chrom}: LD SNP order differs from regression index")
            ld_parts.append(ld[a].to_numpy(dtype=np.float64))
            annot = read_task_annot(work, task)
            n_ref = int((ref.CHR == chrom).sum())
            if len(annot) != n_ref or a not in annot.columns:
                raise WorkerError("EXPORT_INVALID", f"{a} chr{chrom}: annotation rows {len(annot)} != {n_ref}")
            annot_parts.append(annot[a].to_numpy())
            rec = records[task.name]
            idx = list(task.annotations).index(a)
            m += rec["M"][idx]
            m50 += rec["M_5_50"][idx]
        col = np.concatenate(ld_parts)
        member = np.concatenate(annot_parts).astype(np.int32)
        if not set(np.unique(member)) <= {0, 1}:
            raise WorkerError("EXPORT_INVALID", f"{a}: non-binary annotation")
        if int(member.sum()) != m or int(member[ref.common.to_numpy()].sum()) != m50 or m50 > m:
            raise WorkerError("EXPORT_INVALID", f"{a}: M/M_5_50 inconsistent with annotation membership",
                              {"m": m, "m50": m50, "member": int(member.sum())})
        if m == 0:
            raise WorkerError("INPUT_INVALID", f"{a}: annotation covers no reference SNPs genome-wide")
        ld_cols[a], ref_cols[a] = col, member
        counts[a] = {"m": m, "m50": m50,
                     "annotated_regression_snps": int(member[ref.SNP.isin(reg.SNP).to_numpy()].sum())}

    out_tmp.mkdir(parents=True, exist_ok=True)
    ld_table = pa.table({"SNP": pa.array(reg.SNP.to_numpy(), pa.string()),
                         **{a: pa.array(ld_cols[a], pa.float64()) for a in annotation_ids}})
    annot_table = pa.table({"annot": pa.array(annotation_ids, pa.string()),
                            "m": pa.array([float(counts[a]["m"]) for a in annotation_ids], pa.float64()),
                            "m50": pa.array([float(counts[a]["m50"]) for a in annotation_ids], pa.float64())})
    ref_table = pa.table({"SNP": pa.array(ref.SNP.to_numpy(), pa.string()),
                          **{a: pa.array(ref_cols[a], pa.int32()) for a in annotation_ids}})
    pq.write_table(ld_table, out_tmp / "ld.parquet")
    pq.write_table(annot_table, out_tmp / "annot.parquet")
    pq.write_table(ref_table, out_tmp / "annot_ref.parquet")
    return counts


def validate_ldsr_dir(d: Path, bundle: Bundle, annotation_ids: list[str] | None = None,
                      run_r: bool = True) -> dict:
    """Independent re-read of an exported directory plus the ldsR reader gate."""
    ld = pq.read_table(d / "ld.parquet")
    annot = pq.read_table(d / "annot.parquet").to_pandas()
    ref = pq.read_table(d / "annot_ref.parquet")
    names = ld.column_names[1:]
    problems = []
    if ld.column_names[0] != "SNP" or ref.column_names[0] != "SNP":
        problems.append("first column must be SNP")
    if annotation_ids is not None and names != annotation_ids:
        problems.append(f"ld columns {names} != requested {annotation_ids}")
    if annot.annot.tolist() != names or ref.column_names[1:] != names:
        problems.append("annotation names differ between ld/annot/annot_ref")
    snp = ld.column("SNP").to_numpy(zero_copy_only=False)
    if sha256_strings(snp) != bundle.meta["snp_universe"]["regression_order_sha256"]:
        problems.append("ld.parquet SNP order != bundle regression index")
    rsnp = ref.column("SNP").to_numpy(zero_copy_only=False)
    if sha256_strings(rsnp) != bundle.meta["snp_universe"]["reference_order_sha256"]:
        problems.append("annot_ref.parquet SNP order != bundle reference index")
    common = bundle.reference_snps(columns=["common"]).common.to_numpy()
    for i, a in enumerate(names):
        v = ld.column(a).to_numpy()
        if not np.isfinite(v).all():
            problems.append(f"{a}: non-finite LD scores")
        mem = ref.column(a).to_numpy()
        if int(mem.sum()) != annot.m[i] or int(mem[common].sum()) != annot.m50[i] or annot.m50[i] > annot.m[i]:
            problems.append(f"{a}: m/m50 do not match annot_ref membership and common mask")
    result = {"python_checks": "pass" if not problems else "fail", "problems": problems,
              "regression_snps": len(snp), "reference_snps": len(rsnp), "annotations": names}
    if run_r:
        mask_compatible = "ldsr_common_mask_identical" in bundle.meta.get("checks", {})
        mask = "bundled" if mask_compatible else str(bundle.path("reference_snps"))
        result["ldsR"] = run_r_validation(d, mask)
        result["ldsr_mask_compatible"] = mask_compatible
        if result["ldsR"].get("status") != "pass":
            problems.append("ldsR validation failed")
    if problems:
        raise WorkerError("EXPORT_INVALID", f"{d}: export validation failed", result)
    return result


def run_r_validation(d: Path, mask: str = "bundled") -> dict:
    rscript = shutil.which("Rscript")
    if rscript is None:
        raise WorkerError("EXPORT_INVALID", "Rscript not found; ldsR gate cannot run")
    proc = subprocess.run([rscript, str(VALIDATE_R), str(d), mask], capture_output=True, text=True)
    lines = [l for l in proc.stdout.splitlines() if l.startswith("{")]
    if proc.returncode != 0 or not lines:
        return {"status": "fail", "stderr": proc.stderr[-3000:]}
    return json.loads(lines[-1])


def file_hashes(d: Path) -> dict:
    return {p.name: sha256_file(p) for p in sorted(d.iterdir()) if p.is_file()}
