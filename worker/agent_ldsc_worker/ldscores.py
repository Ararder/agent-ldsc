"""One LD-score task = (annotation batch, chromosome), run with the pinned legacy ldsc.py.

A task writes into a temporary directory, validates the outputs against the reference
bundle and its own annotation, then publishes the directory with COMPLETE.json. A task is
reused on resume only if COMPLETE.json's input key matches and output hashes still verify.
"""

from __future__ import annotations

import gzip
import os
import re
import shlex
import shutil
import subprocess
import threading
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

from .intervals import snp_membership
from .refs import Bundle, read_bim
from .util import (WorkerError, canonical_json, now, read_json, sha256_bytes, sha256_file,
                   write_json_atomic)

LDSC_CMD = shlex.split(os.environ.get("AGENT_LDSC_LDSC_CMD", "/opt/envs/ldsc/bin/python /opt/ldsc/ldsc.py"))
LDSC_COMMIT = os.environ.get("AGENT_LDSC_LDSC_COMMIT", "unknown")
OUTPUT_SUFFIXES = (".l2.ldscore.gz", ".l2.M", ".l2.M_5_50", ".log")

_children: set[subprocess.Popen] = set()
_children_lock = threading.Lock()


def terminate_children() -> None:
    with _children_lock:
        for p in list(_children):
            if p.poll() is None:
                p.terminate()


@dataclass(frozen=True)
class Task:
    key: str                      # directory name for the annotation batch
    annotations: tuple[str, ...]
    chrom: int

    @property
    def name(self) -> str:
        return f"{self.key}/chr{self.chrom}"


def task_dir(work: Path, task: Task) -> Path:
    return work / "tasks" / task.key / f"chr{task.chrom}"


def input_key(task: Task, ann_hashes: dict, bundle: Bundle, ld_params: dict) -> str:
    return sha256_bytes(canonical_json({
        "annotations": [[a, ann_hashes[a]] for a in task.annotations],
        "chrom": task.chrom,
        "bundle_id": bundle.id,
        "reference_order_sha256": bundle.meta["snp_universe"]["reference_order_sha256"],
        "regression_order_sha256": bundle.meta["snp_universe"]["regression_order_sha256"],
        "plink_sha256": [bundle.meta["files"][f"{bundle.meta['layout']['plink_prefix']}.{task.chrom}.{e}"]["sha256"]
                         for e in ("bed", "bim", "fam")],
        "ldsc_commit": LDSC_COMMIT,
        "ld_params": ld_params,
    }))


def is_complete(work: Path, task: Task, key: str) -> bool:
    d = task_dir(work, task)
    rec_path = d / "COMPLETE.json"
    if not rec_path.exists():
        return False
    try:
        rec = read_json(rec_path)
    except ValueError:
        return False
    if rec.get("input_key") != key:
        return False
    return all((d / f).is_file() and sha256_file(d / f) == h for f, h in rec["outputs"].items())


def run_task(work: Path, task: Task, intervals: dict, bundle: Bundle, ld_params: dict,
             key: str, threads: int = 1) -> dict:
    if not bundle.root.is_absolute():
        raise ValueError("bundle root must be absolute (ldsc runs inside the task directory)")
    final = task_dir(Path(work).resolve(), task)
    tmp = final.parent / f".tmp-chr{task.chrom}-{os.getpid()}-{threading.get_ident()}"
    shutil.rmtree(tmp, ignore_errors=True)
    tmp.mkdir(parents=True)
    try:
        bim = read_bim(Path(f"{bundle.plink(task.chrom)}.bim"))
        membership = np.column_stack([snp_membership(intervals[a], task.chrom, bim.BP.to_numpy())
                                      for a in task.annotations])
        annot_path = tmp / f"annot.{task.chrom}.annot.gz"
        _write_thin_annot(annot_path, task.annotations, membership)

        out_prefix = tmp / f"ld.{task.chrom}"
        cmd = LDSC_CMD + [
            "--bfile", str(bundle.plink(task.chrom)),
            "--ld-wind-cm", str(ld_params["ld_wind_cm"]),
            "--annot", str(annot_path), "--thin-annot",
            "--out", str(out_prefix),
            "--print-snps", str(bundle.path("print_snps")),
        ]
        env = dict(os.environ, OPENBLAS_NUM_THREADS=str(threads), OMP_NUM_THREADS=str(threads),
                   MKL_NUM_THREADS=str(threads))
        with open(tmp / "ldsc.stdout", "w") as so:
            p = subprocess.Popen(cmd, stdout=so, stderr=subprocess.STDOUT, env=env, cwd=tmp)
            with _children_lock:
                _children.add(p)
            try:
                rc = p.wait()
            finally:
                with _children_lock:
                    _children.discard(p)
        if rc != 0:
            raise WorkerError("LDSC_TASK_FAILED", f"{task.name}: ldsc exited {rc}",
                              _tail(tmp / "ldsc.stdout"))
        result = validate_task_outputs(tmp, task, bim, membership, bundle)
        outputs = {f.name: sha256_file(f) for f in sorted(tmp.iterdir()) if f.is_file()}
        record = {"task": task.name, "annotations": list(task.annotations), "chrom": task.chrom,
                  "input_key": key, "command": cmd, "outputs": outputs, "finished": now(), **result}
        write_json_atomic(tmp / "COMPLETE.json", record)
        if final.exists():
            shutil.rmtree(final)
        os.replace(tmp, final)
        return record
    except BaseException:
        shutil.rmtree(tmp, ignore_errors=True)
        raise


def _write_thin_annot(path: Path, names, membership: np.ndarray) -> None:
    with gzip.GzipFile(path, "wb", mtime=0) as f:
        f.write(("\t".join(names) + "\n").encode())
        f.write("".join("\t".join(map(str, r)) + "\n" for r in membership.tolist()).encode())


def validate_task_outputs(tmp: Path, task: Task, bim: pd.DataFrame, membership: np.ndarray,
                          bundle: Bundle) -> dict:
    prefix = tmp / f"ld.{task.chrom}"
    for suf in OUTPUT_SUFFIXES:
        if not Path(f"{prefix}{suf}").is_file():
            raise WorkerError("LDSC_TASK_FAILED", f"{task.name}: missing output {prefix.name}{suf}")
    log = Path(f"{prefix}.log").read_text()
    if "Error parsing .annot file" in log:
        raise WorkerError("LDSC_TASK_FAILED", f"{task.name}: ldsc could not parse the annotation")
    m = re.search(r"Read (\d+) annotations for (\d+) SNPs", log)
    if not m or int(m.group(1)) != len(task.annotations) or int(m.group(2)) != len(bim):
        raise WorkerError("LDSC_TASK_FAILED", f"{task.name}: annotation read mismatch in ldsc log")
    stdout = (tmp / "ldsc.stdout").read_text()
    kept = re.search(r"After filtering, (\d+) SNPs remain", stdout)
    if not kept or int(kept.group(1)) != len(bim):
        raise WorkerError("REFERENCE_MISMATCH",
                          f"{task.name}: ldsc kept {kept.group(1) if kept else 'unknown'} of {len(bim)} SNPs "
                          "(monomorphic filter?)")

    ld = pd.read_csv(f"{prefix}.l2.ldscore.gz", sep="\t")
    reg = bundle.regression_snps()
    reg = reg[reg.CHR == task.chrom]
    if len(ld) != len(reg) or not (ld.SNP.to_numpy() == reg.SNP.to_numpy()).all():
        raise WorkerError("LDSC_TASK_FAILED", f"{task.name}: output SNPs differ from the regression index")
    cols = _ld_columns(task.annotations)
    if list(ld.columns) != ["CHR", "SNP", "BP"] + cols:
        raise WorkerError("LDSC_TASK_FAILED", f"{task.name}: unexpected columns {list(ld.columns)}")
    vals = ld[cols].to_numpy(dtype=float)
    if not np.isfinite(vals).all():
        raise WorkerError("LDSC_TASK_FAILED", f"{task.name}: non-finite LD scores")

    M = _read_counts(Path(f"{prefix}.l2.M"))
    M550 = _read_counts(Path(f"{prefix}.l2.M_5_50"))
    common = bundle.reference_snps(columns=["CHR", "common"])
    common = common.loc[common.CHR == task.chrom, "common"].to_numpy()
    exp_m = membership.sum(axis=0)
    exp_m50 = membership[common].sum(axis=0)
    if len(M) != len(task.annotations) or not np.array_equal(M, exp_m):
        raise WorkerError("LDSC_TASK_FAILED", f"{task.name}: .l2.M {M} != annotated SNPs {exp_m.tolist()}")
    if not np.array_equal(M550, exp_m50):
        raise WorkerError("REFERENCE_MISMATCH",
                          f"{task.name}: .l2.M_5_50 {M550} != frq-based common count {exp_m50.tolist()}")
    return {"n_regression_snps": len(ld), "n_reference_snps": len(bim),
            "M": [int(x) for x in M], "M_5_50": [int(x) for x in M550]}


def _ld_columns(annotations) -> list[str]:
    # ldsc names a single-annotation column "L2", otherwise "<name>L2"
    return ["L2"] if len(annotations) == 1 else [f"{a}L2" for a in annotations]


def _read_counts(path: Path) -> np.ndarray:
    vals = path.read_text().split()
    return np.array([int(round(float(v))) for v in vals], dtype=np.int64)


def _tail(path: Path, n: int = 4000) -> str:
    try:
        return path.read_text()[-n:]
    except FileNotFoundError:
        return ""


def read_task_ld(work: Path, task: Task) -> pd.DataFrame:
    ld = pd.read_csv(task_dir(work, task) / f"ld.{task.chrom}.l2.ldscore.gz", sep="\t")
    return ld.rename(columns=dict(zip(_ld_columns(task.annotations), task.annotations)))


def read_task_annot(work: Path, task: Task) -> pd.DataFrame:
    return pd.read_csv(task_dir(work, task) / f"annot.{task.chrom}.annot.gz", sep="\t", dtype=np.int8)
