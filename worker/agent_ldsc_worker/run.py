"""Run orchestration: validate -> normalize -> per-chromosome LD tasks -> export -> publish.

Layout under the run directory (see docs/input-output-contract.md):
  input/request.json, input/<files>   staged by the launcher, read-only here
  resolved-request.json, state.json, provenance/, normalized/, tasks/, logs/, output/ldsR/
  COMPLETE.json                        written last; its presence means output/ldsR is final
"""

from __future__ import annotations

import json
import os
import platform
import shutil
import signal
import socket
import subprocess
import traceback
from concurrent.futures import FIRST_COMPLETED, ThreadPoolExecutor, wait
from pathlib import Path

import jsonschema
import pandas as pd

from . import __version__
from .annotations import NormalizedAnnotation, normalize
from .export import export_ldsr, file_hashes, validate_ldsr_dir
from .ldscores import LDSC_COMMIT, Task, input_key, is_complete, run_task, terminate_children
from .refs import Bundle, load_bundle
from .util import (AUTOSOMES, WorkerError, canonical_json, dir_lock, now, read_json, sha256_bytes,
                   sha256_file, write_json_atomic)

HOME = Path(os.environ.get("AGENT_LDSC_HOME", Path(__file__).resolve().parents[2]))
RESERVED = {"snp", "chr", "bp", "cm", "annot", "m", "m50", "maf", "common", "base", "l2"}


class Interrupted(Exception):
    pass


def _raise_interrupt(signum, frame):
    raise Interrupted(signal.Signals(signum).name)


class State:
    def __init__(self, work: Path):
        self.path = work / "state.json"
        self.data = read_json(self.path) if self.path.exists() else {"created": now(), "attempts": []}

    def update(self, **kw) -> None:
        self.data.update(kw, updated=now())
        write_json_atomic(self.path, self.data)


class Log:
    def __init__(self, path: Path):
        path.parent.mkdir(parents=True, exist_ok=True)
        self.f = open(path, "a", buffering=1)

    def __call__(self, msg: str) -> None:
        line = f"[{now()}] {msg}"
        print(line, flush=True)
        self.f.write(line + "\n")


def load_profile(profile_id: str) -> dict:
    p = HOME / "profiles" / "scientific" / f"{profile_id}.json"
    if not p.exists():
        raise WorkerError("INPUT_INVALID", f"unknown scientific profile {profile_id}")
    prof = read_json(p)
    if prof.get("profile_id") != profile_id:
        raise WorkerError("INPUT_INVALID", f"profile file {p} has id {prof.get('profile_id')}")
    return prof


def validate_request(request: dict, input_dir: Path) -> None:
    schema = read_json(HOME / "schemas" / "request.schema.json")
    try:
        jsonschema.validate(request, schema)
    except jsonschema.ValidationError as e:
        raise WorkerError("INPUT_INVALID", f"request schema: {e.message}", list(e.absolute_path)) from None
    ids = [a["id"] for a in request["annotations"]]
    if len(set(ids)) != len(ids):
        raise WorkerError("INPUT_INVALID", "duplicate annotation ids")
    if len({i.lower() for i in ids}) != len(ids):
        raise WorkerError("INPUT_INVALID", "annotation ids differ only by case")
    bad = [i for i in ids if i.lower() in RESERVED or i.endswith("L2")]
    if bad:
        raise WorkerError("INPUT_INVALID", f"reserved annotation ids: {bad}")
    root = input_dir.resolve()
    for a in request["annotations"]:
        p = (input_dir / a["file"]).resolve()
        if root not in p.parents or not p.is_file():
            raise WorkerError("INPUT_INVALID", f"{a['id']}: input file {a['file']} not staged under input/")


def run(request_path: Path, bundle_dir: Path, work: Path, jobs: int = 1,
        annotations_per_task: int = 1, force_unlock: bool = False) -> int:
    # ldsc runs with cwd = task directory, so every path handed down must be absolute
    request_path, bundle_dir, work = Path(request_path).resolve(), Path(bundle_dir).resolve(), Path(work).resolve()
    work.mkdir(parents=True, exist_ok=True)
    lock = work / ".run.lock"
    if force_unlock and lock.exists():
        shutil.rmtree(lock)
    log = Log(work / "logs" / "worker.log")
    try:
        with dir_lock(lock, timeout=0, stale_after=0):
            return _run_locked(request_path, bundle_dir, work, jobs, annotations_per_task, log)
    except WorkerError as e:
        if e.code == "RUN_LOCKED":
            log(f"ERROR {e}")
            return 3
        raise


def _run_locked(request_path, bundle_dir, work, jobs, per_task, log) -> int:
    state = State(work)
    attempt = {"started": now(), "host": socket.gethostname(), "slurm_job_id": os.environ.get("SLURM_JOB_ID"),
               "jobs": jobs, "image": os.environ.get("AGENT_LDSC_IMAGE", "unknown")}
    state.data["attempts"].append(attempt)
    state.update(status="initializing", stage="validate", error=None)
    old = {s: signal.signal(s, _raise_interrupt) for s in (signal.SIGTERM, signal.SIGINT, signal.SIGUSR1)}
    try:
        code = _pipeline(request_path, bundle_dir, work, jobs, per_task, state, log)
        attempt["finished"], attempt["result"] = now(), "completed"
        state.update(status="completed", stage="done")
        return code
    except WorkerError as e:
        log(f"ERROR {e}")
        write_json_atomic(work / "logs" / "error.json", e.to_json())
        attempt["finished"], attempt["result"] = now(), e.code
        state.update(status="failed", error=e.to_json())
        return 2
    except Interrupted as e:
        terminate_children()
        log(f"interrupted by {e}; completed tasks are kept for resume")
        attempt["finished"], attempt["result"] = now(), "interrupted"
        state.update(status="cancelled", error={"code": "RESOURCE_LIMIT", "message": f"interrupted by {e}"})
        return 4
    except Exception as e:  # unexpected: record and surface
        terminate_children()
        log("ERROR unexpected\n" + traceback.format_exc())
        attempt["finished"], attempt["result"] = now(), "crash"
        state.update(status="failed", error={"code": "STATE_INVALID", "message": repr(e)})
        return 5
    finally:
        for s, h in old.items():
            signal.signal(s, h)


def _pipeline(request_path, bundle_dir, work, jobs, per_task, state, log) -> int:
    input_dir = request_path.parent
    request = read_json(request_path)
    validate_request(request, input_dir)
    profile = load_profile(request["scientific_profile"])
    log(f"verifying reference bundle {bundle_dir}")
    bundle = load_bundle(bundle_dir, verify=True)
    if bundle.id != request["reference_id"]:
        raise WorkerError("REFERENCE_MISMATCH", f"request wants {request['reference_id']}, bundle is {bundle.id}")
    ids = [a["id"] for a in request["annotations"]]

    state.update(stage="normalize")
    normalized: dict[str, NormalizedAnnotation] = {}
    input_hashes = {}
    for spec in request["annotations"]:
        path = input_dir / spec["file"]
        input_hashes[spec["id"]] = sha256_file(path)
        log(f"normalizing {spec['id']} ({spec['type']})")
        ann = normalize(spec, path, bundle, profile)
        normalized[spec["id"]] = ann
        _write_normalized(work / "normalized", ann)
        log(f"  {ann.summary['merged_intervals']} merged intervals, {ann.summary['merged_bp']} bp; "
            f"dropped {ann.summary['records_dropped']}")
    ann_hashes = {a: normalized[a].content_hash for a in ids}
    ld_params = {"ld_wind_cm": profile["ld"]["ld_wind_cm"], "thin_annot": True,
                 "print_snps": bundle.meta["layout"]["print_snps"]}
    software = software_inventory()
    sci_key = sha256_bytes(canonical_json({
        "annotations": [[a, ann_hashes[a]] for a in ids], "bundle_id": bundle.id,
        "snp_universe": {k: bundle.meta["snp_universe"][k] for k in ("reference_order_sha256", "regression_order_sha256")},
        "profile": profile, "ldsc_commit": LDSC_COMMIT, "image": software["image"],
    }))
    resolved = {"request": request, "scientific_profile": profile, "input_sha256": input_hashes,
                "annotation_sha256": ann_hashes, "bundle_id": bundle.id, "scientific_key": sci_key,
                "ld_params": ld_params, "software": software, "resolved": now()}
    write_json_atomic(work / "resolved-request.json", resolved)

    complete = work / "COMPLETE.json"
    if complete.exists():
        rec = read_json(complete)
        out = work / "output" / "ldsR"
        if rec.get("scientific_key") == sci_key and out.is_dir() and \
                sha256_file(out / "manifest.json") == rec.get("manifest_sha256"):
            log("run already complete with identical scientific key; nothing to do")
            return 0
        log("existing COMPLETE.json is stale; removing it and the previous output")
        complete.unlink()
        shutil.rmtree(out, ignore_errors=True)

    # ---- LD tasks
    state.update(status="running", stage="ldscores")
    intervals = {a: normalized[a].intervals for a in ids}
    batches = [tuple(ids[i:i + per_task]) for i in range(0, len(ids), per_task)]
    tasks, by_ann = [], {a: {} for a in ids}
    for b in batches:
        key = b[0] if len(b) == 1 else "batch-" + sha256_bytes("\n".join(b).encode())[:12]
        for c in sorted(AUTOSOMES, key=lambda c: c):  # chr1 is largest; start big tasks first
            t = Task(key, b, c)
            tasks.append(t)
            for a in b:
                by_ann[a][c] = t
    keys = {t: input_key(t, ann_hashes, bundle, ld_params) for t in tasks}
    pending = [t for t in tasks if not is_complete(work, t, keys[t])]
    log(f"{len(tasks)} tasks, {len(tasks) - len(pending)} already complete, {len(pending)} to run with {jobs} workers")
    failures = {}
    counts = {"total": len(tasks), "completed": len(tasks) - len(pending), "failed": 0, "running": 0}
    state.update(tasks=counts)
    with ThreadPoolExecutor(max_workers=max(1, jobs)) as pool:
        queue = list(pending)
        running = {}
        while queue or running:
            while queue and len(running) < max(1, jobs):
                t = queue.pop(0)
                running[pool.submit(run_task, work, t, intervals, bundle, ld_params, keys[t])] = t
                counts["running"] = len(running)
                state.update(tasks=counts)
            done, _ = wait(list(running), timeout=30, return_when=FIRST_COMPLETED)
            for f in done:
                t = running.pop(f)
                try:
                    f.result()
                    counts["completed"] += 1
                    log(f"task {t.name} complete ({counts['completed']}/{counts['total']})")
                except WorkerError as e:
                    counts["failed"] += 1
                    failures[t.name] = e.to_json()
                    log(f"task {t.name} FAILED: {e}")
                counts["running"] = len(running)
                state.update(tasks=counts)
    if failures:
        write_json_atomic(work / "logs" / "task-failures.json", failures)
        raise WorkerError("LDSC_TASK_FAILED", f"{len(failures)} tasks failed; completed tasks are kept for resume",
                          dict(list(failures.items())[:10]))
    records = {t.name: read_json(work / "tasks" / t.key / f"chr{t.chrom}" / "COMPLETE.json") for t in tasks}

    # ---- export
    state.update(stage="export")
    out_root = work / "output"
    tmp = out_root / f".tmp-ldsR-{os.getpid()}"
    shutil.rmtree(tmp, ignore_errors=True)
    log("exporting ldsR parquet")
    counts_by_ann = export_ldsr(work, tmp, ids, by_ann, records, bundle)
    log("validating export (python + ldsR)")
    validation = validate_ldsr_dir(tmp, bundle, ids, run_r=True)
    qc = _qc(ids, normalized, counts_by_ann, records, by_ann)
    write_json_atomic(tmp / "qc.json", qc)
    pd.DataFrame(qc["annotations"]).drop(columns=["records_dropped", "zero_coverage_chromosomes"]) \
        .assign(records_dropped=[json.dumps(r["records_dropped"], sort_keys=True) for r in qc["annotations"]],
                zero_coverage_chromosomes=[",".join(map(str, r["zero_coverage_chromosomes"])) for r in qc["annotations"]]) \
        .to_csv(tmp / "qc.tsv", sep="\t", index=False)
    for a in ids:
        shutil.copyfile(work / "normalized" / f"{a}.bed", tmp / f"{a}.canonical.bed")
    manifest = {
        "schema_version": 1, "kind": "agent-ldsc/ldsR-ldscores", "created": now(),
        "scientific_key": sci_key, "request": request, "scientific_profile": profile,
        "reference": {"bundle_id": bundle.id, "genome_build": bundle.meta["genome_build"],
                      "ancestry": bundle.meta["ancestry"], "created": bundle.meta["created"],
                      "sources": {k: {f: v[f] for f in ("url", "sha256") if f in v} for k, v in bundle.meta["sources"].items()},
                      "snp_universe": bundle.meta["snp_universe"]},
        "ld_params": ld_params, "software": software,
        "annotations": [{"id": a, "input_sha256": input_hashes[a], "canonical_bed_sha256": ann_hashes[a],
                         **counts_by_ann[a], "normalization": normalized[a].summary} for a in ids],
        "validation": validation,
        "ldsr_compatible": bool(validation.get("ldsr_mask_compatible")) and validation["ldsR"]["status"] == "pass",
        "outputs": file_hashes(tmp),
    }
    write_json_atomic(tmp / "manifest.json", manifest)
    final = out_root / "ldsR"
    if final.exists():
        shutil.rmtree(final)
    os.replace(tmp, final)
    write_json_atomic(work / "COMPLETE.json", {"completed": now(), "scientific_key": sci_key,
                                              "manifest_sha256": sha256_file(final / "manifest.json"),
                                              "ldsr_compatible": manifest["ldsr_compatible"]})
    log(f"complete: {final}")
    return 0


def _write_normalized(d: Path, ann: NormalizedAnnotation) -> None:
    d.mkdir(parents=True, exist_ok=True)
    tmp = d / f".{ann.id}.bed.tmp"
    tmp.write_bytes(ann.canonical_bed())
    os.replace(tmp, d / f"{ann.id}.bed")
    ann.audit.to_csv(d / f"{ann.id}.audit.tsv.gz", sep="\t", index=False,
                     compression={"method": "gzip", "mtime": 0})
    write_json_atomic(d / f"{ann.id}.summary.json", {**ann.summary, "canonical_bed_sha256": ann.content_hash})


def _qc(ids, normalized, counts, records, by_ann) -> dict:
    rows = []
    for a in ids:
        s = normalized[a].summary
        zero = [c for c in AUTOSOMES
                if records[by_ann[a][c].name]["M"][list(by_ann[a][c].annotations).index(a)] == 0]
        rows.append({"annotation": a, "type": s["type"], "records_in": s["records_in"],
                     "records_dropped": s["records_dropped"], "merged_intervals": s["merged_intervals"],
                     "merged_bp": s["merged_bp"], **counts[a], "zero_coverage_chromosomes": zero})
    warnings = [f"{r['annotation']}: zero annotated SNPs on chromosomes {r['zero_coverage_chromosomes']}"
                for r in rows if r["zero_coverage_chromosomes"]]
    return {"annotations": rows, "warnings": warnings,
            "tasks": {name: {k: rec[k] for k in ("M", "M_5_50", "n_regression_snps", "finished")}
                      for name, rec in sorted(records.items())}}


def software_inventory() -> dict:
    def version(cmd):
        try:
            return subprocess.run(cmd, capture_output=True, text=True, timeout=60).stdout.strip()[:200]
        except (OSError, subprocess.TimeoutExpired):
            return None
    return {
        "image": os.environ.get("AGENT_LDSC_IMAGE", "unknown"),
        "image_source_sha": os.environ.get("AGENT_LDSC_SOURCE_SHA", "unknown"),
        "worker_version": __version__,
        "ldsc_commit": LDSC_COMMIT,
        "ldsr_commit": os.environ.get("AGENT_LDSC_LDSR_COMMIT", "unknown"),
        "python": platform.python_version(),
        "r": version(["Rscript", "-e", "cat(R.version.string)"]),
        "liftover": "ucsc-liftover 482 (conda lock)" if shutil.which("liftOver") else None,
        "platform": platform.platform(),
    }


def status(work: Path) -> dict:
    st = read_json(work / "state.json") if (work / "state.json").exists() else {"status": "unknown"}
    st["complete"] = (work / "COMPLETE.json").exists()
    return st
