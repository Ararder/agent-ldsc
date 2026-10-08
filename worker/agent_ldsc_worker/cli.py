"""agent-ldsc-worker command line (runs inside the image)."""

from __future__ import annotations

import argparse
import json
import os
import sys
from pathlib import Path

from .util import WorkerError


def _default_jobs() -> int:
    return int(os.environ.get("SLURM_CPUS_PER_TASK") or 1)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(prog="agent-ldsc-worker")
    sub = ap.add_subparsers(dest="cmd", required=True)

    refs = sub.add_parser("refs", help="reference bundles").add_subparsers(dest="refs_cmd", required=True)
    ri = refs.add_parser("install", help="download, verify, derive and publish a bundle")
    ri.add_argument("--bundle", required=True)
    ri.add_argument("--cache", required=True, type=Path)
    ri.add_argument("--source-dir", action="append", type=Path, default=[],
                    help="directory with pre-downloaded source files (verified by checksum)")
    ri.add_argument("--offline", action="store_true")
    ri.add_argument("--ldsr-mask", type=Path)
    rv = refs.add_parser("verify", help="re-hash every file of an installed bundle")
    rv.add_argument("bundle_dir", type=Path)

    r = sub.add_parser("run", help="execute or resume a run")
    r.add_argument("--request", required=True, type=Path)
    r.add_argument("--refs", required=True, type=Path, help="installed bundle directory")
    r.add_argument("--work", required=True, type=Path, help="run directory")
    r.add_argument("--jobs", type=int, default=_default_jobs())
    r.add_argument("--annotations-per-task", type=int, default=1)
    r.add_argument("--force-unlock", action="store_true",
                   help="remove a stale run lock (launcher uses this only after checking the scheduler)")

    v = sub.add_parser("validate-output", help="re-validate an exported ldsR directory")
    v.add_argument("dir", type=Path)
    v.add_argument("--refs", required=True, type=Path)

    s = sub.add_parser("status", help="print run state")
    s.add_argument("--work", required=True, type=Path)

    sub.add_parser("versions", help="print software inventory")

    a = ap.parse_args(argv)
    try:
        if a.cmd == "refs" and a.refs_cmd == "install":
            from .refs import install_bundle
            print(install_bundle(a.bundle, a.cache, source_dirs=a.source_dir, offline=a.offline,
                                 ldsr_mask=a.ldsr_mask))
            return 0
        if a.cmd == "refs" and a.refs_cmd == "verify":
            from .refs import load_bundle
            b = load_bundle(a.bundle_dir, verify=True)
            print(json.dumps({"bundle_id": b.id, "files": len(b.meta["files"]), "status": "ok"}))
            return 0
        if a.cmd == "run":
            from .run import run
            return run(a.request, a.refs, a.work, a.jobs, a.annotations_per_task, a.force_unlock)
        if a.cmd == "validate-output":
            from .export import validate_ldsr_dir
            from .refs import load_bundle
            print(json.dumps(validate_ldsr_dir(a.dir, load_bundle(a.refs, verify=False)), default=str))
            return 0
        if a.cmd == "status":
            from .run import status
            print(json.dumps(status(a.work), indent=2))
            return 0
        if a.cmd == "versions":
            from .run import software_inventory
            print(json.dumps(software_inventory(), indent=2))
            return 0
    except WorkerError as e:
        print(json.dumps({"error": e.to_json()}, default=str), file=sys.stderr)
        return 2
    return 1
