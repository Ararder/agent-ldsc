"""Shared helpers: errors, hashing, atomic writes, directory locks, JSON I/O."""

from __future__ import annotations

import contextlib
import datetime as dt
import gzip
import hashlib
import json
import os
import shutil
import socket
import tempfile
import time
from pathlib import Path
from typing import Any, Iterator

AUTOSOMES = tuple(range(1, 23))


class WorkerError(Exception):
    """Failure with a stable machine-readable code (see docs/input-output-contract.md)."""

    CODES = {
        "INPUT_INVALID",
        "BUILD_UNSUPPORTED",
        "MAPPING_LOSS",
        "GENE_AMBIGUOUS",
        "REFERENCE_MISSING",
        "REFERENCE_MISMATCH",
        "DOWNLOAD_FAILED",
        "LDSC_TASK_FAILED",
        "RESOURCE_LIMIT",
        "EXPORT_INVALID",
        "RUN_LOCKED",
        "STATE_INVALID",
    }

    def __init__(self, code: str, message: str, details: Any = None):
        if code not in self.CODES:
            raise ValueError(f"unknown error code {code}")
        super().__init__(f"{code}: {message}")
        self.code = code
        self.message = message
        self.details = details

    def to_json(self) -> dict:
        return {"code": self.code, "message": self.message, "details": self.details}


def now() -> str:
    return dt.datetime.now(dt.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def sha256_file(path: Path | str, chunk: int = 1 << 20) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        while block := f.read(chunk):
            h.update(block)
    return h.hexdigest()


def md5_file(path: Path | str, chunk: int = 1 << 20) -> str:
    h = hashlib.md5()
    with open(path, "rb") as f:
        while block := f.read(chunk):
            h.update(block)
    return h.hexdigest()


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def sha256_strings(values) -> str:
    """Order-sensitive hash of a sequence of strings (newline-joined)."""
    h = hashlib.sha256()
    for v in values:
        h.update(str(v).encode())
        h.update(b"\n")
    return h.hexdigest()


def canonical_json(obj: Any) -> bytes:
    return json.dumps(obj, sort_keys=True, separators=(",", ":")).encode()


def read_json(path: Path | str) -> Any:
    with open(path) as f:
        return json.load(f)


def write_json_atomic(path: Path | str, obj: Any) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp = tempfile.mkstemp(dir=path.parent, prefix=f".{path.name}.", suffix=".tmp")
    try:
        with os.fdopen(fd, "w") as f:
            json.dump(obj, f, indent=2, sort_keys=True)
            f.write("\n")
            f.flush()
            os.fsync(f.fileno())
        os.replace(tmp, path)
    except BaseException:
        with contextlib.suppress(FileNotFoundError):
            os.unlink(tmp)
        raise


def open_text(path: Path | str):
    path = str(path)
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path)


def publish_dir(tmp: Path, final: Path) -> None:
    """Atomically move a fully written directory into place (same filesystem)."""
    if final.exists():
        shutil.rmtree(final)
    os.replace(tmp, final)


@contextlib.contextmanager
def dir_lock(path: Path, timeout: float | None = None, poll: float = 5.0,
             stale_after: float | None = None) -> Iterator[None]:
    """mkdir-based lock; works on NFS/Lustre where flock may not.

    The lock directory holds owner.json (host, pid, time). With stale_after set, a lock whose
    owner is on this host and no longer alive is broken; locks from other hosts are never
    broken automatically.
    """
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    start = time.monotonic()
    while True:
        try:
            path.mkdir()
            break
        except FileExistsError:
            owner = _lock_owner(path)
            if stale_after is not None and owner and _owner_dead(owner):
                shutil.rmtree(path, ignore_errors=True)
                continue
            if timeout is not None and time.monotonic() - start >= timeout:
                raise WorkerError("RUN_LOCKED", f"lock held: {path}", owner)
            time.sleep(poll)
    try:
        write_json_atomic(path / "owner.json",
                          {"host": socket.gethostname(), "pid": os.getpid(), "since": now()})
        yield
    finally:
        shutil.rmtree(path, ignore_errors=True)


def _lock_owner(path: Path) -> dict | None:
    try:
        return read_json(path / "owner.json")
    except (FileNotFoundError, json.JSONDecodeError):
        return None


def _owner_dead(owner: dict) -> bool:
    if owner.get("host") != socket.gethostname():
        return False
    try:
        os.kill(int(owner["pid"]), 0)
    except ProcessLookupError:
        return True
    except (PermissionError, KeyError, ValueError):
        return False
    return False


def normalize_chrom(value: str) -> int | None:
    """Recognized autosome aliases ('chr7', '7') -> 7; anything else -> None."""
    v = value[3:] if value.lower().startswith("chr") else value
    if v.isdigit() and str(int(v)) == v and 1 <= int(v) <= 22:
        return int(v)
    return None
