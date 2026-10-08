"""agent-ldsc host launcher: staging, transport and scheduling only (no science, stdlib only).

Profiles are JSON (see profiles/execution/*.example.json). Lookup order for --profile NAME:
an explicit path, $AGENT_LDSC_PROFILE_DIR/NAME.json, ~/.config/agent-ldsc/profiles/NAME.json.
Local run records live in ~/.local/state/agent-ldsc/runs/RUN_ID.json.
"""

from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import json
import os
import re
import shlex
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

STATE_DIR = Path(os.environ.get("AGENT_LDSC_STATE_DIR", Path.home() / ".local/state/agent-ldsc/runs"))
PROFILE_DIRS = [Path(p) for p in filter(None, [os.environ.get("AGENT_LDSC_PROFILE_DIR")])] + \
    [Path.home() / ".config/agent-ldsc/profiles"]
DIGEST_RE = re.compile(r"^[\w.\-/:]+@sha256:[0-9a-f]{64}$")


class LauncherError(Exception):
    pass


# --------------------------------------------------------------------------- profiles


def load_profile(name: str) -> dict:
    cands = [Path(name)] if name.endswith(".json") else [d / f"{name}.json" for d in PROFILE_DIRS]
    for p in cands:
        if p.is_file():
            prof = json.loads(p.read_text())
            prof.setdefault("profile_id", p.stem)
            _check_profile(prof)
            return prof
    raise LauncherError(f"profile {name!r} not found in {[str(c) for c in cands]}")


def _check_profile(p: dict) -> None:
    for key in ("transport", "executor", "runtime", "image", "cache_root", "work_root"):
        if key not in p:
            raise LauncherError(f"profile {p.get('profile_id')}: missing {key}")
    if not DIGEST_RE.match(p["image"]):
        raise LauncherError("profile image must be pinned by digest (name@sha256:...), never a floating tag")
    if p["executor"]["type"] not in ("slurm", "direct"):
        raise LauncherError("executor.type must be slurm or direct")
    if p["runtime"]["type"] not in ("apptainer", "docker"):
        raise LauncherError("runtime.type must be apptainer or docker")
    if p["transport"]["type"] not in ("ssh", "local"):
        raise LauncherError("transport.type must be ssh or local")
    for root in ("cache_root", "work_root"):
        if not p[root].startswith("/"):
            raise LauncherError(f"{root} must be an absolute path")


# --------------------------------------------------------------------------- transport


def sh(profile: dict, script: str, check: bool = True, quiet: bool = False) -> subprocess.CompletedProcess:
    """Run a bash login-shell script on the execution site (stdin-fed, so the script itself needs
    no remote quoting; login shell so site module systems such as Lmod are initialized).
    stdout is captured; stderr streams to the terminal unless quiet."""
    if profile["transport"]["type"] == "ssh":
        cmd = ["ssh", "-o", "BatchMode=yes", profile["transport"]["host"], "bash", "-l", "-s"]
    else:
        cmd = ["bash", "-l", "-s"]
    proc = subprocess.run(cmd, input="set -euo pipefail\n" + script, text=True, stdout=subprocess.PIPE,
                          stderr=subprocess.PIPE if quiet else None)
    if check and proc.returncode != 0:
        raise LauncherError(f"site command failed with exit {proc.returncode}" +
                            (f":\n{proc.stderr[-3000:]}" if quiet and proc.stderr else " (see stderr above)"))
    return proc


def push(profile: dict, local_dir: Path, remote_dir: str) -> None:
    if profile["transport"]["type"] == "ssh":
        dest = f"{profile['transport']['host']}:{remote_dir}/"
        subprocess.run(["rsync", "-a", "--chmod=Fu=rw,Fgo=r,Du=rwx,Dgo=rx", f"{local_dir}/", dest], check=True)
    else:
        shutil.copytree(local_dir, remote_dir, dirs_exist_ok=True)


def pull(profile: dict, remote_dir: str, local_dir: Path) -> None:
    local_dir.mkdir(parents=True, exist_ok=True)
    if profile["transport"]["type"] == "ssh":
        subprocess.run(["rsync", "-a", f"{profile['transport']['host']}:{remote_dir}/", f"{local_dir}/"], check=True)
    else:
        shutil.copytree(remote_dir, local_dir, dirs_exist_ok=True)


# --------------------------------------------------------------------------- runtime


def _q(s) -> str:
    return shlex.quote(str(s))


def runtime_prelude(profile: dict) -> str:
    rt = profile["runtime"]
    lines = [rt.get("setup", "")]
    if rt["type"] == "apptainer":
        cache = profile["cache_root"]
        lines += [f"export APPTAINER_CACHEDIR={_q(cache + '/apptainer-cache')}",
                  f"export APPTAINER_TMPDIR={_q(cache + '/apptainer-tmp')}",
                  'mkdir -p "$APPTAINER_CACHEDIR" "$APPTAINER_TMPDIR"']
    return "\n".join(lines) + "\n"


def image_ref(profile: dict) -> str:
    """Path of the runnable image on the site (SIF for apptainer, digest ref for docker)."""
    if profile["runtime"]["type"] == "apptainer":
        digest = profile["image"].split("@sha256:")[1]
        return f"{profile['cache_root']}/images/{digest}.sif"
    return profile["image"]


def container_exec(profile: dict, binds: list[tuple[str, str, str]], args: list[str], env: dict | None = None) -> str:
    """Shell command running `args` in the image. binds: (host, container, 'ro'|'rw')."""
    env = env or {}
    if profile["runtime"]["type"] == "apptainer":
        b = ",".join(f"{h}:{c}:{m}" for h, c, m in binds)
        envs = " ".join(f"--env {_q(f'{k}={v}')}" for k, v in env.items())
        return f"apptainer exec --cleanenv --containall --no-home {envs} --bind {_q(b)} {_q(image_ref(profile))} " + \
            " ".join(map(_q, args))
    b = " ".join(f"-v {_q(f'{h}:{c}:{m}')}" for h, c, m in binds)
    envs = " ".join(f"-e {_q(f'{k}={v}')}" for k, v in env.items())
    return f"docker run --rm --network none --user $(id -u):$(id -g) {envs} {b} {_q(profile['image'])} " + \
        " ".join(map(_q, args[1:] if args[0] == "agent-ldsc-worker" else args))


def stage(profile: dict, bundle_id: str) -> dict:
    """Ensure the image and reference bundle exist on the site (network-enabled host)."""
    cache = profile["cache_root"]
    rt = profile["runtime"]["type"]
    sif = image_ref(profile)
    script = runtime_prelude(profile) + f"mkdir -p {_q(cache)}/images {_q(cache)}/references\n"
    if rt == "apptainer":
        script += f"""
if [ ! -s {_q(sif)} ]; then
  tmp={_q(sif)}.tmp.$$
  apptainer pull --disable-cache "$tmp" {_q('docker://' + profile['image'])}
  mv "$tmp" {_q(sif)}
  sha256sum {_q(sif)} | cut -d' ' -f1 > {_q(sif)}.sha256
fi
test "$(sha256sum {_q(sif)} | cut -d' ' -f1)" = "$(cat {_q(sif)}.sha256)" || {{ echo "SIF checksum mismatch" >&2; exit 1; }}
"""
    else:
        script += f"docker pull {_q(profile['image'])} >/dev/null\n"
    install = container_exec(profile, [(cache, "/cache", "rw")],
                             ["agent-ldsc-worker", "refs", "install", "--bundle", bundle_id, "--cache", "/cache"])
    if rt == "docker":
        install = install.replace("--network none ", "")  # staging is the only networked step
    script += f"""
if [ ! -s {_q(cache)}/references/{_q(bundle_id)}/BUNDLE.json ]; then
  {install}
fi
echo "SIF_SHA256=$(cat {_q(sif)}.sha256 2>/dev/null || echo n/a)"
"""
    out = sh(profile, script).stdout
    sif_sha = re.search(r"SIF_SHA256=(\S+)", out)
    return {"image": profile["image"], "sif": sif if rt == "apptainer" else None,
            "sif_sha256": sif_sha.group(1) if sif_sha else None,
            "bundle_dir": f"{cache}/references/{bundle_id}"}


# --------------------------------------------------------------------------- runs


def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def prepare_inputs(request_path: Path, dest: Path) -> tuple[dict, str]:
    """Copy request + annotation files into dest/ with sanitized names; return staged request, hash."""
    req = json.loads(request_path.read_text())
    if not isinstance(req.get("annotations"), list) or not req["annotations"]:
        raise LauncherError("request has no annotations")
    h = hashlib.sha256()
    (dest / "files").mkdir(parents=True)
    for i, a in enumerate(req["annotations"]):
        src = (request_path.parent / a["file"]).resolve()
        if not src.is_file():
            raise LauncherError(f"annotation {a.get('id')}: file not found: {src}")
        safe = re.sub(r"[^A-Za-z0-9._-]", "_", src.name)
        rel = f"files/{i:02d}-{safe}"
        shutil.copyfile(src, dest / rel)
        a["file"] = rel
        h.update(a.get("id", "").encode() + _sha256(src).encode())
    text = json.dumps(req, indent=2, sort_keys=True)
    (dest / "request.json").write_text(text + "\n")
    h.update(text.encode())
    return req, h.hexdigest()


def job_script(profile: dict, run_dir: str, bundle_dir: str, run_id: str, force_unlock: bool) -> str:
    ex = profile["executor"]
    lines = ["#!/bin/bash -l"]
    if ex["type"] == "slurm":
        lines += [f"#SBATCH --job-name=agent-ldsc-{run_id}", f"#SBATCH --output={run_dir}/logs/slurm-%j.out",
                  f"#SBATCH --cpus-per-task={ex['cpus']}", f"#SBATCH --mem={ex['mem']}", "#SBATCH --nodes=1",
                  "#SBATCH --ntasks=1"]
        for k in ("account", "partition", "time"):
            if ex.get(k):
                lines.append(f"#SBATCH --{k}={ex[k]}")
        lines += [f"#SBATCH {x}" for x in ex.get("extra_sbatch", [])]
    jobs = "${SLURM_CPUS_PER_TASK}" if ex["type"] == "slurm" else str(ex.get("cpus", 1))
    args = ["agent-ldsc-worker", "run", "--request", "/work/input/request.json", "--refs", "/refs",
            "--work", "/work", "--jobs", "JOBS"]
    if ex.get("annotations_per_task"):
        args += ["--annotations-per-task", str(int(ex["annotations_per_task"]))]
    if force_unlock:
        args.append("--force-unlock")
    env = {"AGENT_LDSC_IMAGE": profile["image"], "SLURM_JOB_ID": "${SLURM_JOB_ID:-}"}
    cmd = container_exec(profile, [(run_dir, "/work", "rw"), (bundle_dir, "/refs", "ro")], args, env)
    cmd = cmd.replace("JOBS", jobs).replace(_q("SLURM_JOB_ID=${SLURM_JOB_ID:-}"), '"SLURM_JOB_ID=${SLURM_JOB_ID:-}"')
    lines += ["set -euo pipefail", runtime_prelude(profile).strip(), f'echo "agent-ldsc {run_id} on $(hostname) at $(date -u +%FT%TZ)"', cmd]
    return "\n".join(l for l in lines if l) + "\n"


def record_path(run_id: str) -> Path:
    return STATE_DIR / f"{run_id}.json"


def save_record(rec: dict) -> None:
    STATE_DIR.mkdir(parents=True, exist_ok=True)
    tmp = record_path(rec["run_id"]).with_suffix(".tmp")
    tmp.write_text(json.dumps(rec, indent=2, sort_keys=True) + "\n")
    os.replace(tmp, record_path(rec["run_id"]))


def load_record(run_id: str) -> dict:
    p = record_path(run_id)
    if not p.exists():
        raise LauncherError(f"no local record for run {run_id} ({p})")
    return json.loads(p.read_text())


def active_jobs(profile: dict, run_id: str) -> list[str]:
    if profile["executor"]["type"] != "slurm":
        return []
    out = sh(profile, f"squeue -h --me --name={_q('agent-ldsc-' + run_id)} -o '%i %T'").stdout
    return [l.split()[0] for l in out.splitlines() if l.strip()]


def submit(profile: dict, run_dir: str, bundle_dir: str, run_id: str, force_unlock: bool) -> dict:
    script = job_script(profile, run_dir, bundle_dir, run_id, force_unlock)
    attempt = dt.datetime.now(dt.timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    job_file = f"{run_dir}/job-{attempt}.sh"
    write = f"mkdir -p {_q(run_dir)}/logs\ncat > {_q(job_file)} <<'AGENT_LDSC_EOF'\n{script}AGENT_LDSC_EOF\nchmod +x {_q(job_file)}\n"
    if profile["executor"]["type"] == "slurm":
        out = sh(profile, write + f"sbatch --parsable {_q(job_file)}").stdout.strip()
        job_id = out.split(";")[0].splitlines()[-1]
        return {"attempt": attempt, "job_file": job_file, "executor": "slurm", "job_id": job_id}
    out = sh(profile, write + f"nohup bash {_q(job_file)} > {_q(run_dir)}/logs/direct-{attempt}.out 2>&1 &\necho $!").stdout
    return {"attempt": attempt, "job_file": job_file, "executor": "direct", "pid": out.strip().splitlines()[-1]}


# --------------------------------------------------------------------------- commands


def cmd_stage(a) -> int:
    prof = load_profile(a.profile)
    print(json.dumps(stage(prof, a.bundle), indent=2))
    return 0


def cmd_submit(a) -> int:
    prof = load_profile(a.profile)
    req = json.loads(Path(a.request).read_text())
    with tempfile.TemporaryDirectory() as td:
        staged, h = prepare_inputs(Path(a.request), Path(td))
        run_id = a.run_id or dt.datetime.now().strftime("%Y%m%d-%H%M%S-") + h[:8]
        run_dir = f"{prof['work_root']}/{run_id}"
        if sh(prof, f"test -e {_q(run_dir)} && echo exists || true").stdout.strip() == "exists":
            raise LauncherError(f"run directory already exists: {run_dir} (use resume)")
        site = stage(prof, req["reference_id"])
        sh(prof, f"mkdir -p {_q(run_dir)}/input {_q(run_dir)}/logs")
        push(prof, Path(td), f"{run_dir}/input")
    rec = {"run_id": run_id, "profile": prof, "run_dir": run_dir, "request_sha256": h, "site": site,
           "submissions": [], "created": dt.datetime.now(dt.timezone.utc).isoformat()}
    save_record(rec)  # before submission, so an ambiguous failure can be reconciled
    if a.dry_run:
        print(job_script(prof, run_dir, site["bundle_dir"], run_id, False))
        return 0
    sub = submit(prof, run_dir, site["bundle_dir"], run_id, force_unlock=False)
    rec["submissions"].append(sub)
    save_record(rec)
    print(json.dumps({"run_id": run_id, "run_dir": run_dir, **sub}, indent=2))
    return 0


def cmd_resume(a) -> int:
    rec = load_record(a.run_id)
    prof = load_profile(a.profile) if a.profile else rec["profile"]
    running = active_jobs(prof, a.run_id)
    if running:
        raise LauncherError(f"run {a.run_id} still has active jobs {running}; not resubmitting")
    if rec["profile"]["image"] != prof["image"]:
        print(f"note: image changed {rec['profile']['image']} -> {prof['image']}; incompatible tasks will be recomputed",
              file=sys.stderr)
    site = stage(prof, json.loads(sh(prof, f"cat {_q(rec['run_dir'])}/input/request.json").stdout)["reference_id"])
    sub = submit(prof, rec["run_dir"], site["bundle_dir"], a.run_id, force_unlock=True)
    rec["submissions"].append(sub)
    rec["profile"] = prof
    save_record(rec)
    print(json.dumps({"run_id": a.run_id, **sub}, indent=2))
    return 0


def cmd_status(a) -> int:
    rec = load_record(a.run_id)
    prof = rec["profile"]
    rd = _q(rec["run_dir"])
    out = sh(prof, f"cat {rd}/state.json 2>/dev/null || echo '{{}}'; echo; echo '@@'; "
                   f"cat {rd}/COMPLETE.json 2>/dev/null || echo null; echo '@@'; "
                   f"cat {rd}/logs/error.json 2>/dev/null || echo null").stdout.split("@@")
    st = {"run_id": a.run_id, "run_dir": rec["run_dir"], "state": json.loads(out[0]),
          "complete": json.loads(out[1]), "error": json.loads(out[2]), "submissions": rec["submissions"]}
    if prof["executor"]["type"] == "slurm" and rec["submissions"]:
        ids = ",".join(s["job_id"] for s in rec["submissions"] if s.get("job_id"))
        st["scheduler"] = sh(prof, f"sacct -n -P -j {_q(ids)} -X -o JobID,State,Elapsed,MaxRSS,ExitCode 2>/dev/null || true",
                             check=False, quiet=True).stdout.strip().splitlines()
    if a.json:
        print(json.dumps(st, indent=2))
    else:
        s = st["state"]
        print(f"{a.run_id}: {s.get('status', 'unknown')} (stage {s.get('stage')}) tasks {s.get('tasks')}")
        if st["error"]:
            print(f"error: {st['error']['code']}: {st['error']['message']}")
        for line in st.get("scheduler", []):
            print("  slurm:", line)
    return 0


def cmd_fetch(a) -> int:
    rec = load_record(a.run_id)
    prof = rec["profile"]
    complete = sh(prof, f"cat {_q(rec['run_dir'])}/COMPLETE.json 2>/dev/null || true").stdout.strip()
    if not complete:
        raise LauncherError(f"run {a.run_id} is not complete; nothing final to fetch")
    dest = Path(a.destination) / a.run_id
    pull(prof, f"{rec['run_dir']}/output/ldsR", dest / "ldsR")
    pull(prof, f"{rec['run_dir']}/normalized", dest / "normalized")
    (dest / "COMPLETE.json").write_text(complete + "\n")
    manifest_sha = json.loads(complete)["manifest_sha256"]
    if _sha256(dest / "ldsR/manifest.json") != manifest_sha:
        raise LauncherError("fetched manifest does not match COMPLETE.json")
    manifest = json.loads((dest / "ldsR/manifest.json").read_text())
    bad = [f for f, h in manifest["outputs"].items() if _sha256(dest / "ldsR" / f) != h]
    if bad:
        raise LauncherError(f"fetched files fail checksum: {bad}")
    print(dest)
    return 0


def cmd_cancel(a) -> int:
    rec = load_record(a.run_id)
    prof = rec["profile"]
    for j in active_jobs(prof, a.run_id):
        sh(prof, f"scancel --signal=TERM --full {_q(j)}")
        print(f"cancelled {j}")
    return 0


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(prog="agent-ldsc")
    sub = ap.add_subparsers(dest="cmd", required=True)
    s = sub.add_parser("stage", help="ensure image + reference bundle on the execution site")
    s.add_argument("--profile", required=True)
    s.add_argument("--bundle", default="eur-phase3-grch37-v1")
    s.set_defaults(fn=cmd_stage)
    s = sub.add_parser("submit", help="stage inputs and submit a new run")
    s.add_argument("request")
    s.add_argument("--profile", required=True)
    s.add_argument("--run-id")
    s.add_argument("--dry-run", action="store_true", help="stage and print the job script without submitting")
    s.set_defaults(fn=cmd_submit)
    s = sub.add_parser("resume", help="resubmit an incomplete run (completed tasks are reused)")
    s.add_argument("run_id")
    s.add_argument("--profile")
    s.set_defaults(fn=cmd_resume)
    s = sub.add_parser("status")
    s.add_argument("run_id")
    s.add_argument("--json", action="store_true")
    s.set_defaults(fn=cmd_status)
    s = sub.add_parser("fetch", help="copy a completed, checksum-verified ldsR result")
    s.add_argument("run_id")
    s.add_argument("--destination", default="results")
    s.set_defaults(fn=cmd_fetch)
    s = sub.add_parser("cancel")
    s.add_argument("run_id")
    s.set_defaults(fn=cmd_cancel)
    a = ap.parse_args(argv)
    try:
        return a.fn(a)
    except LauncherError as e:
        print(f"agent-ldsc: {e}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
