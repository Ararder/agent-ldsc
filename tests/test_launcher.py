import json
from pathlib import Path

import pytest

from agent_ldsc_launcher.cli import LauncherError, _check_profile, job_script, prepare_inputs

ROOT = Path(__file__).resolve().parents[1]
DIGEST = "ghcr.io/ararder/agent-ldsc@sha256:" + "a" * 64


def profile(name):
    p = json.loads((ROOT / f"profiles/execution/{name}.example.json").read_text())
    p["image"] = DIGEST
    return p


def test_profiles_require_digest_pinned_image():
    p = profile("dardel")
    _check_profile(p)
    p["image"] = "ghcr.io/ararder/agent-ldsc:latest"
    with pytest.raises(LauncherError):
        _check_profile(p)


def test_slurm_job_script_runs_worker_offline_with_readonly_refs():
    p = profile("dardel")
    s = job_script(p, "/w/r1", "/c/references/b", "r1", force_unlock=False)
    assert "#SBATCH --account=naiss2026-4-1187" in s and "#SBATCH --partition=shared" in s
    assert "--containall" in s and "/c/references/b:/refs:ro" in s and "/w/r1:/work:rw" in s
    assert "--jobs ${SLURM_CPUS_PER_TASK}" in s and "--force-unlock" not in s
    assert s.index("ml PDCOLD") < s.index("apptainer exec")
    assert "--force-unlock" in job_script(p, "/w/r1", "/c/references/b", "r1", force_unlock=True)


def test_lab_profile_has_no_account_or_modules():
    s = job_script(profile("lab"), "/w/r1", "/c/b", "r1", force_unlock=False)
    assert "--account" not in s and "ml " not in s and "#SBATCH --partition=main" in s


def test_docker_direct_runs_without_network():
    s = job_script(profile("local-docker"), "/w/r1", "/c/b", "r1", force_unlock=False)
    assert "docker run --rm --network none" in s and "#SBATCH" not in s


def test_prepare_inputs_sanitizes_names_and_rewrites_paths(tmp_path):
    src = tmp_path / "client"
    src.mkdir()
    (src / "my peaks;rm -rf.bed").write_text("chr1\t1\t2\n")
    req = {"schema_version": 1, "reference_id": "x", "scientific_profile": "strict-v1",
           "annotations": [{"id": "a", "type": "bed", "file": "my peaks;rm -rf.bed", "species": "human",
                            "build": "GRCh37"}]}
    (src / "req.json").write_text(json.dumps(req))
    out = tmp_path / "staged"
    staged, h1 = prepare_inputs(src / "req.json", out)
    assert staged["annotations"][0]["file"] == "files/00-my_peaks_rm_-rf.bed"
    assert (out / staged["annotations"][0]["file"]).read_text() == "chr1\t1\t2\n"
    _, h2 = prepare_inputs(src / "req.json", tmp_path / "staged2")
    assert h1 == h2
    (src / "my peaks;rm -rf.bed").write_text("chr1\t1\t3\n")
    _, h3 = prepare_inputs(src / "req.json", tmp_path / "staged3")
    assert h3 != h1


def test_missing_annotation_file_fails(tmp_path):
    (tmp_path / "req.json").write_text(json.dumps({"annotations": [{"id": "a", "file": "nope.bed"}]}))
    with pytest.raises(LauncherError):
        prepare_inputs(tmp_path / "req.json", tmp_path / "out")
