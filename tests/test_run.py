import json
import shutil
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow.parquet as pq
import pytest

from agent_ldsc_worker.run import run
from agent_ldsc_worker.util import read_json
from oracle import ld_scores, read_plink_bed

pytestmark = [pytest.mark.ldsc, pytest.mark.rgate]


def make_run(tmp_path, annotations, files):
    work = tmp_path / "run"
    (work / "input").mkdir(parents=True)
    for name, text in files.items():
        (work / "input" / name).write_text(text)
    req = {"schema_version": 1, "reference_id": "fixture-grch37-v1", "scientific_profile": "strict-v1",
           "annotations": annotations}
    (work / "input" / "request.json").write_text(json.dumps(req))
    return work


# chr1 cluster 0 spans BP 100001..198001 (step 2000), chr2 cluster 1 spans 900001..998001
BED_A = "chr1\t100000\t150000\nchr1\t140000\t160000\nchr2\t900000\t901000\n"   # overlapping peaks merge
BED_B = "chr6\t100000\t300000\nchr1\t1700000\t1750000\n"                       # chr6 only + chr1 cluster 2


def oracle_for(bundle_dir: Path, annots: dict) -> dict:
    ref = pq.read_table(bundle_dir / "snps/reference_snps.parquet").to_pandas()
    reg = pq.read_table(bundle_dir / "snps/regression_snps.parquet").to_pandas()
    out = {}
    for name, merged in annots.items():
        cols = []
        for c in range(1, 23):
            r = ref[ref.CHR == c]
            g = read_plink_bed(bundle_dir / f"plink/fixture.{c}.bed", 61, len(r))
            a = np.zeros(len(r))
            for s, e in merged.get(c, []):
                a += ((r.BP.to_numpy() - 1 >= s) & (r.BP.to_numpy() - 1 < e))
            cluster = r.SNP.str.split("_").str[1].astype(int).to_numpy()
            l2 = ld_scores(g, cluster, a[:, None])[:, 0]
            cols.append(pd.Series(l2, index=r.SNP.to_numpy()))
        out[name] = pd.concat(cols).loc[reg.SNP].to_numpy()
    return out


def test_end_to_end_matches_independent_oracle(tmp_path, fixture_bundle):
    anns = [{"id": "peaksA", "type": "bed", "file": "a.bed", "species": "human", "build": "GRCh37"},
            {"id": "peaksB", "type": "bed", "file": "b.bed", "species": "human", "build": "GRCh37"}]
    work = make_run(tmp_path, anns, {"a.bed": BED_A, "b.bed": BED_B})
    assert run(work / "input/request.json", fixture_bundle, work, jobs=4) == 0

    out = work / "output/ldsR"
    ld = pq.read_table(out / "ld.parquet").to_pandas()
    annot = pq.read_table(out / "annot.parquet").to_pandas()
    ref = pq.read_table(out / "annot_ref.parquet").to_pandas()
    assert list(ld.columns) == ["SNP", "peaksA", "peaksB"]
    assert annot.annot.tolist() == ["peaksA", "peaksB"]
    assert str(ref.dtypes["peaksA"]) == "int32"

    exp = oracle_for(fixture_bundle, {"peaksA": {1: [(100000, 160000)], 2: [(900000, 901000)]},
                                      "peaksB": {6: [(100000, 300000)], 1: [(1700000, 1750000)]}})
    for a in ("peaksA", "peaksB"):
        err = np.abs(ld[a].to_numpy() - exp[a])
        assert err.max() <= 5e-4 + 1e-9, (a, err.max())   # ldsc prints %.3f

    # exact discrete counts: chr1 cluster0 SNPs with BP-1 in [100000,160000): BP 100001..159001 -> 30 SNPs
    bundle_ref = pq.read_table(fixture_bundle / "snps/reference_snps.parquet").to_pandas()
    m_a = int(((bundle_ref.CHR == 1) & (bundle_ref.BP - 1 >= 100000) & (bundle_ref.BP - 1 < 160000)).sum() +
              ((bundle_ref.CHR == 2) & (bundle_ref.BP - 1 >= 900000) & (bundle_ref.BP - 1 < 901000)).sum())
    assert annot.m[0] == m_a == ref.peaksA.sum() == 31
    common = bundle_ref.common.to_numpy()
    assert annot.m50.tolist() == [ref.peaksA[common].sum(), ref.peaksB[common].sum()]

    qc = read_json(out / "qc.json")
    zero_b = qc["annotations"][1]["zero_coverage_chromosomes"]
    assert 2 in zero_b and 6 not in zero_b                  # zero coverage on a chromosome is valid
    manifest = read_json(out / "manifest.json")
    assert manifest["ldsr_compatible"] is False             # fixture is not ldsR-mask compatible
    assert manifest["validation"]["ldsR"]["status"] == "pass"
    assert read_json(work / "COMPLETE.json")["manifest_sha256"]
    assert read_json(work / "state.json")["status"] == "completed"


def test_batched_tasks_equal_single_annotation_tasks(tmp_path, fixture_bundle):
    anns = [{"id": "peaksA", "type": "bed", "file": "a.bed", "species": "human", "build": "GRCh37"},
            {"id": "peaksB", "type": "bed", "file": "b.bed", "species": "human", "build": "GRCh37"}]
    w1 = make_run(tmp_path / "one", anns, {"a.bed": BED_A, "b.bed": BED_B})
    w2 = make_run(tmp_path / "two", anns, {"a.bed": BED_A, "b.bed": BED_B})
    assert run(w1 / "input/request.json", fixture_bundle, w1, jobs=4, annotations_per_task=1) == 0
    assert run(w2 / "input/request.json", fixture_bundle, w2, jobs=4, annotations_per_task=2) == 0
    for f in ("ld.parquet", "annot.parquet", "annot_ref.parquet"):
        a = pq.read_table(w1 / "output/ldsR" / f).to_pandas()
        b = pq.read_table(w2 / "output/ldsR" / f).to_pandas()
        pd.testing.assert_frame_equal(a, b)


def test_resume_skips_completed_tasks_and_redoes_tampered(tmp_path, fixture_bundle):
    anns = [{"id": "peaksA", "type": "bed", "file": "a.bed", "species": "human", "build": "GRCh37"}]
    work = make_run(tmp_path, anns, {"a.bed": BED_A})
    assert run(work / "input/request.json", fixture_bundle, work, jobs=4) == 0
    first = read_json(work / "tasks/peaksA/chr5/COMPLETE.json")["finished"]
    # simulate interruption: drop the final output + one task, corrupt another
    (work / "COMPLETE.json").unlink()
    shutil.rmtree(work / "output/ldsR")
    shutil.rmtree(work / "tasks/peaksA/chr3")
    with open(work / "tasks/peaksA/chr4/ld.4.l2.M", "a") as f:
        f.write("tampered\n")
    assert run(work / "input/request.json", fixture_bundle, work, jobs=2) == 0
    assert read_json(work / "tasks/peaksA/chr5/COMPLETE.json")["finished"] == first
    assert read_json(work / "tasks/peaksA/chr3/COMPLETE.json")["finished"] != first
    log = (work / "logs/worker.log").read_text()
    assert "22 tasks, 20 already complete, 2 to run" in log


def test_changed_input_invalidates_stale_output(tmp_path, fixture_bundle):
    anns = [{"id": "peaksA", "type": "bed", "file": "a.bed", "species": "human", "build": "GRCh37"}]
    work = make_run(tmp_path, anns, {"a.bed": BED_A})
    assert run(work / "input/request.json", fixture_bundle, work, jobs=4) == 0
    old_key = read_json(work / "COMPLETE.json")["scientific_key"]
    (work / "input/a.bed").write_text("chr1\t100000\t120000\n")
    assert run(work / "input/request.json", fixture_bundle, work, jobs=4) == 0
    assert read_json(work / "COMPLETE.json")["scientific_key"] != old_key
    assert pq.read_table(work / "output/ldsR/annot.parquet").to_pandas().m[0] == 10


def test_genome_wide_zero_coverage_fails_without_output(tmp_path, fixture_bundle):
    anns = [{"id": "empty", "type": "bed", "file": "a.bed", "species": "human", "build": "GRCh37"}]
    work = make_run(tmp_path, anns, {"a.bed": "chr1\t2500000\t2600000\n"})   # between clusters
    assert run(work / "input/request.json", fixture_bundle, work, jobs=4) == 2
    assert not (work / "COMPLETE.json").exists()
    assert not (work / "output/ldsR").exists()
    assert read_json(work / "state.json")["status"] == "failed"


def test_failed_task_returns_nonzero_and_publishes_nothing(tmp_path, fixture_bundle, monkeypatch):
    anns = [{"id": "peaksA", "type": "bed", "file": "a.bed", "species": "human", "build": "GRCh37"}]
    work = make_run(tmp_path, anns, {"a.bed": BED_A})
    monkeypatch.setattr("agent_ldsc_worker.ldscores.LDSC_CMD", ["false"])
    assert run(work / "input/request.json", fixture_bundle, work, jobs=4) == 2
    assert read_json(work / "logs/error.json")["code"] == "LDSC_TASK_FAILED"
    assert not (work / "COMPLETE.json").exists()


def test_corrupted_reference_is_rejected(tmp_path, fixture_bundle):
    bad = tmp_path / "bundle"
    shutil.copytree(fixture_bundle, bad)
    with open(bad / "plink/fixture.7.bim", "a") as f:
        f.write("7\trs_extra\t9\t2900000\tA\tG\n")
    anns = [{"id": "peaksA", "type": "bed", "file": "a.bed", "species": "human", "build": "GRCh37"}]
    work = make_run(tmp_path, anns, {"a.bed": BED_A})
    assert run(work / "input/request.json", bad, work) == 2
    assert read_json(work / "logs/error.json")["code"] == "REFERENCE_MISMATCH"


def test_concurrent_run_is_locked(tmp_path, fixture_bundle):
    anns = [{"id": "peaksA", "type": "bed", "file": "a.bed", "species": "human", "build": "GRCh37"}]
    work = make_run(tmp_path, anns, {"a.bed": BED_A})
    (work / ".run.lock").mkdir()
    (work / ".run.lock/owner.json").write_text(json.dumps({"host": "elsewhere", "pid": 1}))
    assert run(work / "input/request.json", fixture_bundle, work) == 3
    assert run(work / "input/request.json", fixture_bundle, work, force_unlock=True, jobs=4) == 0
