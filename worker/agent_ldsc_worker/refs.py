"""Reference bundle: download, verify, derive, validate, publish atomically, and load.

A bundle directory is valid only when BUNDLE.json exists (written last) and every file it
lists has the recorded size and sha256. Builders never modify a published bundle.
"""

from __future__ import annotations

import gzip
import os
import re
import shutil
import tarfile
import time
import urllib.error
import urllib.request
import zipfile
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq

from . import __version__
from .util import (AUTOSOMES, WorkerError, dir_lock, md5_file, now, read_json, sha256_file,
                   sha256_strings, write_json_atomic)

HOME = Path(os.environ.get("AGENT_LDSC_HOME", Path(__file__).resolve().parents[2]))
DEFAULT_LOCK = HOME / "references" / "references.lock.json"
LDSR_MASK_CANDIDATES = [
    Path(os.environ.get("CONDA_PREFIX", "/opt/envs/worker")) / "lib/R/library/ldsR/extdata/common_snps.parquet",
    Path("/opt/envs/worker/lib/R/library/ldsR/extdata/common_snps.parquet"),
]
PLINK_PREFIX = "plink/1000G.EUR.QC"
COMMON_MAF = 0.05


# --------------------------------------------------------------------------- loading


@dataclass
class Bundle:
    root: Path
    meta: dict

    @property
    def id(self) -> str:
        return self.meta["bundle_id"]

    def plink(self, chrom: int) -> Path:
        return self.root / f"{self.meta['layout']['plink_prefix']}.{chrom}"

    def path(self, key: str) -> Path:
        return self.root / self.meta["layout"][key]

    def reference_snps(self, columns=None) -> pd.DataFrame:
        return pq.read_table(self.path("reference_snps"), columns=columns).to_pandas()

    def regression_snps(self) -> pd.DataFrame:
        return pq.read_table(self.path("regression_snps")).to_pandas()


def load_bundle(root: Path | str, verify: bool = True) -> Bundle:
    root = Path(root).resolve()
    meta_path = root / "BUNDLE.json"
    if not meta_path.exists():
        raise WorkerError("REFERENCE_MISSING", f"no BUNDLE.json in {root}")
    meta = read_json(meta_path)
    if verify:
        verify_bundle_files(root, meta)
    return Bundle(root, meta)


def verify_bundle_files(root: Path, meta: dict) -> None:
    bad = []
    for rel, rec in sorted(meta["files"].items()):
        p = root / rel
        if not p.is_file() or p.stat().st_size != rec["bytes"] or sha256_file(p) != rec["sha256"]:
            bad.append(rel)
    if bad:
        raise WorkerError("REFERENCE_MISMATCH", f"{len(bad)} bundle files fail checksum", bad[:20])


# --------------------------------------------------------------------------- download


def fetch_source(name: str, src: dict, downloads: Path, source_dirs: list[Path],
                 offline: bool, retries: int = 8) -> Path:
    """Return a verified local copy of a locked source, downloading only if needed."""
    downloads.mkdir(parents=True, exist_ok=True)
    basename = src["url"].rsplit("/", 1)[-1]
    final = downloads / f"{src['sha256'][:16]}-{basename}"
    if final.exists() and _source_ok(final, src):
        return final
    for d in source_dirs:
        for cand in (d / basename, d / final.name):
            if cand.is_file() and _source_ok(cand, src):
                shutil.copyfile(cand, final.with_suffix(".part"))
                os.replace(final.with_suffix(".part"), final)
                return final
    if offline:
        raise WorkerError("REFERENCE_MISSING",
                          f"source {name} ({basename}) not available offline", src["url"])
    part = final.with_name(final.name + ".part")
    last = None
    for attempt in range(retries):
        try:
            req = urllib.request.Request(src["url"], headers={"User-Agent": f"agent-ldsc/{__version__}"})
            with urllib.request.urlopen(req, timeout=120) as r, open(part, "wb") as f:
                shutil.copyfileobj(r, f, 1 << 20)
            if _source_ok(part, src):
                os.replace(part, final)
                return final
            last = f"checksum/size mismatch ({part.stat().st_size} bytes)"
        except (urllib.error.URLError, TimeoutError, ConnectionError, OSError) as e:
            last = repr(e)
        time.sleep(min(300, 15 * 2 ** attempt))
    part.unlink(missing_ok=True)
    raise WorkerError("DOWNLOAD_FAILED", f"could not obtain verified {name}", {"url": src["url"], "last_error": last})


def _source_ok(path: Path, src: dict) -> bool:
    if path.stat().st_size != src["bytes"]:
        return False
    if sha256_file(path) != src["sha256"]:
        return False
    return "md5" not in src or md5_file(path) == src["md5"]


def _safe_extract_tar(archive: Path, dest: Path) -> None:
    with tarfile.open(archive) as tf:
        for m in tf.getmembers():
            if m.issym() or m.islnk() or m.name.startswith("/") or ".." in Path(m.name).parts:
                raise WorkerError("REFERENCE_MISMATCH", f"unsafe archive member {m.name}")
        tf.extractall(dest, filter="data")


def _safe_extract_zip(archive: Path, dest: Path) -> None:
    with zipfile.ZipFile(archive) as zf:
        for n in zf.namelist():
            if n.startswith("/") or ".." in Path(n).parts:
                raise WorkerError("REFERENCE_MISMATCH", f"unsafe archive member {n}")
        zf.extractall(dest)


# --------------------------------------------------------------------------- install


def install_bundle(bundle_id: str, cache: Path, lock_path: Path = DEFAULT_LOCK,
                   source_dirs: list[Path] | None = None, offline: bool = False,
                   ldsr_mask: Path | None = None, log=print) -> Path:
    lock = read_json(lock_path)
    if bundle_id not in lock["bundles"]:
        raise WorkerError("REFERENCE_MISSING", f"bundle {bundle_id} not in {lock_path}")
    spec = lock["bundles"][bundle_id]
    refs_dir = cache / "references"
    final = refs_dir / bundle_id
    with dir_lock(refs_dir / f".{bundle_id}.lock", timeout=6 * 3600, stale_after=0):
        if (final / "BUNDLE.json").exists():
            log(f"bundle {bundle_id} already installed; verifying")
            load_bundle(final, verify=True)
            return final
        tmp = refs_dir / f".tmp-{bundle_id}-{os.getpid()}"
        shutil.rmtree(tmp, ignore_errors=True)
        tmp.mkdir(parents=True)
        try:
            paths = {}
            for name in spec["sources"]:
                log(f"fetching {name}")
                paths[name] = fetch_source(name, lock["sources"][name], cache / "downloads",
                                           source_dirs or [], offline)
            meta = _build_eur_phase3(tmp, paths, spec, lock, bundle_id, ldsr_mask, log)
            write_json_atomic(tmp / "BUNDLE.json", meta)
            load_bundle(tmp, verify=True)
            if final.exists():
                shutil.rmtree(final)
            os.replace(tmp, final)
        except BaseException:
            shutil.rmtree(tmp, ignore_errors=True)
            raise
    log(f"installed {final}")
    return final


def _build_eur_phase3(out: Path, paths: dict, spec: dict, lock: dict, bundle_id: str,
                      ldsr_mask: Path | None, log) -> dict:
    exp = spec["expected"]
    stage = out / ".extract"
    stage.mkdir()
    log("extracting sldsc_ref")
    _safe_extract_tar(paths["sldsc_ref"], stage)
    src = stage / "sldsc_ref"

    (out / "plink").mkdir()
    (out / "frq").mkdir()
    (out / "weights").mkdir()
    for c in AUTOSOMES:
        for ext in ("bed", "bim", "fam"):
            os.replace(src / f"1000G_EUR_Phase3_plink/1000G.EUR.QC.{c}.{ext}", out / f"{PLINK_PREFIX}.{c}.{ext}")
        os.replace(src / f"1000G_Phase3_frq/1000G.EUR.QC.{c}.frq", out / f"frq/1000G.EUR.QC.{c}.frq")
        os.replace(src / f"1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC.{c}.l2.ldscore.gz",
                   out / f"weights/weights.hm3_noMHC.{c}.l2.ldscore.gz")
    (out / "snps").mkdir()
    os.replace(src / "hm_snp.txt", out / "snps/print_snps.txt")

    log("deriving SNP indexes")
    bim = pd.concat([read_bim(out / f"{PLINK_PREFIX}.{c}.bim") for c in AUTOSOMES], ignore_index=True)
    frq = pd.concat([pd.read_csv(out / f"frq/1000G.EUR.QC.{c}.frq", sep=r"\s+") for c in AUTOSOMES],
                    ignore_index=True)
    checks = {}
    _require(len(bim) == exp["reference_snps"], "reference SNP count", len(bim))
    _require(not bim.SNP.duplicated().any(), "duplicate reference SNP ids")
    _require(bool((bim.SNP.values == frq.SNP.values).all()), "frq/bim SNP order")
    _require(bool(((bim.A1.values == frq.A1.values) & (bim.A2.values == frq.A2.values)).all()), "frq/bim alleles")
    _require(all(bool((g.BP.diff().dropna() > 0).all()) for _, g in bim.groupby("CHR")),
             "BIM positions strictly increasing within chromosome")
    _require(int((frq.MAF.values == 0).sum()) == 0, "monomorphic SNPs present")
    fam_n = sum(1 for _ in open(out / f"{PLINK_PREFIX}.1.fam"))
    _require(fam_n == exp["individuals"], "individual count", fam_n)
    for c in AUTOSOMES[1:]:
        _require(sum(1 for _ in open(out / f"{PLINK_PREFIX}.{c}.fam")) == fam_n, f"fam chr{c} count")
    common = frq.MAF.values > COMMON_MAF
    _require(int(common.sum()) == exp["common_snps_maf_gt_0.05"], "common SNP count", int(common.sum()))
    checks["frq_matches_bim"] = True

    ref = pd.DataFrame({"CHR": bim.CHR.astype("int8"), "SNP": bim.SNP, "BP": bim.BP.astype("int64"),
                        "CM": bim.CM.astype("float64"), "A1": bim.A1, "A2": bim.A2,
                        "MAF": frq.MAF.astype("float64"), "common": common})
    pq.write_table(pa.Table.from_pandas(ref, preserve_index=False), out / "snps/reference_snps.parquet")

    printed = set(pd.read_csv(out / "snps/print_snps.txt", header=None, dtype=str)[0])
    reg = ref.loc[ref.SNP.isin(printed), ["CHR", "SNP", "BP"]].reset_index(drop=True)
    _require(len(reg) == exp["regression_snps"], "regression SNP count", len(reg))
    pq.write_table(pa.Table.from_pandas(reg, preserve_index=False), out / "snps/regression_snps.parquet")

    log("checking ldsR compatibility targets")
    mask_path = ldsr_mask or next((p for p in LDSR_MASK_CANDIDATES if p.exists()), None)
    _require(mask_path is not None, "ldsR common_snps.parquet not found (pass --ldsr-mask)")
    mask = pq.read_table(mask_path).column("common").to_numpy()
    _require(len(mask) == len(ref) and bool((mask == common).all()), "ldsR common mask positional identity")
    checks["ldsr_common_mask_identical"] = {"path": str(mask_path), "sha256": sha256_file(mask_path)}

    _safe_extract_zip(paths["ldsr_baseline_v1.2"], stage)
    base_dir = out / "ldsr_baseline_v1.2"
    os.replace(stage / "baseline_model-v1.2", base_dir)
    base_ld = pq.read_table(base_dir / "ld.parquet", columns=["SNP"]).column(0).to_pylist()
    base_ref = pq.read_table(base_dir / "annot_ref.parquet", columns=["SNP"]).column(0).to_pylist()
    _require(base_ld == reg.SNP.tolist(), "ldsR baseline ld.parquet SNP order == regression index")
    _require(base_ref == ref.SNP.tolist(), "ldsR baseline annot_ref.parquet SNP order == reference index")
    checks["ldsr_baseline_order_identical"] = True

    log("deriving gene tables")
    (out / "genes").mkdir()
    genes = parse_gencode_genes(paths["gencode_v50lift37_basic"])
    hgnc = parse_hgnc(paths["hgnc_2026_10_06"])
    ens2hgnc = (hgnc.loc[hgnc.ensembl_gene_id != "", ["ensembl_gene_id", "hgnc_id"]]
                .drop_duplicates("ensembl_gene_id"))
    genes = genes.merge(ens2hgnc, left_on="gene_id", right_on="ensembl_gene_id", how="left") \
                 .drop(columns="ensembl_gene_id").fillna({"hgnc_id": ""})
    genes.to_csv(out / "genes/genes.tsv.gz", sep="\t", index=False, compression={"method": "gzip", "mtime": 0})
    hgnc.to_csv(out / "genes/hgnc.tsv.gz", sep="\t", index=False, compression={"method": "gzip", "mtime": 0})

    (out / "chains").mkdir()
    shutil.copyfile(paths["ucsc_hg38ToHg19"], out / "chains/hg38ToHg19.over.chain.gz")
    shutil.copyfile(paths["ucsc_hg19ToHg38"], out / "chains/hg19ToHg38.over.chain.gz")
    sizes = pd.read_csv(paths["ucsc_hg19_chrom_sizes"], sep="\t", header=None, names=["chrom", "length"])
    sizes = sizes[sizes.chrom.isin([f"chr{c}" for c in AUTOSOMES])]
    sizes["chrom_num"] = sizes.chrom.str[3:].astype(int)
    sizes = sizes.sort_values("chrom_num")[["chrom", "length"]]
    _require(len(sizes) == 22, "hg19 autosome sizes")
    (out / "chrom_sizes").mkdir()
    sizes.to_csv(out / "chrom_sizes/GRCh37.autosomes.tsv", sep="\t", index=False)
    shutil.rmtree(stage)

    files = {str(p.relative_to(out)): {"bytes": p.stat().st_size, "sha256": sha256_file(p)}
             for p in sorted(out.rglob("*")) if p.is_file()}
    return {
        "schema_version": 1,
        "bundle_id": bundle_id,
        "builder": {"agent_ldsc_worker": __version__, "source_sha": os.environ.get("AGENT_LDSC_SOURCE_SHA", "unknown")},
        "created": now(),
        "genome_build": spec["genome_build"],
        "ancestry": spec["ancestry"],
        "description": spec["description"],
        "sources": {n: lock["sources"][n] for n in spec["sources"]},
        "layout": {
            "plink_prefix": PLINK_PREFIX,
            "reference_snps": "snps/reference_snps.parquet",
            "regression_snps": "snps/regression_snps.parquet",
            "print_snps": "snps/print_snps.txt",
            "genes": "genes/genes.tsv.gz",
            "hgnc": "genes/hgnc.tsv.gz",
            "chrom_sizes": "chrom_sizes/GRCh37.autosomes.tsv",
            "chain_hg38_to_ref": "chains/hg38ToHg19.over.chain.gz",
            "ldsr_baseline": "ldsr_baseline_v1.2",
            "weights_prefix": "weights/weights.hm3_noMHC",
        },
        "snp_universe": {
            "reference_snps": len(ref),
            "regression_snps": len(reg),
            "common_snps": int(common.sum()),
            "common_definition": "MAF > 0.05 from 1000G_Phase3_frq (identical to ldsR common_snps mask)",
            "individuals": fam_n,
            "reference_order_sha256": sha256_strings(ref.SNP),
            "regression_order_sha256": sha256_strings(reg.SNP),
            "common_mask_sha256": sha256_strings(common.astype(int)),
            "per_chromosome": {str(c): {"reference": int((ref.CHR == c).sum()), "regression": int((reg.CHR == c).sum())}
                               for c in AUTOSOMES},
        },
        "genes": {"n_genes": len(genes), "n_autosomal": int(genes.chrom_num.between(1, 22).sum()),
                  "coordinate_source": "gencode_v50lift37_basic", "identifier_source": "hgnc_2026_10_06"},
        "checks": checks,
        "files": files,
    }


def _require(ok: bool, what: str, observed=None) -> None:
    if not ok:
        raise WorkerError("REFERENCE_MISMATCH", f"reference check failed: {what}", observed)


def read_bim(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep=r"\s+", header=None, names=["CHR", "SNP", "CM", "BP", "A1", "A2"],
                       dtype={"CHR": int, "SNP": str, "CM": float, "BP": int, "A1": str, "A2": str})


_ATTR = re.compile(r'(\w+) "?([^";]*)"?;')


def parse_gencode_genes(gtf: Path) -> pd.DataFrame:
    rows = []
    with gzip.open(gtf, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            p = line.rstrip("\n").split("\t")
            if p[2] != "gene":
                continue
            a = dict(_ATTR.findall(p[8]))
            gid_v = a["gene_id"]
            rows.append({
                "gene_id": gid_v.split(".")[0] + ("_PAR_Y" if gid_v.endswith("_PAR_Y") else ""),
                "gene_id_version": gid_v,
                "gene_name": a.get("gene_name", ""),
                "gene_type": a.get("gene_type", ""),
                "chrom": p[0],
                "start": int(p[3]),
                "end": int(p[4]),
                "strand": p[6],
                "level": a.get("level", ""),
                "remap_status": a.get("remap_status", ""),
                "remap_num_mappings": a.get("remap_num_mappings", ""),
            })
    df = pd.DataFrame(rows)
    df["chrom_num"] = [c if c is not None else 0 for c in map(_chrom_num, df.chrom)]
    return df


def _chrom_num(chrom: str):
    v = chrom[3:] if chrom.startswith("chr") else chrom
    return int(v) if v.isdigit() else None


def parse_hgnc(path: Path) -> pd.DataFrame:
    cols = ["hgnc_id", "symbol", "status", "locus_type", "alias_symbol", "prev_symbol", "ensembl_gene_id"]
    df = pd.read_csv(path, sep="\t", usecols=cols, dtype=str, keep_default_na=False)
    return df[cols]
