"""Generate a tiny, transparent 22-chromosome reference bundle for tests.

Design (so an independent numpy oracle can reproduce ldsc exactly, see tests/oracle.py):
- Every chromosome has 2-3 LD clusters of exactly 50 SNPs (ldsc's chunk size), aligned to
  chunk boundaries; each cluster spans 0.441 cM and clusters are 2 cM apart. ldsc's chunked
  window then equals "all SNP pairs within a cluster".
- 61 individuals (122 alleles, so no SNP has MAF exactly 0.05); no missing genotypes; no
  monomorphic SNPs; a few rare SNPs (MAF < 0.05) per chromosome to exercise M_5_50.
- Gene/HGNC tables and a small hg38->hg19 chain with a deleted segment for mapping tests.

Usage: python make_fixture_bundle.py OUT_DIR
"""

from __future__ import annotations

import gzip
import hashlib
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq

N_IND = 61
CLUSTER = 50
CHROM_LEN = 3_000_000
SEED = 20261008
PREFIX = "plink/fixture"


def sha(values) -> str:
    h = hashlib.sha256()
    for v in values:
        h.update(str(v).encode())
        h.update(b"\n")
    return h.hexdigest()


def simulate_cluster(rng: np.random.Generator, n_snps: int) -> np.ndarray:
    """Allele-1 counts (0/1/2) for N_IND individuals from 4 noisy founder haplotypes."""
    while True:
        founders = rng.random((4, n_snps)) < rng.uniform(0.15, 0.6, n_snps)
        hap = founders[rng.integers(0, 4, 2 * N_IND)]
        hap ^= rng.random(hap.shape) < 0.04
        rare = rng.choice(n_snps, 3, replace=False)
        hap[:, rare] = False
        for j in rare:
            hap[rng.choice(2 * N_IND, rng.integers(1, 5), replace=False), j] = True
        g = hap[0::2].astype(int) + hap[1::2].astype(int)
        freq = g.sum(axis=0) / (2 * N_IND)
        if ((freq > 0) & (freq < 1)).all():
            return g


def write_bed(path: Path, geno: np.ndarray) -> None:
    """PLINK 1 SNP-major .bed; geno = counts of allele A1 (n_ind x n_snp)."""
    code = {2: 0b00, 1: 0b10, 0: 0b11}  # 00 hom A1, 10 het, 11 hom A2 (01 = missing)
    out = bytearray(b"\x6c\x1b\x01")
    for j in range(geno.shape[1]):
        col = geno[:, j]
        for i in range(0, N_IND, 4):
            byte = 0
            for k, gi in enumerate(col[i:i + 4]):
                byte |= code[int(gi)] << (2 * k)
            out.append(byte)
    path.write_bytes(bytes(out))


def main(out: Path) -> None:
    rng = np.random.default_rng(SEED)
    (out / "plink").mkdir(parents=True)
    (out / "snps").mkdir()
    ref_rows, genos = [], {}
    for c in range(1, 23):
        n_clusters = 3 if c in (1, 2, 6) else 2
        g = np.hstack([simulate_cluster(rng, CLUSTER) for _ in range(n_clusters)])
        genos[c] = g
        for k in range(n_clusters):
            for j in range(CLUSTER):
                bp = 100_000 + k * 800_000 + j * 2_000 + 1
                ref_rows.append({"CHR": c, "SNP": f"rs{c}_{k}_{j}", "CM": round(2.0 * k + 0.009 * j, 6),
                                 "BP": bp, "A1": "A", "A2": "G"})
        bim = pd.DataFrame([r for r in ref_rows if r["CHR"] == c])
        bim[["CHR", "SNP", "CM", "BP", "A1", "A2"]].to_csv(out / f"{PREFIX}.{c}.bim", sep="\t", header=False, index=False)
        with open(out / f"{PREFIX}.{c}.fam", "w") as f:
            for i in range(N_IND):
                f.write(f"ID{i} ID{i} 0 0 0 -9\n")
        write_bed(out / f"{PREFIX}.{c}.bed", g)
    ref = pd.DataFrame(ref_rows)
    p1 = np.concatenate([genos[c].sum(axis=0) / (2 * N_IND) for c in range(1, 23)])
    ref["MAF"] = np.minimum(p1, 1 - p1)
    ref["common"] = ref.MAF > 0.05
    ref = ref[["CHR", "SNP", "BP", "CM", "A1", "A2", "MAF", "common"]].astype({"CHR": "int8"})
    pq.write_table(pa.Table.from_pandas(ref, preserve_index=False), out / "snps/reference_snps.parquet")

    # print list: ~70% of reference SNPs plus ids absent from the panel
    keep = rng.random(len(ref)) < 0.7
    printed = list(ref.SNP[keep]) + ["rs_not_in_panel_1", "rs_not_in_panel_2"]
    (out / "snps/print_snps.txt").write_text("\n".join(printed) + "\n")
    reg = ref.loc[keep, ["CHR", "SNP", "BP"]].reset_index(drop=True)
    pq.write_table(pa.Table.from_pandas(reg, preserve_index=False), out / "snps/regression_snps.parquet")

    (out / "chrom_sizes").mkdir()
    pd.DataFrame({"chrom": [f"chr{c}" for c in range(1, 23)], "length": CHROM_LEN}) \
        .to_csv(out / "chrom_sizes/GRCh37.autosomes.tsv", sep="\t", index=False)

    # genes (1-based inclusive spans) and HGNC symbols
    (out / "genes").mkdir()
    genes = pd.DataFrame([
        ("ENSG00000000001", "ENSG00000000001.5_1", "GENEA", "protein_coding", "chr1", 100_001, 140_000, "+", "2", "full_contig", "1", "HGNC:1"),
        ("ENSG00000000002", "ENSG00000000002.3_1", "GENEB", "protein_coding", "chr2", 900_001, 920_000, "-", "2", "full_contig", "1", "HGNC:2"),
        ("ENSG00000000003", "ENSG00000000003.1_1", "GENEC", "lncRNA", "chr3", 100_001, 100_500, "+", "2", "partial", "1", "HGNC:3"),
        ("ENSG00000000004", "ENSG00000000004.1_1", "GENED", "protein_coding", "chrX", 100_001, 200_000, "+", "2", "full_contig", "1", "HGNC:4"),
        ("ENSG00000000005", "ENSG00000000005.2_1", "GENEE", "protein_coding", "chr22", 2_850_001, 2_990_000, "+", "2", "full_contig", "1", "HGNC:5"),
    ], columns=["gene_id", "gene_id_version", "gene_name", "gene_type", "chrom", "start", "end", "strand",
                "level", "remap_status", "remap_num_mappings", "hgnc_id"])
    genes["chrom_num"] = [int(c[3:]) if c[3:].isdigit() else 0 for c in genes.chrom]
    genes.to_csv(out / "genes/genes.tsv.gz", sep="\t", index=False, compression={"method": "gzip", "mtime": 0})
    hgnc = pd.DataFrame([
        ("HGNC:1", "GENEA", "Approved", "gene with protein product", "AMBIG|OLDA", "", "ENSG00000000001"),
        ("HGNC:2", "GENEB", "Approved", "gene with protein product", "AMBIG", "OLDB", "ENSG00000000002"),
        ("HGNC:3", "GENEC", "Approved", "RNA, long non-coding", "", "", "ENSG00000000003"),
        ("HGNC:4", "GENED", "Approved", "gene with protein product", "", "", "ENSG00000000004"),
        ("HGNC:5", "GENEE", "Approved", "gene with protein product", "", "", "ENSG00000000005"),
        ("HGNC:6", "GENEF", "Approved", "gene with protein product", "", "", ""),
    ], columns=["hgnc_id", "symbol", "status", "locus_type", "alias_symbol", "prev_symbol", "ensembl_gene_id"])
    hgnc.to_csv(out / "genes/hgnc.tsv.gz", sep="\t", index=False, compression={"method": "gzip", "mtime": 0})

    # hg38 -> hg19 chain: chr1 [0,500k) -> +1000; [500k,510k) deleted; [510k,1.01M) -> [501k,1.001M);
    # chr2..chr22 identity over the whole fixture length.
    (out / "chains").mkdir()
    lines = [f"chain 1000 chr1 {CHROM_LEN} + 0 1010000 chr1 {CHROM_LEN} + 1000 1001000 1",
             "500000\t10000\t0", "500000", ""]
    for c in range(2, 23):
        lines += [f"chain 1000 chr{c} {CHROM_LEN} + 0 {CHROM_LEN} chr{c} {CHROM_LEN} + 0 {CHROM_LEN} {c}",
                  str(CHROM_LEN), ""]
    with gzip.GzipFile(out / "chains/hg38ToHg19.over.chain.gz", "wb", mtime=0) as f:
        f.write("\n".join(lines).encode())

    files = {str(p.relative_to(out)): {"bytes": p.stat().st_size, "sha256": hashlib.sha256(p.read_bytes()).hexdigest()}
             for p in sorted(out.rglob("*")) if p.is_file()}
    meta = {
        "schema_version": 1, "bundle_id": "fixture-grch37-v1", "builder": {"agent_ldsc_worker": "fixture"},
        "created": "2026-10-08T00:00:00Z", "genome_build": "GRCh37", "ancestry": "synthetic",
        "description": "Synthetic test fixture; not a scientific reference.", "sources": {},
        "layout": {"plink_prefix": PREFIX, "reference_snps": "snps/reference_snps.parquet",
                   "regression_snps": "snps/regression_snps.parquet", "print_snps": "snps/print_snps.txt",
                   "genes": "genes/genes.tsv.gz", "hgnc": "genes/hgnc.tsv.gz",
                   "chrom_sizes": "chrom_sizes/GRCh37.autosomes.tsv",
                   "chain_hg38_to_ref": "chains/hg38ToHg19.over.chain.gz"},
        "snp_universe": {"reference_snps": len(ref), "regression_snps": len(reg), "common_snps": int(ref.common.sum()),
                         "common_definition": "MAF > 0.05", "individuals": N_IND,
                         "reference_order_sha256": sha(ref.SNP), "regression_order_sha256": sha(reg.SNP),
                         "common_mask_sha256": sha(ref.common.astype(int))},
        "genes": {"n_genes": len(genes), "coordinate_source": "fixture", "identifier_source": "fixture"},
        "checks": {}, "files": files,
    }
    (out / "BUNDLE.json").write_text(json.dumps(meta, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main(Path(sys.argv[1]))
