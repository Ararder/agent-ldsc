# Provenance of the recovered method

The LD-score method reproduces the lab's historical Dardel workflow (2023–2026), recovered
from untracked scripts in `generate-ldscores` and from `stratified-ldscore-regression/scripts/compute_ldscores_container.sh`.

| Historical step | Historical implementation | Here |
|---|---|---|
| hg38 → hg19 | `liftOver` (image `vanallenlab/liftover:latest`, default `-minMatch=0.95`), UCSC `hg38ToHg19.over.chain.gz` | bioconda `ucsc-liftover` 482, same parameters. The chain is byte-identical (sha256 `14a712e8…b1f505`) |
| BED → annotation | `make_annot.py --bed-file … --bimfile …` (merge, `[BP-1, BP)` intersect) | `intervals.snp_membership`, same semantics, tested at the boundaries |
| LD scores | `ldsc.py --bfile 1000G.EUR.QC.<chr> --ld-wind-cm 1 --annot … --thin-annot --print-snps hm_snp.txt` in image `arvhar/ldsc:latest` (ldsc `aa33296`, Python 2.7.13, numpy 1.16.0, scipy 0.18.0, pandas 0.20.0, bitarray 0.8.3) | same command and commit. pandas is 0.20.3 (0.20.0 has no Linux wheel) |
| Reference | Zenodo 8367200 `sldsc_ref.tar.gz` (MD5 `442741559abab2680ad431021244e1f7`) | same archive, now also pinned by sha256 and with its contents audited (`references/references.lock.json`) |
| ldsR conversion | `ldsc_to_parquet` + `reorder_data` to `ld_order`/`annot_ref_order` | strict exporter. Orders are proven identical to the ldsR baseline-v1.2 deposit (Zenodo 15194124) |

Archive audit (2026-10-08): 9,997,231 autosomal SNPs, 489 EUR individuals, no duplicate rsIDs,
no monomorphic SNPs, and BIM positions strictly increasing. frq rows equal BIM rows (SNP and
alleles). 5,961,159 SNPs have MAF > 0.05, which matches ldsR's `common_snps` mask position by
position. 1,190,321 of the 1,217,311 HM3 print-list SNPs are in the panel. Their BIM order
equals the baseline `.l2.ldscore` and ldsR `ld.parquet` order.

Comparison target for real data: `Data/ldsR_ldscores/ss3x_fixed500_19` on Dardel, built from
`generate-ldscores/workflow/ss3x_fixed500_all/<celltype>.bed` (hg38). Spot check, `inh_SNCG_2` chr22:
the canonical intervals equal the historical lifted BED, LD scores are identical (max |Δ| = 0 over
17,489 SNPs), and M equals the historical `annot_ref` count.

The historical `topx` array workflow had chromosome 1 produce the lifted BED while other array
tasks polled for it. Here, normalization is a completed stage before any LD task starts.
