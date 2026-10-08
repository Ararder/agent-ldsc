# Input/output contract

## Request (`schemas/request.schema.json`)

```json
{
  "schema_version": 1,
  "reference_id": "eur-phase3-grch37-v1",
  "scientific_profile": "strict-v1",
  "annotations": [
    {"id": "atac_neurons", "type": "bed", "file": "peaks.bed", "species": "human", "build": "GRCh38",
     "allow_drop": ["liftover_unmapped", "non_autosomal"]},
    {"id": "marker_genes", "type": "genes", "file": "genes.txt", "species": "human",
     "id_type": "symbol", "gene_model": "gene_span", "flank_bp": 100000}
  ]
}
```

- `id`: `^[A-Za-z][A-Za-z0-9_.-]{0,63}$`, unique (also case-insensitively), not a reserved name
  (`SNP`, `CHR`, `BP`, `CM`, `annot`, `m`, `m50`, `MAF`, `common`, `base`, `L2`) and not ending in `L2`.
  Column order in every output file follows the request order.
- `file`: relative to the request file. The launcher copies it into the run's `input/files/`.
- BED: tab or whitespace separated, at least chrom/start/end, 0-based half-open. Lines starting
  with `#`, `track` or `browser` are skipped; extra columns are ignored, except that column 4
  is kept in the audit table as the record name. `build` is required (`GRCh37`/`hg19` or
  `GRCh38`/`hg38`). It is never guessed.
- Genes: one identifier per line. Blank lines and `#` lines are skipped. `header: true` skips the
  first remaining line. `id_type` is required (`symbol` or `ensembl_gene`). `flank_bp` is
  required: there is no default.
- `allow_drop` lists the loss categories you accept. Any other loss fails the run with `MAPPING_LOSS`:
  - `non_autosomal`: input on chrX/Y/M/alt contigs, or genes located there.
  - `liftover_unmapped`: liftOver deleted, split, partially deleted, or mapped to a non-autosomal contig.
  - `gene_unmapped`: symbol not found, no Ensembl ID, not in the gene model, invalid ID, or a
    GENCODE remap status excluded by the profile.
  - Ambiguous gene symbols always fail with `GENE_AMBIGUOUS`, and the error lists the candidates.

## Scientific profile `strict-v1` (`profiles/scientific/strict-v1.json`)

1 cM LD window; thin binary annotations; HM3 print list; liftOver `-minMatch=0.95`, no
multiple mappings. Genes use the GENCODE v50lift37 `gene` record span, symmetrically flanked
and clipped to chromosome ends. Symbols match exact HGNC approved symbols, or a unique
alias/previous symbol. GENCODE `remap_status` partial, or multiple remaps, counts as unmapped.
Any change requires a new profile id.

SNP membership: a SNP at 1-based BIM position `p` is in BED `[s, e)` iff `s <= p-1 < e`
(identical to `make_annot.py`). Intervals within an annotation are unioned. Overlapping peaks
therefore never count twice.

## Run directory

```text
<work_root>/<run-id>/
  input/request.json, input/files/*     staged by the launcher
  resolved-request.json                 request + profile + input/canonical hashes + scientific key
  state.json                            status, stage, task counts, attempts (host, Slurm job id)
  normalized/<id>.bed                   canonical merged intervals (reference build, chrN, sorted 1..22)
  normalized/<id>.audit.tsv.gz          every input record: status, reason, mapped coordinates
  normalized/<id>.summary.json          counts and base pairs in/kept/dropped/merged
  tasks/<id>/chr<N>/                    thin annot, ldsc outputs + log, COMPLETE.json (input key + output hashes)
  logs/worker.log, logs/error.json, logs/slurm-<job>.out
  output/ldsR/                          final, published atomically
  COMPLETE.json                         written last: scientific key + manifest sha256
```

Status values: `initializing`, `running`, `failed`, `cancelled`, `completed`.
Exit codes: 0 complete, 2 failed (`logs/error.json`), 3 run locked, 4 interrupted (resumable),
5 unexpected crash.

Error codes: `INPUT_INVALID`, `BUILD_UNSUPPORTED`, `MAPPING_LOSS`, `GENE_AMBIGUOUS`,
`REFERENCE_MISSING`, `REFERENCE_MISMATCH`, `DOWNLOAD_FAILED`, `LDSC_TASK_FAILED`,
`RESOURCE_LIMIT`, `EXPORT_INVALID`, `RUN_LOCKED`, `STATE_INVALID`.

## Output `output/ldsR/`

| File | Contract |
|---|---|
| `ld.parquet` | `SNP` (string) plus one float64 column per annotation; 1,190,321 regression SNPs in BIM order chr1→22 (equal to the ldsR baseline-v1.2 order) |
| `annot.parquet` | `annot` (string), `m`, `m50` (float64); one row per annotation in column order |
| `annot_ref.parquet` | `SNP` plus one int32 0/1 column per annotation; all 9,997,231 reference SNPs in BIM order (positionally aligned with ldsR's `common_snps` mask) |
| `<id>.canonical.bed` | exact intervals that defined each annotation |
| `manifest.json` | request, profile, reference bundle identity and SNP-order hashes, software (image digest, ldsc/ldsR commits), per-annotation normalization summary and counts, validation results, `ldsr_compatible`, output sha256 |
| `qc.json`, `qc.tsv` | records in/dropped by reason, merged intervals/bp, m, m50, annotated regression SNPs, zero-coverage chromosomes, per-task M/M_5_50 |

Checks before publication: the SNP order of every task equals the bundle index. LD scores
are finite. `.l2.M` equals the annotated SNP count and `.l2.M_5_50` equals the frq-based common
count (so genotype MAF and the ldsR mask agree). `m50 <= m`. 22 chromosomes are present. The
ldsR gate (`R/validate_ldsR.R`) reads the directory with `parse_parquet_dir(read_ref = TRUE)` and
recomputes `m`/`m50` from `annot_ref` and ldsR's bundled mask. `ldsr_compatible` is true only
when the bundle passed the mask/baseline identity checks and this gate passed. Zero coverage on
a single chromosome is valid (and reported). Zero coverage genome-wide fails.

## Numerical tolerances

- Mapping, membership and counts: exact.
- Docker and Apptainer: exact semantic equality of the parquet tables (`container/compare_ldsr.py`).
- Comparison with an independent numpy oracle (`tests/oracle.py`): |Δ| ≤ 5e-4, because ldsc writes
  LD scores with `%.3f`. The plan's provisional 1e-8 + 1e-6·|b| applies only between runs of
  the same tool, where values are in fact identical.

## LDSC window note

ldsc rounds each SNP's window up to its 50-SNP chunk (`ceil(size / c) * c`), and it computes
whole-chromosome LD when the first out-of-window SNP falls past the chromosome end. This is
the recovered historical method and is reproduced, not corrected. The test fixture uses
chunk-aligned LD clusters so that the oracle and ldsc define the same pair set.
