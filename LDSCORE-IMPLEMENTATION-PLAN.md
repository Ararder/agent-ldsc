# LD-score calculation: implementation handoff draft

2026-10-08. Scope: implement annotation-to-LD-score generation and validated ldsR export. This document takes precedence over the broader DESIGN-DRAFT.md for this component. Proposed commands and schemas below do not yet exist. No container was built or job submitted during planning.

## 1. Outcome and scope

A user submits one request containing human genes or human genomic intervals. The selected local or remote execution site stages the input and references, creates the run structure, maps annotations into the reference build, calculates LD scores for chromosomes 1–22, and writes an ldsR-compatible artifact with provenance and QC. The same scientific worker runs locally and on Slurm.

In scope: human gene symbols/Ensembl gene IDs, BED intervals (including ATAC peaks), GRCh37 and GRCh38 input intervals, binary annotations, multiple named annotation sets, container build/release, pinned reference cache, local execution, SSH transport, one Slurm job per request, resumability, ldsR conversion.

Out of scope: GWAS acquisition/preparation, regressions or enrichment reporting, automated scientific interpretation, quantitative/continuous annotations, sex chromosomes, other ancestries/reference panels, mouse input. Mouse orthology and cross-species interval mapping are a separate second phase. SCZ/BIP/height files already retrieved are not dependencies of this component or its CI.

## 2. Architecture decision

Use one OCI image containing all scientific and workflow-worker software. A small host launcher handles container startup, SSH/rsync, and scheduler submission only. It must not perform biological transformations.

Separate three concepts:

- Transport: current machine or an explicitly configured SSH alias.
- Executor: direct process or Slurm allocation.
- Runtime: Docker or Apptainer.

A remote workstation can use direct execution; Slurm is not synonymous with remote. The host requires the runtime and, as appropriate, shell, SSH/rsync and Slurm client. The worker requires no host R, Python, conda, bedtools, or liftover installation. No container-in-container execution.

For v1, a Slurm request receives one allocation. Inside it, the worker schedules independent chromosome/annotation tasks with bounded concurrency. Preparation finishes before LD tasks start; conversion waits for all tasks to validate. This is the same process graph as direct execution. Avoid nested sbatch calls and separate job arrays initially. Persist per-task completion so a wall-time interruption can resume in a new allocation. Arrays can be added later without changing the scientific task interface.

## 3. Software image and GitHub Actions

Yes: build and publish through GitHub Actions to `ghcr.io/ararder/agent-ldsc` (proposed package name; repository ownership must be established during implementation).

Image contents:

- Pinned Linux base by digest.
- Modern Python for worker orchestration, request validation, cache management and structured status.
- A separately isolated, locked legacy Python environment for a pinned original bulik/ldsc commit, initially preserving the recovered calculation method. Lock exact dependency artifacts/checksums; never solve a floating conda environment on a user's machine. A maintained Python 3 replacement is allowed only after numerical-equivalence validation, not as an untested substitution.
- Pinned UCSC liftOver binary or reproducible source build, plus bedtools if used for interval normalization.
- Pinned R, Arrow and the packages required for export; pinned ldsR commit for compatibility tests. Use the existing inspected commit as the initial candidate: e880755004946770b0d8a048f970d513fafc16c8.
- This project's CLI, scripts, schemas and tool-version inventory.

All software installations happen during image build. Runtime performs no pip/conda/R package installation. Reference panels, gene tables and chain files are versioned data outside the image. Do not embed lab GWAS files, credentials, local datasets or source-snapshots wholesale; use an explicit .dockerignore and selected build context. Verify redistribution terms for included binaries.

Initial supported image platform: linux/amd64. Docker on Intel Linux and Apptainer on Dardel are initial execution targets. Apple Silicon may run amd64 emulation, but do not claim native ARM support or comparable performance until separately tested.

CI workflow:

1. Pull requests: build image, run schema/unit/fixture tests and a small container end-to-end test; do not publish release images or expose publishing credentials to fork code.
2. Main branch: optionally publish development images tagged by source SHA.
3. Version release: build once, test that image, publish version and source-SHA tags, and record the immutable OCI digest. Downstream execution uses the digest, never `latest`.
4. Generate SBOM/build provenance where supported. Pin third-party Actions by commit SHA; use minimal permissions (`contents: read`, `packages: write` for publishing, appropriate attestation permissions only when used).
5. Run Docker/Apptainer fixture parity in a suitable Linux CI job. Apptainer converts the exact OCI digest to a cached SIF; record both source OCI digest and resulting SIF SHA256. A prebuilt SIF is optional, not a second scientific build recipe.
6. Full reference/22-chromosome acceptance is a separately triggered release validation; tiny CI fixtures must not be presented as full-panel validation. No automatic submission to Dardel from GitHub Actions.

A full toolchain may be too large for a naive legacy build on a hosted runner. The first implementation spike must prove a clean build, measure disk/time/image size and resolve unavailable legacy dependencies before coding the rest of the workflow around that recipe.

Sources: [GitHub image publishing](https://docs.github.com/en/actions/tutorials/publish-packages/publish-docker-images), [Docker build attestations](https://docs.docker.com/build/ci/github-actions/attestations/), [Apptainer OCI registries](https://apptainer.org/docs/user/main/registry.html).

## 4. Reference lock and staging

Preserve the recovered source pin:

- https://zenodo.org/records/8367200/files/sldsc_ref.tar.gz
- 712902644 bytes; MD5 `442741559abab2680ad431021244e1f7`.
- Metadata agrees with the historical Dardel script. Archive contents have not yet been audited in this session.

First implementation task: download and verify the archive; inventory actual contents, genome build, panel individuals/SNP filters, frequencies, genetic positions, HM3 output list and chromosome coverage. Establish the exact ordered reference SNP identities. Add SHA256 for the archive and extracted assets. Do not infer these from directory names or assume the archive contains every required asset.

A committed `references.lock.json` must describe:

- Exact source URL/version, byte size, MD5 where supplied, SHA256 and extraction paths.
- PLINK BED/BIM/FAM files for all 22 chromosomes; panel ancestry and genome build.
- Explicit ordered full-reference SNP index and ordered regression-output SNP index, with hashes.
- Allele-frequency source and the precise common-SNP definition used for M_5_50.
- Human gene coordinates, aliases and chromosome lengths from a pinned release. For v1 resolve genes directly in the reference build rather than lifting gene spans. Identifier mapping release and coordinate assembly are separate metadata fields.
- Supported chain files with source and checksums; copied hg38ToHg19 and hg19ToHg38 chains are candidates, not evidence of a pinned tool version.
- Compatibility target: pinned ldsR and the baseline-1.2 reference deposit at https://zenodo.org/records/15194124. Verify full annotation SNP order and common-SNP mask against this baseline before declaring overlap compatibility. Baseline LD scores are not needed to calculate a custom annotation, but this compatibility gate is needed for downstream use.

Data staging is automatic on a configured execution site before job submission when missing. Use one site cache, file locks, temporary download/extraction directories, checksum enforcement and atomic publication of completed bundles. Validate archive paths before extraction. A cached bundle has an immutable ID and validation manifest; existence alone is insufficient. Fail on corruption rather than silently mixing assets or switching references.

Compute nodes must be able to run offline: stage the image and all data on a network-enabled host first. If that is unavailable, support explicit offline cache installation. No task should download its own copy. On Dardel, caches and work directories belong on project storage, not home. Existing reference directories can be imported after validation rather than downloaded again.

## 5. Request and CLI contract

Proposed user interface:

```sh
agent-ldsc submit request.json --profile local
agent-ldsc submit request.json --profile dardel
agent-ldsc status RUN_ID --json
agent-ldsc resume RUN_ID --profile dardel
agent-ldsc fetch RUN_ID --destination results/
```

Worker interface inside the image:

```sh
agent-ldsc-worker run --request /input/request.json --refs /refs --work /work
agent-ldsc-worker validate-output /work/output/ldsR
```

Request example (profile and bundle IDs are proposed names):

```json
{
  "schema_version": 1,
  "reference_id": "eur-phase3-grch37-v1",
  "annotations": [
    {"id": "atac_neurons", "type": "bed", "file": "inputs/peaks.bed", "species": "human", "build": "GRCh38"},
    {"id": "marker_genes", "type": "genes", "file": "inputs/genes.txt", "species": "human", "id_type": "ensembl_gene", "gene_model": "gene_span", "flank_bp": 100000}
  ]
}
```

The example's 100 kb flank is explicit illustrative configuration, not an agreed scientific default. Require it or a named versioned annotation profile; never choose it based on agent judgment. Input BED semantics are 0-based half-open. Gene lists accept one identifier per line with explicitly documented header handling; symbols and Ensembl IDs are declared, not guessed. Reserve annotation names such as SNP and reject duplicate/unsafe names. Stage files under input/ with immutable hashes, not arbitrary client paths carried into the worker.

A named scientific profile pins LD window (initially the recovered 1 cM), input normalization, mapping ambiguity/loss policies and zero-coverage rules. Initial proposal: exact current symbols or unique aliases; uniquely resolvable IDs only; no fuzzy corrections. Ambiguous genes fail with candidates; unmapped genes and lost intervals fail unless the request explicitly permits dropping them. Always retain a complete audit table. These strict defaults can be relaxed in a named lab profile without changing worker code.

Execution profile contains SSH alias, executor, runtime, image digest, shared cache/work roots, CPU/memory/time limits and Slurm account/partition. Site paths/account are not scientific provenance or hardcoded analysis defaults. The launcher stages inputs, runs preflight, reserves a run directory and submits the worker. It does not silently choose hosts or fall back to a login node. New/resumed submission records durable IDs; do not blindly retry sbatch after an ambiguous SSH disconnect—reconcile the run token with scheduler records first.

## 6. Worker stages and scientific invariants

### A. Validate and initialize

Validate request/schema, readable staged inputs, available reference, memory/concurrency configuration, unique annotation IDs and supported species/build. Create resolved request and input hashes. Enforce autosomes 1–22 and normalize only recognized chromosome aliases. Unsupported contigs are reported; dropping them must follow the explicit policy.

### B. Normalize genes or intervals

Genes: map IDs against the pinned lookup, strip Ensembl version suffixes only by documented policy, deduplicate, resolve a pinned gene-span definition, apply symmetric flank, and clip to chromosome bounds. Convert 1-based inclusive source coordinates to BED deliberately. Do not interpret gene symbols as genomic coordinates.

BED: parse at least chr/start/end with declared handling of extra columns; reject non-integers, negative starts and end <= start. Preserve input record IDs. For GRCh38 input and GRCh37 reference, apply the pinned hg38->hg19 chain and explicit liftOver parameters. Reject multiple/split mappings unless a future profile explicitly supports them; retain unlifted records and reasons. Record mapped interval/base-pair loss, not just row counts. No silent build guessing. Same-build input skips liftover.

Both routes produce binary union intervals per annotation. Merge overlapping intervals within an annotation so overlapping peaks do not become unintended counts. Do not merge different annotation sets. Emit input-to-output mapping and canonical BED files.

### C. Build SNP annotation and calculate per chromosome

For each annotation and chromosome, generate binary membership against the pinned full LD-reference SNP order. Reuse make_annot.py where verified against the pinned commit. Match thin versus full annotation format explicitly: a thin annotation is only annotation columns; keep its ordered SNP sidecar and require exact agreement with the filtered genotype rows LDSC consumes. Coordinate intersections must account for BIM positions being 1-based.

Run the pinned equivalent of the recovered calculation:

```text
ldsc.py --bfile <chr-prefix> --ld-wind-cm 1
        --annot <thin-annotation> --thin-annot
        --out <task-prefix> --print-snps <pinned-output-list>
```

Validate the committed command's exact behavior against the chosen LDSC commit. Calculate LD against the full eligible panel; the print list restricts output, not the LD neighborhood. Preserve `.annot.gz`, `.l2.ldscore.gz`, `.l2.M`, `.l2.M_5_50` and logs. Explicitly check which genotype variants LDSC excludes (e.g. monomorphic SNPs) so membership/count/order contracts reflect the actual reference universe.

Zero annotation coverage on an individual chromosome is valid and must still yield complete zero annotation/LD-score outputs and counts; handle/validate the pinned tool's behavior. Genome-wide zero coverage is an invalid annotation and fails. Do not reject negative estimated LD scores merely for being negative; finite numerical results and method-specific checks are required.

Each task writes to a temporary directory, validates output, then publishes a completion record with input and output hashes. Limit task concurrency by allocated resources, not detected host CPU count; avoid BLAS oversubscription. Start with one annotation per task for simplicity; multi-column batching is a later optimization.

### D. Convert to ldsR and validate

Write a strict converter rather than blindly invoking to_celltype_dataset: the inspected helper combines columns positionally and can drop incomplete sets. Parse chromosomes numerically 1–22; never lexically sort 1,10,11,...,2. Preserve requested annotation column order and validate SNP identity across every task before combining.

Output directory:

| File | Contract |
|---|---|
| ld.parquet | SNP string + one numeric LD-score column per annotation; regression-output SNP rows in canonical genomic order |
| annot.parquet | annot string, m numeric, m50 numeric; exactly one row per annotation in LD-column order |
| annot_ref.parquet | SNP string + binary columns for every annotation; full eligible reference-SNP universe, not only HM3 output rows |
| snp_freq.parquet | SNP and MAF for that exact reference universe/order, if required by the selected compatibility interface |
| manifest.json | Schema, source/build/reference/software hashes, distinct SNP-order hashes, annotation mapping and LD parameters, output checksums |
| qc.json / qc.tsv | Mapping losses, intervals/base pairs, annotated SNP counts, warnings, all task results |

Sum M/M_5_50 across chromosomes, verify against membership and pinned frequency definition, and assert m50 <= m. Do not recompute these from regression-output SNPs. Verify no duplicate/missing SNPs, finite LD scores, matching annotation columns, complete 22-chromosome outputs and correct common-SNP mask alignment.

The ldsR compatibility gate runs parse_parquet_dir(read_ref=TRUE) inside the image plus stricter identity/order validation against the pinned baseline and bundled mask. Passing the parser alone does not prove overlap compatibility. If the existing baseline/mask cannot be matched, fail this gate and document a targeted ldsR/reference fix; never publish a misleading compatible=true flag.

Publish final outputs atomically only after all annotations pass. Partial raw work remains resumable but is not presented as a complete ldsR dataset. Write COMPLETE.json last with manifest checksum.

## 7. Run layout, cache and recovery

```text
<site-work>/runs/<run-id>/
  input/                 # exact submitted files and request
  resolved-request.json
  state.json
  provenance/
  normalized/            # BED files and mapping tables
  tasks/<annotation>/chr<N>/
  logs/
  output/ldsR/
  COMPLETE.json
<site-cache>/
  images/<digest>/
  references/<bundle-id>/
```

Scientific cache key: normalized annotation hash + input transformation provenance + reference identity + scientific parameters + worker/tool image digest. Executor, host and Slurm job ID do not alter the scientific key. Resume validates hashes and only repeats invalid/incomplete stages; changed input/reference must not reuse stale output. A run-level lock prevents simultaneous resume. Stop/cancellation leaves completed tasks intact. State exposes initializing, staged, queued, running, failed, cancelled, completed with timestamps, task counts, log locations and error codes.

Use explicit error classes such as INPUT_INVALID, BUILD_UNSUPPORTED, MAPPING_LOSS, REFERENCE_MISMATCH, LDSC_TASK_FAILED, RESOURCE_LIMIT, EXPORT_INVALID. Never repair failure by silently changing flanks, reference or SNP filters. Return nonzero exit status for incomplete final output.

## 8. Suggested repository modules

```text
container/Dockerfile
container/software-lock.*
.github/workflows/container.yml
schemas/request.schema.json
schemas/output-manifest.schema.json
references/references.lock.json
profiles/local.example.*
profiles/slurm.example.*
launcher/                  # runtime, transport, scheduler; no science
worker/                    # CLI, state, DAG, bounded processes, cache
worker/annotations/        # genes, BED, liftover, audit tables
worker/ldscores/            # per-chromosome tasks and checks
R/export_ldsR.R
R/validate_ldsR.R
tests/fixtures/
tests/integration/
docs/installation.md
docs/input-output-contract.md
```

Implementation can adjust filenames, but keep scheduler/transport separate from science and conversion. Use argument vectors, not shell interpolation of filenames. R scripts follow the user's tidyverse/native-pipe conventions. Extract only reviewed logic from historical source snapshots; do not inherit their absolute paths, floating image tags or polling-based dependencies.

## 9. Implementation sequence and acceptance

1. **Lock and build spike:** audit Zenodo archive, establish coordinate/mask identity, choose exact gene release and tool commits, create reference/software locks, build image in Actions and pull it through Docker and Apptainer. Resolve legacy build feasibility before progressing.
2. **BED vertical slice:** human GRCh37 BED -> complete raw LD outputs -> strict ldsR export, direct executor. Include a tiny transparent genotype fixture with independently calculated expected interval membership and LD-score checks.
3. **Human mappings:** gene IDs/symbols and GRCh38 liftover, explicit flanks/loss policies, boundary and chain fixtures; equivalent genes/BED inputs should generate the same canonical annotation and scores.
4. **Execution:** launcher local profile, then SSH + Slurm single-allocation profile; offline worker, bounded resources, durable job state, logs, resume and result retrieval. Dardel job submission follows the user's per-job authorization workflow during development.
5. **Release validation:** full 22-chromosome reference comparison against the historical/pinned LDSC command; Docker/Apptainer agreement; ldsR reader and baseline/mask identity checks; publish digest and acceptance evidence.
6. **Later:** mouse input adapters, chromosome arrays if measured demand warrants them, continuous annotations, native ARM.

Required tests: BED off-by-one boundaries; merged overlapping peaks; gene dedup/ambiguous aliases; declared-build mismatch; unmapped and multiply mapped intervals; corrupted reference; chromosome 2 versus 10 ordering; shuffled SNP rows; incomplete chromosome; monomorphic genotype filtering; zero coverage on one chromosome; global zero coverage; multi-annotation column/count consistency; interrupted task; concurrent cache installation; duplicate submission recovery; changed flank/reference invalidation; stale COMPLETE marker; failed task returns nonzero and does not publish completed export.

Discrete mapping/membership/count results must match exactly. Provisional floating comparison tolerance: abs(a-b) <= 1e-8 + 1e-6*abs(b), with observed maximum errors recorded and any adjustment justified before release. Compare semantic parquet data, not byte identity, across runtimes. Independent small expected-value fixtures are necessary because simply running the same implementation twice is not validation.

Definition of done: a fresh user with only a supported container runtime and configured execution profile can submit human genes or BED input, receive a complete validated ldsR directory, and reproduce it from the manifest. A failed/interrupted run can resume without repeating successful chromosome work. No host scientific package installation, no network access from compute tasks, no guessed build and no GWAS dependency.
