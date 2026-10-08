# Installation and operation

## Requirements

- Execution site: Apptainer (≥ 1.3) or Docker on linux/amd64, plus Slurm client tools if
  `executor.type` is `slurm`. No R, Python, conda, bedtools or liftOver on the host.
- Client (where you type `agent-ldsc`): Python ≥ 3.10 (stdlib only), `ssh` and `rsync`.

```sh
uv tool install git+https://github.com/Ararder/agent-ldsc   # or: pip install --user .
mkdir -p ~/.config/agent-ldsc/profiles
cp profiles/execution/dardel.example.json ~/.config/agent-ldsc/profiles/dardel.json   # edit image digest
```

## Image

Built by `.github/workflows/container.yml`. Every push to `main` publishes
`ghcr.io/ararder/agent-ldsc:sha-<commit>` together with SBOM and provenance attestations, then
tests that exact digest under Docker and Apptainer and checks parity. A `v*` tag re-tags the
tested digest, with no rebuild. Profiles must reference the image by digest.

On Apptainer sites, `agent-ldsc stage` converts the digest once to
`<cache_root>/images/<digest>.sif`. It records the SIF sha256 next to the file and checks it
before every use. Apptainer cache/tmp go under `<cache_root>`, never `$HOME`.

Software pins: `container/locks/*.lock` (explicit conda URLs + md5),
`container/envs/ldsc-requirements.txt` (PyPI wheels by sha256), bulik/ldsc `aa33296`,
ldsR `e880755`, micromamba 2.9.0 (sha256), Debian base by digest.

## Reference bundle

`references/references.lock.json` pins every source (URL, size, md5 where published, sha256).
`agent-ldsc stage --profile P` runs `agent-ldsc-worker refs install` inside the image on the
site's network-enabled host. That step downloads only missing sources, verifies them, derives
the SNP indexes and gene tables, checks identity against ldsR's common-SNP mask and the
baseline-v1.2 deposit, and atomically publishes `<cache_root>/references/<bundle-id>/` with
`BUNDLE.json`. Each run re-verifies every file hash before use. Compute tasks never download.

Offline install: copy the source files to the site and run
`agent-ldsc-worker refs install --bundle eur-phase3-grch37-v1 --cache /cache --source-dir /sources --offline`.

## Running

```sh
agent-ldsc submit request.json --profile dardel     # stages image/bundle if missing, then sbatch
agent-ldsc status RUN_ID [--json]
agent-ldsc resume RUN_ID                            # refuses while a job for the run is active
agent-ldsc cancel RUN_ID
agent-ldsc fetch RUN_ID --destination results/      # only complete runs; verifies checksums
```

Measured on the full 1000G panel: about 0.5 GB RSS per concurrent chromosome task. On Dardel's
`shared` partition, memory is tied to cores (about 0.8 GB per CPU), so a larger `--mem` silently
raises the CPU allocation, and with it billing and concurrency. Keep `mem` ≈ 0.75 GB × `cpus`.

One Slurm allocation per run. Inside it, the worker runs up to `--cpus-per-task` chromosome
tasks at once, each with single-threaded BLAS. A wall-time kill leaves completed tasks in
place, and `resume` recomputes only invalid or missing ones.

## Licensing note

The image redistributes UCSC `liftOver` (bioconda `ucsc-liftover`). UCSC Genome Browser
binaries are free for academic, non-profit and personal use; commercial use requires a license
from UCSC. bulik/ldsc is GPL-3.0 and ldsR is MIT.

agent-ldsc itself is MIT-licensed (LICENSE). The image additionally contains GPL-3.0 bulik/ldsc
and UCSC liftOver under its own terms.
