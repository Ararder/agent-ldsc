# agent-ldsc

Human genes or genomic intervals (e.g. ATAC peaks) → per-chromosome S-LDSC LD scores → a
validated [ldsR](https://github.com/Ararder/ldsR) parquet directory, with full provenance.
All science runs in one pinned container image; a small stdlib launcher stages, submits and
fetches runs on a local machine or over SSH (Slurm or direct).

- Scope and design: [LDSCORE-IMPLEMENTATION-PLAN.md](LDSCORE-IMPLEMENTATION-PLAN.md)
- Request/output contract, error codes, tolerances: [docs/input-output-contract.md](docs/input-output-contract.md)
- Installation, image, reference bundle, operation: [docs/installation.md](docs/installation.md)
- Recovered historical method and reference audit: [docs/provenance.md](docs/provenance.md)

```sh
agent-ldsc submit examples/atac-peaks.request.json --profile dardel
agent-ldsc status <run-id>
agent-ldsc fetch <run-id> --destination results/
```

Layout: `container/` (Dockerfile, locks, CI helpers), `worker/agent_ldsc_worker/` (in-image
pipeline), `launcher/agent_ldsc_launcher/` (host side, no science), `R/validate_ldsR.R`
(ldsR gate), `references/references.lock.json`, `profiles/{scientific,execution}/`,
`schemas/`, `tests/` (synthetic 22-chromosome fixture + independent numpy oracle).

Local development: `uv venv && uv pip install -e '.[dev]' && pytest` (tests needing ldsc,
liftOver or ldsR skip unless available; CI runs all of them inside the image).
