# agent-ldsc

LD-score component of LDSCORE-IMPLEMENTATION-PLAN.md is implemented (worker, launcher, image, reference bundle); see README.md and docs/. The broader design draft for the not-yet-implemented regression part is local-only in notes/DESIGN-DRAFT.md (with notes/DARDEL-DISCOVERY.md).

Remote hosts: Dardel (SSH alias `dardel`) and the lab server (SSH alias `lab`). Read the matching skill before operating there. Proposed site roots (cache + runs): Dardel `/cfs/klemming/projects/supr/ki-pgi-storage/shared/arvhar/agent-ldsc/`, lab `/mnt/data/user/arvid/agent-ldsc/`. Existing source projects on Dardel: `/cfs/klemming/projects/supr/ki-pgi-storage/shared/arvhar/generate-ldscores` and `stratified-ldscore-regression`. Shared inputs: `/cfs/klemming/projects/supr/ki-pgi-storage/Data`. Historical ldsR outputs used as comparison targets: `Data/ldsR_ldscores/` (e.g. `ss3x_fixed500_19`).

Local `data/` and `source-snapshots/` are ignored. Keep lab GWAS files out of Git. Source snapshots are historical evidence, not new implementation. Method provenance: docs/provenance.md. Local GWAS copies and their manifest: data/lab-gwas/.
