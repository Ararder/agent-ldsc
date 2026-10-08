"""Independent LD-score oracle for the fixture bundle (no ldsc code).

For the fixture's chunk-aligned clusters, ldsc's L2 for SNP j equals
    sum over SNPs k in j's cluster of  a_k * (r_jk^2 - (1 - r_jk^2) / (n - 2))
with r computed from genotypes standardized to mean 0 and population SD 1.
"""

from __future__ import annotations

import numpy as np


def read_plink_bed(path, n_ind: int, n_snp: int) -> np.ndarray:
    raw = np.fromfile(path, dtype=np.uint8)
    assert raw[:3].tolist() == [0x6C, 0x1B, 0x01]
    per = (n_ind + 3) // 4
    blocks = raw[3:].reshape(n_snp, per)
    codes = np.stack([(blocks >> (2 * k)) & 0b11 for k in range(4)], axis=2).reshape(n_snp, per * 4)[:, :n_ind]
    dosage = np.select([codes == 0, codes == 2, codes == 3], [2.0, 1.0, 0.0], default=np.nan)
    assert not np.isnan(dosage).any()
    return dosage.T  # n_ind x n_snp


def ld_scores(geno: np.ndarray, cluster_ids: np.ndarray, annot: np.ndarray) -> np.ndarray:
    n = geno.shape[0]
    x = (geno - geno.mean(axis=0)) / geno.std(axis=0)
    out = np.zeros((geno.shape[1], annot.shape[1]))
    for cl in np.unique(cluster_ids):
        idx = np.flatnonzero(cluster_ids == cl)
        r = x[:, idx].T @ x[:, idx] / n
        r2 = r ** 2 - (1 - r ** 2) / (n - 2)
        out[idx] = r2 @ annot[idx]
    return out
