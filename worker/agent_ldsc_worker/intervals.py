"""Interval primitives. Intervals are BED-style: 0-based, half-open [start, end)."""

from __future__ import annotations

import numpy as np
import pandas as pd

from .util import AUTOSOMES


def merge_intervals(df: pd.DataFrame) -> pd.DataFrame:
    """Union of intervals per chromosome (overlapping or book-ended intervals are joined).

    Expects columns chrom (int 1-22), start, end. Returns sorted chrom/start/end.
    """
    out = []
    for chrom, g in df.sort_values(["chrom", "start", "end"]).groupby("chrom", sort=True):
        s = g.start.to_numpy(np.int64)
        e = g.end.to_numpy(np.int64)
        run_end = np.maximum.accumulate(e)
        new = np.ones(len(s), dtype=bool)
        new[1:] = s[1:] > run_end[:-1]
        idx = np.flatnonzero(new)
        ends = np.append(idx[1:], len(s)) - 1
        out.append(pd.DataFrame({"chrom": chrom, "start": s[idx], "end": run_end[ends]}))
    if not out:
        return pd.DataFrame({"chrom": pd.Series(dtype=int), "start": pd.Series(dtype=np.int64),
                             "end": pd.Series(dtype=np.int64)})
    return pd.concat(out, ignore_index=True)


def total_bp(df: pd.DataFrame) -> int:
    return int((df.end - df.start).sum()) if len(df) else 0


def snp_membership(merged: pd.DataFrame, chrom: int, bp: np.ndarray) -> np.ndarray:
    """Binary membership of SNPs (1-based BIM positions) in merged intervals on one chromosome.

    A SNP at 1-based position p occupies the BED base [p-1, p); it is inside [s, e) iff
    s <= p-1 < e. This is exactly make_annot.py's bedtools intersect of [BP-1, BP).
    """
    g = merged[merged.chrom == chrom]
    if len(g) == 0:
        return np.zeros(len(bp), dtype=np.int8)
    s = g.start.to_numpy(np.int64)
    e = g.end.to_numpy(np.int64)
    pos0 = np.asarray(bp, dtype=np.int64) - 1
    i = np.searchsorted(s, pos0, side="right") - 1
    inside = (i >= 0) & (pos0 < e[np.clip(i, 0, None)])
    return inside.astype(np.int8)


def per_chrom_summary(merged: pd.DataFrame) -> dict:
    return {str(c): {"intervals": int((merged.chrom == c).sum()),
                     "bp": total_bp(merged[merged.chrom == c])} for c in AUTOSOMES}
