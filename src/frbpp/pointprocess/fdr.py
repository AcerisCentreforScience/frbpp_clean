"""Benjamini-Hochberg FDR across sources."""
from __future__ import annotations

import numpy as np


def bh_fdr(pvalues):
    p = np.asarray(pvalues, dtype=float)
    n = p.size
    q = np.empty(n, dtype=float)
    order = np.argsort(p)
    prev = 1.0
    for rank in range(n, 0, -1):
        i = order[rank - 1]
        prev = min(prev, p[i] * n / rank)
        q[i] = prev
    return q