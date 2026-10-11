"""Aliasing and phase-degeneracy guards for period searches."""
from __future__ import annotations

import numpy as np

F_LO, F_HI, N_FREQ = 1.0 / 200.0, 1.0, 4000
RANK_FLAG = 0.995
W_REF = 0.05
N_BINS = 200


def capture_ratio(grid_t, grid_w, p_grid, w_ref=W_REF, n_bins=N_BINS):
    t = np.asarray(grid_t, dtype=float)
    w = np.asarray(grid_w, dtype=float)
    tot = float(w.sum())
    nb = int(n_bins)
    win = max(int(round(w_ref * nb)), 1)
    out = np.empty(len(p_grid), dtype=float)
    for i, P in enumerate(p_grid):
        b = np.minimum((np.mod(t / P, 1.0) * nb).astype(int), nb - 1)
        hist = np.bincount(b, weights=w, minlength=nb)
        full = np.concatenate([hist, hist[:win]])
        csum = np.concatenate([[0.0], np.cumsum(full)])
        rolling = csum[win:win + nb] - csum[:nb]
        out[i] = float(rolling.max()) / tot
    return out


def alias_test(grid_t, grid_w, p_hat, n_freq=N_FREQ, f_lo=F_LO, f_hi=F_HI,
               w_ref=W_REF):
    f_grid = np.linspace(f_lo, f_hi, n_freq)
    p_grid = 1.0 / f_grid
    rho = capture_ratio(grid_t, grid_w, p_grid, w_ref=w_ref)
    rho_hat = float(capture_ratio(grid_t, grid_w,
                                  np.asarray([float(p_hat)], dtype=float),
                                  w_ref=w_ref)[0])
    rank = float((rho < rho_hat).mean())
    return {"rank": rank, "rho_hat": rho_hat,
            "rho_q995": float(np.quantile(rho, 0.995)),
            "flagged": bool(rank >= RANK_FLAG)}