"""Surrogate null generators for the point-process tests."""
from __future__ import annotations

import numpy as np


def sample_uniform_given_exposure(grid_t, grid_w, n, rng):
    w = np.asarray(grid_w, dtype=float)
    p = w / w.sum()
    idx = rng.choice(len(grid_t), size=int(n), replace=False, p=p)
    return np.sort(np.asarray(grid_t, dtype=float)[idx])


def sample_periodic_signal(grid_t, grid_w, t0, P, phi0, width_frac, contrast,
                           n_target, rng):
    t = np.asarray(grid_t, dtype=float)
    w = np.asarray(grid_w, dtype=float)
    contrast = float(contrast)
    if contrast < 1.0:
        raise ValueError("contrast 必须 ≥ 1")
    accepted = []
    got = 0
    while got < n_target:
        idx = rng.choice(len(t), size=max(2 * (n_target - got), 64), replace=True,
                         p=w / w.sum())
        cand = t[idx]
        ph = np.mod((cand - t0) / P, 1.0)
        in_win = np.mod(ph - phi0, 1.0) < width_frac
        accept_p = np.where(in_win, 1.0, 1.0 / contrast)
        keep = cand[rng.random(cand.size) < accept_p]
        accepted.append(keep)
        got += keep.size
    out = np.concatenate(accepted)
    if out.size > n_target:
        out = rng.choice(out, size=n_target, replace=False)
    return np.sort(out)


def scramble_episodes(episode_times, P, rng):
    t = np.asarray(episode_times, dtype=float)
    return t + rng.uniform(0.0, 1.0, size=t.size) * float(P)


def max_stat_null(scanner, grid_t, grid_w, n_events, n_surr, rng):
    lams = np.empty(int(n_surr), dtype=float)
    for i in range(int(n_surr)):
        ev = sample_uniform_given_exposure(grid_t, grid_w, n_events, rng)
        lam, _ = scanner.scan(ev)
        lams[i] = float(lam.max())
    return lams


def empirical_p(observed, null_samples):
    null = np.asarray(null_samples, dtype=float)
    return float((1 + int(np.sum(null >= observed))) / (1 + null.size))