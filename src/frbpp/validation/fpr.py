"""False-positive-rate battery over synthetic nulls."""
from __future__ import annotations

import numpy as np

from frbpp.pointprocess.decluster import decluster
from frbpp.pointprocess.null import (empirical_p, max_stat_null,
                                     sample_uniform_given_exposure,
                                     scramble_episodes)


def fpr_battery_n1(scanner, grid_t, grid_w, n_events, n_real, n_null, alpha, rng,
                   null_samples=None):
    if null_samples is None:
        null_samples = max_stat_null(scanner, grid_t, grid_w, n_events, n_null, rng)
    lams = np.empty(int(n_real), dtype=float)
    for i in range(int(n_real)):
        ev = sample_uniform_given_exposure(grid_t, grid_w, n_events, rng)
        lam, _ = scanner.scan(ev)
        lams[i] = float(lam.max())
    p = np.array([empirical_p(l, null_samples) for l in lams], dtype=float)
    n_false = int(np.sum(p < alpha))
    return {"null_samples": np.asarray(null_samples), "lambdas": lams, "p_values": p,
            "n_false": n_false, "fpr": n_false / int(n_real),
            "n_real": int(n_real), "alpha": float(alpha),
            "null_median": float(np.median(null_samples)),
            "null_q99": float(np.quantile(null_samples, 0.99))}


def fpr_battery_n2(scanner, grid_t, grid_w, n_events, n_real, n_surr, alpha, rng,
                   t_cluster_days: float = 1.0):
    ps = np.empty(int(n_real), dtype=float)
    for i in range(int(n_real)):
        ev = sample_uniform_given_exposure(grid_t, grid_w, n_events, rng)
        lam, _ = scanner.scan(ev)
        lam_hat = float(lam.max())
        p_hat = float(scanner.p_grid[int(np.argmax(lam))])
        ep, _ = decluster(ev, t_cluster_days=t_cluster_days)
        scr = np.empty(int(n_surr), dtype=float)
        for j in range(int(n_surr)):
            lam_s, _ = scanner.scan(scramble_episodes(ep, p_hat, rng))
            scr[j] = float(lam_s.max())
        ps[i] = empirical_p(lam_hat, scr)
    n_false = int(np.sum(ps < alpha))
    return {"p_values": ps, "n_false": n_false, "fpr": n_false / int(n_real),
            "n_real": int(n_real), "alpha": float(alpha), "n_surrogate": int(n_surr),
            "null_mode": "full_scan_max"}