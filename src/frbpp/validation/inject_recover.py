"""Injection-recovery completeness surface."""
from __future__ import annotations

import numpy as np

from frbpp.pointprocess.null import sample_periodic_signal


def sample_gaussian_window_signal(grid_t, grid_w, t0, P, phi_c, sigma_frac,
                                  contrast, n_target, rng):
    t = np.asarray(grid_t, dtype=float)
    w = np.asarray(grid_w, dtype=float)
    accepted, got = [], 0
    while got < n_target:
        idx = rng.choice(len(t), size=max(2 * (n_target - got), 64), replace=True,
                         p=w / w.sum())
        cand = t[idx]
        ph = np.mod((cand - t0) / P, 1.0)
        dphi = np.mod(ph - phi_c + 0.5, 1.0) - 0.5
        g = np.exp(-0.5 * (dphi / sigma_frac) ** 2)
        accept_p = (1.0 + (contrast - 1.0) * g) / contrast
        keep = cand[rng.random(cand.size) < accept_p]
        accepted.append(keep)
        got += keep.size
    out = np.concatenate(accepted)
    if out.size > n_target:
        out = rng.choice(out, size=n_target, replace=False)
    return np.sort(out)


def recovery_trial(scanner, grid_t, grid_w, t0, P_true, width, contrast, n_ep,
                   rng, gaussian_sigma=None):
    phi0 = float(rng.uniform(0.0, 1.0))
    if gaussian_sigma is None:
        ev = sample_periodic_signal(grid_t, grid_w, t0, P_true, phi0, width,
                                    contrast, n_ep, rng)
    else:
        ev = sample_gaussian_window_signal(grid_t, grid_w, t0, P_true, phi0,
                                           gaussian_sigma, contrast, n_ep, rng)
    lam, _ = scanner.scan(ev)
    p_hat = float(scanner.p_grid[int(np.argmax(lam))])
    return abs(p_hat - P_true) / P_true


def completeness_surface(scanner, grid_t, grid_w, t0, configs, n_rep,
                         contrast, tol_frac, rng):
    from scipy.stats import beta

    rows = []
    for cfg in configs:
        hits = 0
        for _ in range(int(n_rep)):
            dev = recovery_trial(scanner, grid_t, grid_w, t0, cfg["P"],
                                 cfg["width"], contrast, cfg["n_ep"], rng,
                                 gaussian_sigma=cfg.get("gaussian_sigma"))
            hits += int(dev <= tol_frac)
        lo, hi = beta.ppf([0.025, 0.975], hits + 0.5, n_rep - hits + 0.5)
        rows.append({"P": cfg["P"], "width": cfg["width"], "n_ep": cfg["n_ep"],
                     "gaussian_sigma": cfg.get("gaussian_sigma"),
                     "hits": hits, "n_rep": int(n_rep),
                     "recovery": hits / int(n_rep),
                     "ci95": [float(lo), float(hi)]})
    return rows