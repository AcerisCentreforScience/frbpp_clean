"""Tier-1 model ladder and Tier-2 grammar enumeration under anchor constraints."""
from __future__ import annotations

import numpy as np

WIDTHS_V2 = np.array([0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.40])
PHI0_24 = np.arange(24) / 24.0
_W_MAX_ANCHOR = 0.5
_DBIC_MIN = 6.0


def _phases(t, t0, P):
    return np.mod((np.asarray(t, dtype=float) - t0) / P, 1.0)


def fit_m0(n: int, t_tot: float):
    ll = 0.0 if n == 0 else float(n * np.log(n / t_tot))
    return ll - n, 1, {}


def fit_m1(ep, t0, P, grid_t, grid_w, phi0_grid=PHI0_24, widths=WIDTHS_V2):
    ph_e = _phases(ep, t0, P)
    ph_g = _phases(grid_t, t0, P)
    T_tot = float(grid_w.sum())
    n = np.asarray(ep).size
    widths = np.asarray([min(float(w), _W_MAX_ANCHOR) for w in widths])
    best_ll, best = -np.inf, {}
    for p0 in phi0_grid:
        de = np.mod(ph_e - p0, 1.0)
        dg = np.mod(ph_g - p0, 1.0)
        for w in widths:
            n_in = int((de < w).sum())
            T_in = float(grid_w[dg < w].sum())
            n_out, T_out = n - n_in, T_tot - T_in
            if (n_in > 0 and T_in <= 0) or (n_out > 0 and T_out <= 0):
                continue
            ll = 0.0 if n_in == 0 else n_in * np.log(n_in / T_in)
            ll += 0.0 if n_out == 0 else n_out * np.log(n_out / T_out)
            ll -= n
            if ll > best_ll:
                best_ll, best = ll, {"phi0": float(p0), "width": float(w),
                                     "n_in": n_in, "T_in": T_in}
    return float(best_ll), 4, best


def fit_m2(ep, t0, P, grid_t, grid_w, n_phi=48, n_a=40, a_max=0.95):
    ph_e = _phases(ep, t0, P)
    ph_g = _phases(grid_t, t0, P)
    th_e, th_g = 2 * np.pi * ph_e, 2 * np.pi * ph_g
    T_tot = float(grid_w.sum())
    n = np.asarray(ep).size
    best_ll, best = -np.inf, {}
    for phi in np.linspace(0.0, 2 * np.pi, n_phi, endpoint=False):
        se = np.sin(th_e + phi)
        Sg = float((grid_w * np.sin(th_g + phi)).sum())
        for A in np.linspace(0.05, a_max, n_a):
            m = 1.0 + A * se
            if (m <= 0).any():
                continue
            I = T_tot + A * Sg
            if I <= 0:
                continue
            ll = n * np.log(n / I) + float(np.log(m).sum()) - n
            if ll > best_ll:
                best_ll, best = ll, {"phi": float(phi), "A": float(A)}
    return float(best_ll), 3, best


def _mask_terms(ep, t0, P, grid_t, grid_w, phi0, width):
    ph_e = _phases(ep, t0, P)
    ph_g = _phases(grid_t, t0, P)
    me = np.mod(ph_e - phi0, 1.0) < width
    mg = np.mod(ph_g - phi0, 1.0) < width
    return ph_e, ph_g, me, mg


def fit_m3(ep, t0, P, grid_t, grid_w, phi0, width, t_mid, half_span):
    from scipy.optimize import minimize

    ph_e, ph_g, me, mg = _mask_terms(ep, t0, P, grid_t, grid_w, phi0, width)
    n_in = int(me.sum())
    if n_in < 3:
        return -np.inf, 5, {}
    T_tot = float(grid_w.sum())
    T_in = float(grid_w[mg].sum())
    T_out = T_tot - T_in
    tau_e = (np.asarray(ep)[me] - t_mid) / half_span
    tau_g = (grid_t[mg] - t_mid) / half_span
    Tg = float((grid_w[mg] * tau_g).sum())
    n_out = np.asarray(ep).size - n_in

    def nll(x):
        a, beta, c = x
        if a <= 0 or c <= 0 or abs(beta) >= 0.99:
            return 1e18
        if (1.0 + beta * tau_e <= 0).any():
            return 1e18
        ll = n_in * np.log(a) + float(np.log(1.0 + beta * tau_e).sum())
        ll -= a * (T_in + beta * Tg)
        if n_out:
            ll += n_out * np.log(c)
        ll -= c * T_out
        return -ll

    a0 = n_in / T_in if T_in > 0 else 1e-3
    c0 = n_out / T_out if T_out > 0 else 1e-3
    r = minimize(nll, [max(a0, 1e-3), 0.0, max(c0, 1e-3)], method="L-BFGS-B",
                 bounds=[(1e-9, None), (-0.99, 0.99), (1e-9, None)])
    a_f, beta_f, c_f = (float(v) for v in r.x)
    return -float(r.fun), 5, {"a": a_f, "beta": beta_f, "c": c_f,
                              "b_per_day": a_f * beta_f / max(half_span, 1e-9)}


def fit_s1(ep, t0, P, grid_t, grid_w, phi0, width, n_phi=24, n_b=20, b_max=0.95):
    ph_e, ph_g, me, mg = _mask_terms(ep, t0, P, grid_t, grid_w, phi0, width)
    n_in = int(me.sum())
    if n_in < 3:
        return -np.inf, 6, {}
    T_tot = float(grid_w.sum())
    T_in = float(grid_w[mg].sum())
    T_out = T_tot - T_in
    th_e, th_g = 2 * np.pi * ph_e[me], 2 * np.pi * ph_g[mg]
    n_out = np.asarray(ep).size - n_in
    best_ll, best = -np.inf, {}
    for phis in np.linspace(0.0, 2 * np.pi, n_phi, endpoint=False):
        se = np.sin(th_e + phis)
        Sg = float((grid_w[mg] * np.sin(th_g + phis)).sum())
        for B in np.linspace(0.05, b_max, n_b):
            m = 1.0 + B * se
            if (m <= 0).any():
                continue
            I_in = T_in + B * Sg
            if I_in <= 0:
                continue
            a_in = n_in / I_in
            c_out = (n_out / T_out) if (n_out and T_out > 0) else 0.0
            ll = n_in * np.log(a_in) + float(np.log(m).sum()) - a_in * I_in
            if n_out and c_out > 0:
                ll += n_out * np.log(c_out) - c_out * T_out
            if ll > best_ll:
                best_ll, best = ll, {"phis": float(phis), "B": float(B)}
    return float(best_ll), 6, best


def ladder_fit(ep, t0, P, grid_t, grid_w):
    from frbpp.pointprocess.periodogram import _check_epoch

    ep = np.asarray(ep, dtype=float)
    _check_epoch(ep, np.asarray(grid_t, dtype=float))
    n = ep.size
    t_mid = float((ep.min() + ep.max()) / 2.0)
    half_span = max(float(ep.max() - ep.min()) / 2.0, 1e-6)
    T_tot = float(grid_w.sum())
    ll0, k0, _ = fit_m0(n, T_tot)
    ll1, k1, d1 = fit_m1(ep, t0, P, grid_t, grid_w)
    ll2, k2, d2 = fit_m2(ep, t0, P, grid_t, grid_w)
    if np.isfinite(ll1) and d1:
        ll3, k3, d3 = fit_m3(ep, t0, P, grid_t, grid_w, d1["phi0"], d1["width"],
                             t_mid, half_span)
        ll4, k4, d4 = fit_s1(ep, t0, P, grid_t, grid_w, d1["phi0"], d1["width"])
    else:
        ll3, k3, d3, ll4, k4, d4 = (-np.inf, 5, {}), (-np.inf, 6, {})
    ln_n = float(np.log(n))
    table = {
        "M0": {"logL": ll0, "k": k0, "bic": -2 * ll0 + k0 * ln_n, "det": {}},
        "M1": {"logL": ll1, "k": k1, "bic": -2 * ll1 + k1 * ln_n, "det": d1},
        "M2": {"logL": ll2, "k": k2, "bic": -2 * ll2 + k2 * ln_n, "det": d2},
        "M3": {"logL": ll3, "k": k3, "bic": -2 * ll3 + k3 * ln_n, "det": d3},
        "S1": {"logL": ll4, "k": k4, "bic": -2 * ll4 + k4 * ln_n, "det": d4},
    }
    bic0 = table["M0"]["bic"]
    sup = {}
    sup["M1"] = bool(np.isfinite(table["M1"]["logL"]) and bic0 - table["M1"]["bic"] >= _DBIC_MIN)
    sup["M2"] = bool(np.isfinite(table["M2"]["logL"]) and bic0 - table["M2"]["bic"] >= _DBIC_MIN)
    sup["M3"] = bool(sup["M1"] and np.isfinite(table["M3"]["logL"]) and
                     table["M1"]["bic"] - table["M3"]["bic"] >= _DBIC_MIN)
    sup["S1"] = bool(sup["M1"] and np.isfinite(table["S1"]["logL"]) and
                     table["M1"]["bic"] - table["S1"]["bic"] >= _DBIC_MIN)
    supported = [m for m in ("M1", "M2", "M3", "S1") if sup[m]]
    final = min(supported, key=lambda m: table[m]["bic"]) if supported else "M0"
    degenerate = ("M1" in supported and "M2" in supported and
                  abs(table["M1"]["bic"] - table["M2"]["bic"]) < _DBIC_MIN)
    return {"final": final, "supported": supported, "degenerate": bool(degenerate),
            "table": table, "n_episodes": int(n), "P": float(P)}