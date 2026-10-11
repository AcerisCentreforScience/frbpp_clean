"""DM(t) analysis: empirical sigma floor, inflation factor, model ladder, cross-catalogue offset."""
from __future__ import annotations

import numpy as np

_TREND_DBIC = 6.0
_TREND_P = 0.05
N_PERM = 300


def sigma_eff_jitter(dm: np.ndarray) -> float:
    d = np.asarray(dm, dtype=float)
    if d.size < 3:
        return float("nan")
    return float(1.4826 * np.median(np.abs(np.diff(d) - np.median(np.diff(d)))) / np.sqrt(2))


def inflation_factor(sig_eff: float, formal: np.ndarray) -> float:
    f = np.asarray(formal, dtype=float)
    med = float(np.median(f))
    return float(sig_eff / med) if (np.isfinite(sig_eff) and med > 0) else float("nan")


def wls_poly(t, y, err, deg: int):
    t = np.asarray(t, dtype=float)
    y = np.asarray(y, dtype=float)
    err = np.asarray(err, dtype=float)
    w = 1.0 / err
    tc = t - t.mean()
    g = np.polyfit(tc, y, deg, w=w)
    yhat = np.polyval(g, tc)
    chi2 = float(np.sum(((y - yhat) / err) ** 2))
    return g, chi2, yhat


def wls_sinusoid(t, y, err, P, t0, n_phi=48):
    t = np.asarray(t, dtype=float)
    th = 2 * np.pi * (t - t0) / P
    A = np.column_stack([np.ones_like(th), np.sin(th), np.cos(th)])
    w = 1.0 / np.asarray(err, dtype=float)
    Aw = A * w[:, None]
    bw = np.asarray(y, dtype=float) * w
    coef, *_ = np.linalg.lstsq(Aw, bw, rcond=None)
    yhat = A @ coef
    chi2 = float(np.sum(((np.asarray(y, dtype=float) - yhat) / np.asarray(err, dtype=float)) ** 2))
    amp = float(np.hypot(coef[1], coef[2]))
    return coef, chi2, amp, yhat


def fit_dm_models(t, dm, err, p_burst=None, rng=None, n_perm=N_PERM):
    t = np.asarray(t, dtype=float)
    dm = np.asarray(dm, dtype=float)
    err = np.asarray(err, dtype=float)
    n = t.size
    out = {"n": int(n), "sigma_jit": np.nan, "inflation": np.nan,
           "floor_applied": None, "scatter_gate_ok": None, "err_scale": None,
           "trend_significant": None, "models": {}, "mc_p_trend": np.nan,
           "sinusoid": None}
    if n < 4:
        out["note"] = "样本过少（<4 暴），不出结论"
        return out
    sig = sigma_eff_jitter(dm)
    out["sigma_jit"] = sig
    out["inflation"] = inflation_factor(sig, err)
    if np.isfinite(sig):
        err_use = np.maximum(err, sig)
        out["err_scale"] = "jitter-floored"
    else:
        err_use = err
        out["err_scale"] = "formal"
    out["floor_applied"] = bool(out["err_scale"] == "jitter-floored")

    for deg, name in ((0, "M0"), (1, "M1"), (2, "M2")):
        if n <= deg + 1:
            continue
        g, chi2, _ = wls_poly(t, dm, err_use, deg)
        bic = chi2 + (deg + 1) * np.log(n)
        out["models"][name] = {"coef_high2low": [float(x) for x in g],
                               "chi2": chi2, "bic": float(bic)}

    mc_p = np.nan
    if "M1" in out["models"] and rng is not None:
        dchi_obs = out["models"]["M0"]["chi2"] - out["models"]["M1"]["chi2"]
        cnt = 0
        for _ in range(int(n_perm)):
            dmp = rng.permutation(dm)
            _, c0, _ = wls_poly(t, dmp, err_use, 0)
            _, c1, _ = wls_poly(t, dmp, err_use, 1)
            cnt += (c0 - c1) >= dchi_obs
        mc_p = (1 + cnt) / (n_perm + 1)
    out["mc_p_trend"] = float(mc_p)

    dbic10 = (out["models"]["M0"]["bic"] - out["models"]["M1"]["bic"]
              if "M1" in out["models"] else np.nan)
    trend = bool(np.isfinite(dbic10) and dbic10 >= _TREND_DBIC and mc_p < _TREND_P)
    out["trend_significant"] = trend
    out["dbic_M0_M1"] = float(dbic10)
    out["scatter_gate_ok"] = bool(n > 2 and np.std(dm, ddof=1) >= sig) if np.isfinite(sig) else None

    if p_burst is not None and n >= 6:
        coef, chi2, amp, _ = wls_sinusoid(t, dm, err, float(p_burst), float(t.min()))
        _, chi2_0u, _ = wls_poly(t, dm, err, 0)
        bic0u = chi2_0u + 1 * np.log(n)
        bic3 = chi2 + 3 * np.log(n)
        dbic30 = float(bic0u - bic3)
        mc_p3 = np.nan
        if rng is not None:
            cnt = 0
            for _ in range(int(n_perm)):
                _, chi2_p, _, _ = wls_sinusoid(t, rng.permutation(dm), err,
                                               float(p_burst), float(t.min()))
                cnt += chi2_p <= chi2
            mc_p3 = (1 + cnt) / (n_perm + 1)
        out["sinusoid"] = {"P": float(p_burst), "amp": amp, "chi2": chi2,
                           "dbic_M0_M3": dbic30, "mc_p": float(mc_p3),
                           "err_class": "coherent(no-floor)",
                           "significant": bool(dbic30 >= _TREND_DBIC and mc_p3 < _TREND_P)}
    return out