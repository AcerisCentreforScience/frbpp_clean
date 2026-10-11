"""Tests for the DM floor, trend guards, model ladder, offset time series."""
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from frbpp.analysis.dm import (fit_dm_models, sigma_eff_jitter,   # noqa: E402
                               wls_poly)


def test_sigma_eff_recovers_known_jitter():
    rng = np.random.default_rng(0)
    jitter = 0.05
    dm = 350 + rng.normal(0, jitter, size=400)
    est = sigma_eff_jitter(dm)
    assert abs(est - jitter) / jitter < 0.2


def test_trend_guard_floored_scale():
    rng = np.random.default_rng(1)
    t = np.sort(rng.uniform(0, 400, size=40))
    dm = 350 + 0.01 * t + rng.normal(0, 0.02, size=40)
    err = np.full(40, 0.02)
    out = fit_dm_models(t, dm, err, rng=rng, n_perm=120)
    assert out["err_scale"] == "jitter-floored"
    assert out["floor_applied"] is True
    assert out["trend_significant"] is True
    dm2 = 350 + rng.normal(0, 0.05, size=40)
    out2 = fit_dm_models(t, dm2, np.full(40, 0.01), rng=rng, n_perm=120)
    assert out2["trend_significant"] is False


def test_wls_poly_bic_ordering():
    rng = np.random.default_rng(2)
    t = np.linspace(0, 100, 30)
    y = 5 + 0.1 * t
    g, chi2, _ = wls_poly(t, y, np.full(30, 0.1), 1)
    assert abs(g[0] - 0.1) < 1e-6 and chi2 < 1e-6


def test_sinusoid_recovery_mini():
    rng = np.random.default_rng(3)
    P, A, sig = 16.35, 0.04, 0.02
    t = np.sort(rng.uniform(0, 400, size=40))
    dm = 350 + A * np.sin(2 * np.pi * t / P + rng.uniform(0, 6.28)) + rng.normal(0, sig, 40)
    out = fit_dm_models(t, dm, np.full(40, sig / 3), p_burst=P, rng=rng, n_perm=120)
    assert out["sinusoid"]["significant"] is True
    assert 0.5 * A < out["sinusoid"]["amp"] < 1.5 * A