"""Tests for the false-positive battery and completeness surface."""
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from frbpp.pointprocess.periodogram import StepWindowScanner   # noqa: E402
from frbpp.validation.fpr import fpr_battery_n1, fpr_battery_n2  # noqa: E402
from frbpp.validation.inject_recover import (                  # noqa: E402
    completeness_surface, sample_gaussian_window_signal)


def _small_scanner():
    rng = np.random.default_rng(0)
    centers = np.cumsum(rng.uniform(0.9, 1.1, size=120))
    windows = [(c, c + 1.0 / 24.0) for c in centers]
    grid_t = np.concatenate([np.linspace(a, b, 60, endpoint=False)
                             for a, b in windows])
    grid_w = np.full(grid_t.size, (1.0 / 24.0) / 60.0)
    sc = StepWindowScanner.build(grid_t, grid_w, 0.0, np.linspace(2.0, 6.0, 300),
                                 np.arange(12) / 12.0, np.array([0.1, 0.25]))
    return sc, grid_t, grid_w


def test_fpr_batteries_small_run():
    sc, gt, gw = _small_scanner()
    rng = np.random.default_rng(2)
    b1 = fpr_battery_n1(sc, gt, gw, 8, n_real=12, n_null=40, alpha=0.05, rng=rng)
    assert {"n_false", "fpr", "p_values"}.issubset(b1.keys())
    assert b1["p_values"].size == 12
    b2 = fpr_battery_n2(sc, gt, gw, 8, n_real=6, n_surr=30, alpha=0.05, rng=rng)
    assert b2["p_values"].size == 6


def test_gaussian_sampler_concentration():
    sc, gt, gw = _small_scanner()
    rng = np.random.default_rng(3)
    ev = sample_gaussian_window_signal(gt, gw, 0.0, 4.0, 0.2, 0.05, 50.0, 200, rng)
    dphi = np.mod(np.mod(ev / 4.0, 1.0) - 0.2 + 0.5, 1.0) - 0.5
    assert float(np.mean(np.abs(dphi) < 0.05)) > 0.5


def test_completeness_surface_shapes():
    sc, gt, gw = _small_scanner()
    rng = np.random.default_rng(4)
    rows = completeness_surface(sc, gt, gw, 0.0, [{"P": 4.0, "width": 0.25, "n_ep": 40}],
                                n_rep=3, contrast=50.0, tol_frac=0.05, rng=rng)
    assert len(rows) == 1
    assert rows[0]["hits"] in (0, 1, 2, 3)
    assert len(rows[0]["ci95"]) == 2