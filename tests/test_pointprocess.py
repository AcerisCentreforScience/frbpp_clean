"""Tests for de-clustering, NHPP likelihood, BH FDR, null sampling, periodogram."""
import sys
from pathlib import Path

import numpy as np
import pytest  # noqa: F401

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from frbpp.pointprocess.decluster import decluster              # noqa: E402
from frbpp.pointprocess.fdr import bh_fdr                       # noqa: E402
from frbpp.pointprocess.likelihood import (                     # noqa: E402
    loglike_constant, loglike_step)
from frbpp.pointprocess.null import (                           # noqa: E402
    empirical_p, sample_periodic_signal, sample_uniform_given_exposure)
from frbpp.pointprocess.periodogram import StepWindowScanner    # noqa: E402


def _small_grid(n_windows: int = 10):
    windows = [(i * 1.0, i * 1.0 + 1.0 / 24.0) for i in range(n_windows)]
    grid_t = np.concatenate([np.linspace(a, b, 60, endpoint=False)
                             for a, b in windows])
    grid_w = np.full(grid_t.size, (1.0 / 24.0) / 60.0)
    return windows, grid_t, grid_w


def test_decluster_merges_clusters():
    t = np.array([0.0, 0.01, 0.02, 1.0, 1.03, 5.0])
    ep, meta = decluster(t, t_cluster_days=0.05)
    assert ep.tolist() == [0.0, 1.0, 5.0]
    assert meta["counts"].tolist() == [3, 2, 1]
    assert meta["n_episodes"] == 3


def test_nhpp_likelihood_analytic():
    ll = float(loglike_step(2.0, 1.0, 3.0, 3.0))
    assert np.isclose(ll, 2.0 * np.log(2.0) + 3.0 * np.log(1.0) - 5.0)
    assert np.isclose(float(loglike_constant(5.0, 4.0)), 5.0 * np.log(5.0 / 4.0) - 5.0)
    assert np.isclose(float(loglike_step(0.0, 0.0, 3.0, 3.0)), -3.0)
    assert np.isneginf(float(loglike_step(1.0, 0.0, 3.0, 3.0)))


def test_bh_fdr_properties():
    p = np.array([0.001, 0.01, 0.02, 0.04, 0.5])
    q = bh_fdr(p)
    assert np.all(q >= p - 1e-12)
    assert np.all(np.diff(q[np.argsort(p)]) >= -1e-12)
    assert np.isclose(q[0], 0.001 * 5 / 1)
    q2 = bh_fdr(np.full(4, 0.9))
    assert np.all(q2 <= 1.0)


def test_empirical_p_conservative():
    null = np.arange(1, 101, dtype=float)
    assert np.isclose(empirical_p(101.0, null), 1.0 / 101.0)
    assert np.isclose(empirical_p(50.0, null), 52.0 / 101.0)
    assert empirical_p(0.0, null) == 1.0


def test_scanner_finds_injected_window():
    rng = np.random.default_rng(0)
    centers = np.cumsum(rng.uniform(0.9, 1.1, size=160))
    windows = [(c, c + 1.0 / 24.0) for c in centers]
    grid_t = np.concatenate([np.linspace(a, b, 60, endpoint=False)
                             for a, b in windows])
    grid_w = np.full(grid_t.size, (1.0 / 24.0) / 60.0)
    sc = StepWindowScanner.build(grid_t, grid_w, 0.0, np.linspace(2.0, 6.0, 400),
                                 np.arange(12) / 12.0, np.array([0.1, 0.25]))
    ph = np.mod(grid_t / 4.0, 1.0)
    cand = grid_t[(ph >= 0.05) & (ph < 0.30)]
    ev_sig = np.sort(rng.choice(cand, size=12, replace=False))
    ev_unif = sample_uniform_given_exposure(grid_t, grid_w, 12, rng)
    assert sc.stat_at(ev_sig, 4.0)[0] > 2.0 * sc.stat_at(ev_unif, 4.0)[0]
    assert abs(sc.p_grid[int(np.argmax(sc.scan(ev_sig)[0]))] - 4.0) <= 0.05


def test_null_sampler_window_guard_and_contrast():
    windows, grid_t, grid_w = _small_grid()
    rng = np.random.default_rng(1)
    s = sample_uniform_given_exposure(grid_t, grid_w, 50, rng)
    assert np.all([any(a <= x <= b for a, b in windows) for x in s])
    ev = sample_periodic_signal(grid_t, grid_w, 0.0, 4.0, 0.2, 0.25, 50.0, 100, rng)
    frac = float(np.mean(np.mod(np.mod(ev / 4.0, 1.0) - 0.2, 1.0) < 0.25))
    assert frac > 0.8