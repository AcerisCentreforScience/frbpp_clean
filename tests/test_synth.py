"""Unit tests for the synthetic NHPP generator."""
import numpy as np

from frbpp.synth import (
    WindowExposure,
    constant_intensity,
    sample_nhpp,
    sample_uniform_given_exposure,
    square_window_intensity,
)

SEED = 20261007


def _trapz(y, t):
    y = np.asarray(y, dtype=float)
    return float(np.sum(0.5 * (y[1:] + y[:-1]) * np.diff(t)))


def test_uniform_given_exposure_density_and_confinement():
    expo = WindowExposure(t_start=[0.0, 10.0, 30.0],
                          t_end=[10.0, 20.0, 50.0],
                          sens=[1.0, 0.5, 0.25])
    rng = np.random.default_rng(SEED)
    x = sample_uniform_given_exposure(expo, n=50_000, rng=rng, T0=0.0, T1=50.0)
    assert np.all(expo(x) > 0.0)
    w = (expo.t_end - expo.t_start) * expo.sens
    p = w / w.sum()
    for wi, a, b in zip(p, expo.t_start, expo.t_end):
        frac = float(np.mean((x >= a) & (x < b)))
        tol = 4.0 * float(np.sqrt(wi * (1.0 - wi) / x.size)) + 2e-3
        assert abs(frac - wi) < tol


def test_constant_rate_moments():
    expo = WindowExposure(t_start=[0.0, 200.0, 400.0],
                          t_end=[100.0, 300.0, 600.0],
                          sens=[1.0, 0.8, 0.6])
    lam = constant_intensity(0.5)
    T0, T1 = 0.0, 600.0
    mu_int = 0.5 * expo.integral
    rng = np.random.default_rng(SEED)
    K = 400
    counts = np.array([sample_nhpp(lam, expo, T0, T1, rng=rng).size
                       for _ in range(K)])
    assert abs(counts.mean() - mu_int) < 4.0 * np.sqrt(mu_int / K)
    ratio = counts.var(ddof=1) / counts.mean()
    assert 0.7 < ratio < 1.35
    t = sample_nhpp(lam, expo, T0, T1, rng=np.random.default_rng(SEED + 1))
    assert np.all(expo(t) > 0.0)


def test_periodic_phase_and_count_recovery():
    P, w = 10.0, 0.25
    lam = square_window_intensity(c=0.2, a=1.8, period=P, t0=0.0, width_frac=w)
    expo = WindowExposure(t_start=[0.0, 2000.0], t_end=[1000.0, 4000.0],
                          sens=[1.0, 0.7])
    T0, T1 = 0.0, 4000.0
    grid = np.linspace(T0, T1, 400_001)
    mu = lam(grid) * expo(grid)
    in_window = (np.mod(grid / P, 1.0) < w)
    exp_n = _trapz(mu, grid)
    exp_in = _trapz(mu * in_window, grid)
    rng = np.random.default_rng(SEED)
    t = sample_nhpp(lam, expo, T0, T1, rng=rng)
    assert np.all(expo(t) > 0.0)
    assert not np.any((t > 1000.0) & (t < 2000.0))
    assert abs(t.size - exp_n) < 4.0 * np.sqrt(exp_n) + 2.0
    f_pred = exp_in / exp_n
    f_obs = float(np.mean(np.mod(t / P, 1.0) < w))
    tol = 4.0 * float(np.sqrt(f_pred * (1.0 - f_pred) / t.size)) + 5e-3
    assert abs(f_obs - f_pred) < tol