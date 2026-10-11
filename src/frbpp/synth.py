"""Synthetic NHPP generator used by the acceptance gates."""
from __future__ import annotations

import numpy as np


def sample_uniform_given_exposure(exposure, n, rng=None, T0=0.0, T1=1.0):
    if rng is None:
        rng = np.random.default_rng()
    out, got = [], 0
    while got < n:
        cand = rng.uniform(T0, T1, size=max(2 * (n - got), 64))
        p = np.asarray(exposure(cand), dtype=float)
        if p.size and (p.min() < 0.0 or p.max() > 1.0):
            raise ValueError("exposure 必须取值于 [0, 1]")
        acc = cand[rng.random(cand.size) < p]
        got += acc.size
        out.append(acc)
    pool = np.concatenate(out)
    if pool.size > n:
        pool = rng.choice(pool, size=n, replace=False)
    return np.sort(pool)


def sample_nhpp(intensity, exposure, T0, T1, rate_bound=None, rng=None):
    if rng is None:
        rng = np.random.default_rng()
    if rate_bound is None:
        grid = np.linspace(T0, T1, 8192)
        mu_grid = (np.asarray(intensity(grid), dtype=float)
                   * np.asarray(exposure(grid), dtype=float))
        rate_bound = float(mu_grid.max()) * 1.05 + 1e-12
    n_cand = int(rng.poisson(rate_bound * (T1 - T0)))
    cand = rng.uniform(T0, T1, size=n_cand)
    p = (np.asarray(intensity(cand), dtype=float)
         * np.asarray(exposure(cand), dtype=float) / rate_bound)
    keep = rng.random(cand.size) < np.clip(p, 0.0, 1.0)
    return np.sort(cand[keep])


class WindowExposure:

    def __init__(self, t_start, t_end, sens):
        t_start = np.asarray(t_start, dtype=float)
        t_end = np.asarray(t_end, dtype=float)
        sens = np.asarray(sens, dtype=float)
        if not (len(t_start) == len(t_end) == len(sens)):
            raise ValueError("t_start / t_end / sens 长度必须一致")
        order = np.argsort(t_start)
        self.t_start = t_start[order]
        self.t_end = t_end[order]
        self.sens = sens[order]
        if np.any(self.t_end <= self.t_start):
            raise ValueError("存在非正长度观测窗")
        if np.any(self.sens <= 0.0) or np.any(self.sens > 1.0):
            raise ValueError("sens 必须取值于 (0, 1]")
        if len(self.t_start) > 1 and np.any(self.t_start[1:] < self.t_end[:-1]):
            raise ValueError("观测窗必须互不重叠（窗口守卫精神）")

    def __call__(self, t):
        t = np.asarray(t, dtype=float)
        out = np.zeros(t.shape, dtype=float)
        for a, b, s in zip(self.t_start, self.t_end, self.sens):
            out[(t >= a) & (t < b)] = s
        return out

    @property
    def integral(self) -> float:
        return float(np.sum((self.t_end - self.t_start) * self.sens))


def constant_intensity(rate):
    rate = float(rate)
    if rate <= 0.0:
        raise ValueError("rate 必须为正")

    def intensity(t):
        t = np.asarray(t, dtype=float)
        return np.full(t.shape, rate)

    return intensity


def square_window_intensity(c, a, period, t0=0.0, width_frac=0.25):
    c = float(c)
    a = float(a)
    period = float(period)
    if c <= 0.0 or a < 0.0:
        raise ValueError("要求 c > 0, a >= 0")
    if not 0.0 < width_frac <= 0.5:
        raise ValueError("活动窗宽度 ≤ P/2（物理锚点）")
    if period <= 0.0:
        raise ValueError("period 必须为正")

    def intensity(t):
        ph = np.mod((np.asarray(t, dtype=float) - t0) / period, 1.0)
        w = (ph < width_frac).astype(float)
        return c + a * w

    return intensity