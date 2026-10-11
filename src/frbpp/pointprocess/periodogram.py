"""Periodogram over the observing grid plus efficiency-corrected phase density."""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from frbpp.contracts import ObservationLog

from .likelihood import loglike_constant, loglike_step


def obs_grid_from_obslog(log: ObservationLog, dt_days: float = 60.0 / 86400.0):
    ts, ws = [], []
    for a, b, s in zip(log.t_start, log.t_end, log.sens):
        n = max(int(np.ceil((b - a) / dt_days)), 1)
        t = np.linspace(a, b, n + 1)[:-1]
        ts.append(t)
        ws.append(np.full(n, float(s) * dt_days))
    return np.concatenate(ts), np.concatenate(ws)


def _check_epoch(t, grid_t):
    if t.size and abs(float(np.median(t)) - float(np.median(grid_t))) > 3650.0:
        raise ValueError(
            "事件与观测网格时间基准不一致（疑似 MJD/BJD 混用）: "
            f"median(事件)={float(np.median(t)):.1f}, median(网格)={float(np.median(grid_t)):.1f}")


@dataclass
class StepWindowScanner:

    t0: float
    grid_t: np.ndarray
    grid_w: np.ndarray
    p_grid: np.ndarray
    phi0_grid: np.ndarray
    width_grid: np.ndarray
    t_in_table: np.ndarray
    t_tot: float

    @classmethod
    def build(cls, grid_t, grid_w, t0, p_grid, phi0_grid, width_grid):
        p_grid = np.asarray(p_grid, dtype=float)
        phi0_grid = np.asarray(phi0_grid, dtype=float)
        width_grid = np.asarray(width_grid, dtype=float)
        grid_t = np.asarray(grid_t, dtype=float)
        grid_w = np.asarray(grid_w, dtype=float)
        nphi, nw = len(phi0_grid), len(width_grid)
        tin = np.empty((len(p_grid), nphi, nw), dtype=float)
        for i, P in enumerate(p_grid):
            ph = np.mod((grid_t - t0) / P, 1.0)
            diff = np.mod(ph[None, :] - phi0_grid[:, None], 1.0)
            mask = diff[:, None, :] < width_grid[None, :, None]
            tin[i] = (mask * grid_w[None, None, :]).sum(axis=2)
        return cls(float(t0), grid_t, grid_w, p_grid, phi0_grid, width_grid,
                   tin, float(grid_w.sum()))

    def _mask_counts(self, phases, phi0_grid, width_grid):
        d = np.mod(phases[None, None, :] - phi0_grid[:, None, None], 1.0)
        return (d < width_grid[None, :, None]).sum(axis=2)

    def scan(self, event_times):
        t = np.asarray(event_times, dtype=float)
        _check_epoch(t, self.grid_t)
        n_tot = t.size
        lam = np.empty(len(self.p_grid), dtype=float)
        best = np.empty((len(self.p_grid), 2), dtype=int)
        ll0 = float(loglike_constant(n_tot, self.t_tot))
        for i, P in enumerate(self.p_grid):
            ph = np.mod((t - self.t0) / P, 1.0)
            n_in = self._mask_counts(ph, self.phi0_grid, self.width_grid)
            tin = self.t_in_table[i]
            ll = loglike_step(n_in, tin, n_tot - n_in, self.t_tot - tin)
            j = np.unravel_index(np.argmax(ll), ll.shape)
            lam[i] = 2.0 * (float(ll[j]) - ll0)
            best[i] = j
        return lam, best

    def stat_at(self, event_times, P: float):
        t = np.asarray(event_times, dtype=float)
        _check_epoch(t, self.grid_t)
        n_tot = t.size
        ph = np.mod((t - self.t0) / P, 1.0)
        ph_g = np.mod((self.grid_t - self.t0) / P, 1.0)
        diff = np.mod(ph_g[None, :] - self.phi0_grid[:, None], 1.0)
        mask = diff[:, None, :] < self.width_grid[None, :, None]
        tin = (mask * self.grid_w[None, None, :]).sum(axis=2)
        n_in = self._mask_counts(ph, self.phi0_grid, self.width_grid)
        ll = loglike_step(n_in, tin, n_tot - n_in, self.t_tot - tin)
        j = np.unravel_index(np.argmax(ll), ll.shape)
        return 2.0 * (float(ll[j]) - float(loglike_constant(n_tot, self.t_tot))), j


def kuiper_statistic(event_times, grid_t, grid_w, t0, P):
    t = np.asarray(event_times, dtype=float)
    n = t.size
    if n == 0:
        return 0.0
    ph_e = np.sort(np.mod((t - t0) / P, 1.0))
    ph_g = np.mod((grid_t - t0) / P, 1.0)
    order = np.argsort(ph_g)
    ph_sorted = ph_g[order]
    cw = np.cumsum(np.asarray(grid_w, dtype=float)[order])
    cw = cw / cw[-1]
    idx = np.searchsorted(ph_sorted, ph_e, side="right") - 1
    f_exp = np.where(idx >= 0, cw[np.clip(idx, 0, None)], 0.0)
    f_obs = np.arange(1, n + 1) / n
    return float(np.max(f_obs - f_exp) + np.max(f_exp - f_obs))