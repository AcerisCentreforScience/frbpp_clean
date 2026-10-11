"""Closed-form MLE and log-likelihood for the stepwise activity-window NHPP."""
from __future__ import annotations

import numpy as np


def _term(n, t):
    n = np.asarray(n, dtype=float)
    t = np.asarray(t, dtype=float)
    n, t = np.broadcast_arrays(n, t)
    val = np.zeros(n.shape, dtype=float)
    pos = n > 0
    good = pos & (t > 0)
    val[good] = n[good] * np.log(n[good] / t[good])
    val[pos & (t <= 0)] = -np.inf
    return val


def loglike_step(n_in, t_in, n_out, t_out):
    n_in = np.asarray(n_in, dtype=float)
    n_out = np.asarray(n_out, dtype=float)
    return _term(n_in, t_in) + _term(n_out, t_out) - (n_in + n_out)


def loglike_constant(n, t_tot):
    n = np.asarray(n, dtype=float)
    return _term(n, t_tot) - n