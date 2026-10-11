"""Merge bursts closer than t_cluster into episodes."""
from __future__ import annotations

import numpy as np


def decluster(times, t_cluster_days: float = 1.0 / 24.0):
    t = np.sort(np.asarray(times, dtype=float))
    if t.size == 0:
        return t, {"n_episodes": 0, "counts": np.array([], dtype=int),
                   "t_cluster_days": float(t_cluster_days)}
    breaks = np.where(np.diff(t) >= t_cluster_days)[0]
    edges = np.concatenate(([0], breaks + 1, [t.size]))
    ep_t = np.array([t[i] for i in edges[:-1]], dtype=float)
    meta = {
        "n_episodes": int(ep_t.size),
        "n_events": int(t.size),
        "counts": np.diff(edges).astype(int),
        "t_cluster_days": float(t_cluster_days),
    }
    return ep_t, meta