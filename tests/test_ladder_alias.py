"""Tests for the aliasing guard and the model ladder."""
import json
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from frbpp.models.ladder import (fit_m2, ladder_fit)          # noqa: E402
from frbpp.pointprocess.alias import alias_test               # noqa: E402
from frbpp.pointprocess.null import sample_periodic_signal    # noqa: E402
from frbpp.pointprocess.periodogram import obs_grid_from_obslog  # noqa: E402
from frbpp.contracts import ObservationLog                    # noqa: E402


def _small_grid(n_days=400):
    starts = np.arange(0.0, n_days, 0.99727)
    t = np.concatenate([s + np.arange(6) * (720.0 / 86400.0) for s in starts])
    w = np.full(t.size, 720.0 / 86400.0)
    return starts, t, w


def test_alias_capture_for_sidereal_beat():
    starts, _, _ = _small_grid()
    ones = np.ones_like(starts)
    r = alias_test(starts, ones, p_hat=3 * 0.99727)
    assert r["flagged"] and r["rank"] >= 0.995
    assert 0.2 < r["rho_hat"] < 0.5
    rng = np.random.default_rng(0)
    r2 = alias_test(rng.uniform(0, 400, size=starts.size), ones, p_hat=3 * 0.99727)
    assert not r2["flagged"]


def test_epoch_guard_raises():
    _, t, w = _small_grid()
    ev = np.array([2459000.0, 2459010.0])
    with pytest.raises(ValueError):
        ladder_fit(ev, float(ev[0]), 16.35, t, w)


def test_ladder_m1_form_recovery_and_anchors():
    starts, t, w = _small_grid()
    log = ObservationLog(t_start=t, t_end=t + 720.0 / 86400.0, sens=np.ones(t.size),
                         tgt_ra=10.0, tgt_dec=60.0,
                         beam_ra=np.full(t.size, np.nan), beam_dec=np.full(t.size, np.nan))
    grid_t, grid_w = t, w
    rng = np.random.default_rng(1)
    ev = sample_periodic_signal(grid_t, grid_w, 0.0, 16.35, 0.2, 0.15, 50.0, 40, rng)
    res = ladder_fit(ev, float(ev[0]), 16.35, grid_t, grid_w)
    assert res["final"] in ("M1", "S1")
    assert res["final"] == "M1"
    assert res["table"]["M1"]["det"]["width"] <= 0.5 + 1e-9
    ll2, _, d2 = fit_m2(ev, float(ev[0]), 16.35, grid_t, grid_w)
    assert d2["A"] < 1.0


def test_ladder_deterministic():
    starts, t, w = _small_grid()
    rng = np.random.default_rng(2)
    ev = sample_periodic_signal(t, w, 0.0, 16.35, 0.1, 0.2, 50.0, 30, rng)
    r1 = ladder_fit(ev, float(ev[0]), 16.35, t, w)
    r2 = ladder_fit(ev, float(ev[0]), 16.35, t, w)
    assert json.dumps(r1, sort_keys=True, default=str) == \
        json.dumps(r2, sort_keys=True, default=str)