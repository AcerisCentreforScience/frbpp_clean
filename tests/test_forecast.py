"""Tests for activity-window extrapolation, pre-registration, verification accounting."""
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from frbpp.forecast.predict import (make_prediction, next_windows,   # noqa: E402
                                    register, verify_event, window_occurs)


def test_window_math():
    t0, P, phi0, w = 100.0, 10.0, 0.2, 0.3
    t = np.array([100.0, 102.0, 104.9, 105.1, 113.0])
    occ = window_occurs(t, t0, P, phi0, w)
    assert occ.tolist() == [False, True, True, False, True]
    wins = next_windows(103.0, t0, P, phi0, w, n_out=2)
    assert wins[0] == (102.0, 105.0) and wins[1] == (112.0, 115.0)


def test_register_idempotent(tmp_path):
    pred = make_prediction("SRC", 16.35, 0.1, 0.2, 1000.0, 2000.0, n_out=3,
                           file_stem="PRED_SRC")
    st1, k1, s1 = register(pred, tmp_path)
    st2, k2, s2 = register(pred, tmp_path)
    assert st1 == "registered" and st2 == "exists" and k1 == k2
    reg = (tmp_path / "REGISTRY.jsonl").read_text(encoding="utf-8").splitlines()
    assert len(reg) == 1 and json.loads(reg[0])["sha256"] == s1


def test_register_epoch_distinct(tmp_path):
    p1 = make_prediction("SRC", 16.35, 0.1, 0.2, 1000.0, 2000.0, n_out=3,
                         file_stem="PRED_SRC_g1")
    p2 = make_prediction("SRC", 16.35, 0.1, 0.2, 1000.0, 2060.0, n_out=3,
                         file_stem="PRED_SRC_g2")
    st1, k1, _ = register(p1, tmp_path)
    st2, k2, _ = register(p2, tmp_path)
    st3, k3, _ = register(p2, tmp_path)
    assert (st1, st2, st3) == ("registered", "registered", "exists")
    assert k1 != k2 and k3 == k2
    assert (tmp_path / "PRED_SRC_g1.json").exists()
    assert (tmp_path / "PRED_SRC_g2.json").exists()
    reg = (tmp_path / "REGISTRY.jsonl").read_text(encoding="utf-8").splitlines()
    assert len(reg) == 2


def test_verify_event_statuses(tmp_path):
    pred = make_prediction("SRC", 10.0, 0.2, 0.3, 0.0, 0.0, n_out=2,
                           file_stem="PRED_SRC")
    register(pred, tmp_path)
    hist = tmp_path / "HISTORY.jsonl"
    a, b = pred["windows"][0]
    r1 = verify_event("SRC", 0.5 * (a + b), tmp_path, hist)
    r2 = verify_event("SRC", a - 50.0, tmp_path, hist)
    r3 = verify_event("OTHER", 0.5 * (a + b), tmp_path, hist)
    r4 = verify_event("SRC", 0.5 * (a + b), tmp_path, hist)
    assert r1["status"] == "HIT" and r2["status"] == "MISS"
    assert r3["status"] == "NO_MATCH" and r4["duplicate"] is True
    assert len(hist.read_text(encoding="utf-8").splitlines()) == 3