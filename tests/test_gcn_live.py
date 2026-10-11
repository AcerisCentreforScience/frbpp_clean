"""Tests for GCN alert parsing, forecast verification, upgrade bookkeeping."""
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

import gcn_live                                                    # noqa: E402
from frbpp.forecast.predict import make_prediction, register       # noqa: E402


def _ctx(tmp_path, monkeypatch):
    pred = make_prediction("SRC_A", 16.35, 0.1, 0.2, 2459030.0, 2459030.0,
                           n_out=3, file_stem="PRED_SRC_A")
    register(pred, tmp_path)
    monkeypatch.setattr(gcn_live, "PRED_DIR", tmp_path)
    monkeypatch.setattr(gcn_live, "HIST", tmp_path / "HISTORY.jsonl")
    coords = pd.DataFrame([{"name": "SRC_A", "ra": 10.0, "dec": 20.0, "dm_fitb": 500.0}])
    watch = pd.DataFrame([{"name": "WATCH_1", "ra": 100.0, "dec": -30.0, "dm": 300.0}])
    return {"coords": coords, "watch": watch, "pred_sources": {"SRC_A"}}, pred


def _iso(bjd):
    from astropy.time import Time
    return Time(bjd, format="jd", scale="tdb").utc.iso


def test_subsequent_hit_and_duplicate(tmp_path, monkeypatch):
    ctx, pred = _ctx(tmp_path, monkeypatch)
    a, b = pred["windows"][0]
    v = {"trigger_time": _iso(0.5 * (a + b)), "id": "111", "alert_type": "subsequent",
         "known_source": "SRC_A", "association_probability": 0.99,
         "ra": 10.0, "dec": 20.0, "ra_dec_error": [0.1, 0.1, 0], "dm": 500.0,
         "dm_error": 1.0, "snr": 12.0}
    r = gcn_live.handle_frb_notice(v, ctx)
    assert r["action"] == "PREDICTION_HIT"
    r2 = gcn_live.handle_frb_notice(v, ctx)
    assert r2["detail"]["duplicate"] is True


def test_retraction_ignored(tmp_path, monkeypatch):
    ctx, pred = _ctx(tmp_path, monkeypatch)
    a, b = pred["windows"][0]
    v = {"trigger_time": _iso(0.5 * (a + b)), "id": "112",
         "alert_type": "retraction", "known_source": "SRC_A",
         "ra": 10.0, "dec": 20.0, "dm": 500.0}
    r = gcn_live.handle_frb_notice(v, ctx)
    assert r["action"] == "RETRACTION"
    assert not (tmp_path / "HISTORY.jsonl").exists()


def test_cone_fallback_numeric_id(tmp_path, monkeypatch):
    ctx, pred = _ctx(tmp_path, monkeypatch)
    a, b = pred["windows"][0]
    t = _iso(0.5 * (a + b))
    v = {"trigger_time": t, "id": "113", "alert_type": "subsequent",
         "known_source": "231298202", "ra": 10.05, "dec": 20.0,
         "ra_dec_error": [0.2, 0.2, 0], "dm": 501.0, "snr": 15.0}
    assert gcn_live.handle_frb_notice(v, ctx)["action"] == "PREDICTION_HIT"
    v_bad = dict(v, dm=560.0)
    assert gcn_live.handle_frb_notice(v_bad, ctx)["action"] != "PREDICTION_HIT"


def test_miss_upgrade_and_unmatched(tmp_path, monkeypatch):
    ctx, pred = _ctx(tmp_path, monkeypatch)
    far = _iso(pred["windows"][0][0] - 50.0)
    v = {"trigger_time": far, "id": "114", "alert_type": "subsequent",
         "known_source": "SRC_A", "ra": 10.0, "dec": 20.0, "dm": 500.0}
    assert gcn_live.handle_frb_notice(v, ctx)["action"] == "PREDICTION_MISS"
    v2 = {"trigger_time": far, "id": "115", "alert_type": "initial",
          "known_source": None, "ra": 100.02, "dec": -30.0, "dm": 305.0}
    r = gcn_live.handle_frb_notice(v2, ctx)
    assert r["action"] == "UPGRADE_EVIDENCE"
    assert r["detail"]["watch_name"] == "WATCH_1"
    v3 = dict(v2, dm=900.0)
    assert gcn_live.handle_frb_notice(v3, ctx)["action"] == "unmatched"
    v4 = {"trigger_time": far, "id": "116", "alert_type": "initial",
          "ra": 250.0, "dec": 60.0, "dm": 500.0}
    assert gcn_live.handle_frb_notice(v4, ctx)["action"] == "unmatched"