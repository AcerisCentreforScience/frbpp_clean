"""Activity-window extrapolation, pre-registration records, event verification accounting."""
from __future__ import annotations

import datetime as dt
import hashlib
import json
from pathlib import Path

import numpy as np


def window_occurs(t, t0, P, phi0, w):
    ph = np.mod((np.asarray(t, dtype=float) - t0) / P, 1.0)
    return np.mod(ph - phi0, 1.0) < w


def next_windows(t_now, t0, P, phi0, w, n_out=5):
    k = int(np.floor((t_now - t0) / P)) - 1
    out = []
    while len(out) < n_out:
        st = t0 + (k + phi0) * P
        en = t0 + (k + phi0 + w) * P
        if en > t_now:
            out.append((st, en))
        k += 1
    return out


def make_prediction(source, P, phi0, w, t0, t_now, n_out=5, meta=None,
                    file_stem=None):
    wins = next_windows(t_now, t0, P, phi0, w, n_out)
    tol = max(0.25, 0.5 * w * P)
    pred = {
        "source": str(source), "P_days": float(P), "phi0": float(phi0),
        "width_frac": float(w), "t0_bjd": float(t0),
        "windows": [[float(a), float(b)] for a, b in wins],
        "tolerance_days": float(tol),
        "scoring_rule": "HIT = t ∈ [start−tol, end+tol]，tol = max(0.25 d, 0.5·窗长)",
        "created_utc": dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds"),
    }
    if meta:
        pred.update(meta)
    if file_stem:
        pred["file_stem"] = file_stem
    return pred


def _idem_key(pred):
    s = (f"{pred['source']}|{pred['P_days']:.6f}|{pred['phi0']:.6f}|"
         f"{pred['width_frac']:.4f}|{pred['windows'][0][0]:.5f}")
    return hashlib.sha256(s.encode()).hexdigest()[:16]


def register(pred, registry_dir):
    d = Path(registry_dir)
    d.mkdir(parents=True, exist_ok=True)
    reg = d / "REGISTRY.jsonl"
    key = _idem_key(pred)
    if reg.exists():
        for line in reg.read_text(encoding="utf-8").splitlines():
            if not line.strip():
                continue
            e = json.loads(line)
            if e.get("key") == key:
                return "exists", key, e.get("sha256", "")
    stem = pred.get("file_stem") or f"PRED_{pred['source']}"
    fp = d / f"{stem}.json"
    body = json.dumps(pred, sort_keys=True, ensure_ascii=False)
    sha = hashlib.sha256(body.encode()).hexdigest()
    fp.write_text(json.dumps(pred, indent=2, ensure_ascii=False), encoding="utf-8")
    entry = {"key": key, "file": fp.name, "sha256": sha,
             "created_utc": pred["created_utc"], "source": pred["source"],
             "P_days": pred["P_days"]}
    with reg.open("a", encoding="utf-8") as f:
        f.write(json.dumps(entry, ensure_ascii=False) + "\n")
    return "registered", key, sha


def verify_event(source, t_bjd, registry_dir, history_path):
    d = Path(registry_dir)
    preds = []
    if d.exists():
        for fp in sorted(d.glob("PRED_*.json")):
            try:
                p = json.loads(fp.read_text(encoding="utf-8"))
            except Exception:  # noqa: BLE001
                continue
            if p.get("source") == source:
                preds.append(p)
    status, matched = "NO_MATCH", None
    for p in preds:
        tol = p["tolerance_days"]
        for a, b in p["windows"]:
            if a - tol <= t_bjd <= b + tol:
                status, matched = "HIT", {"file": p.get("file_stem", p["source"]),
                                          "window": [a, b], "P_days": p["P_days"]}
                break
        if status == "HIT":
            break
    if preds and status != "HIT":
        status = "MISS"
    hist = Path(history_path)
    hist.parent.mkdir(parents=True, exist_ok=True)
    ekey = hashlib.sha256(f"{source}|{t_bjd:.4f}".encode()).hexdigest()[:16]
    dup = False
    if hist.exists():
        dup = any(ekey in ln for ln in hist.read_text(encoding="utf-8").splitlines())
    row = {"event_key": ekey, "source": source, "t_bjd": float(t_bjd),
           "status": status, "matched": matched,
           "logged_utc": dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds")}
    if not dup:
        with hist.open("a", encoding="utf-8") as f:
            f.write(json.dumps(row, ensure_ascii=False) + "\n")
    return {**row, "duplicate": dup}