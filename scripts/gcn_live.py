"""Consume the GCN live stream: verification against pre-registered forecasts plus upgrade bookkeeping.

Usage: python scripts/gcn_live.py --minutes 60
"""
from __future__ import annotations

import argparse
import base64
import datetime as dt
import json
import sys
import time
from pathlib import Path

import requests

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from frbpp.forecast.predict import verify_event                      # noqa: E402

CFG = ROOT / "config" / "gcn_kafka.local.json"
PRED_DIR = ROOT / "predictions"
HIST = PRED_DIR / "HISTORY.jsonl"
LEDGER = PRED_DIR / "GCN_LIVE.jsonl"
CAT2 = ROOT / "data" / "interim" / "cat2_bursts.parquet"
C2V1 = ROOT / "data" / "processed" / "upgrade" / "c2_candidates_v1.parquet"

FRB_TOPICS = ("chime.frb", "dsa110.frb")
CONE_WATCH_ARCSEC = 600.0


def _bjd_from_iso(iso: str) -> float:
    from astropy.time import Time
    return float(Time(str(iso), scale="utc").tdb.jd)


TOKEN_URL = "https://auth.gcn.nasa.gov/oauth2/token"


def _fetch_token(cfg: dict) -> tuple[str, float]:
    r = requests.post(TOKEN_URL,
                      data={"grant_type": "client_credentials",
                            "client_id": cfg["client_id"],
                            "client_secret": cfg["client_secret"]},
                      timeout=30)
    r.raise_for_status()
    j = r.json()
    return str(j["access_token"]), time.time() + float(j.get("expires_in", 3600))


def _jwt_sub(token: str) -> str:
    try:
        p = token.split(".")[1]
        p += "=" * (-len(p) % 4)
        return str(json.loads(base64.urlsafe_b64decode(p)).get("sub", ""))
    except Exception:  # noqa: BLE001
        return ""


def make_oauth_cb(cfg: dict):
    cache: dict = {"token": None, "exp": 0.0}

    def oauth_cb(_oauth_config: str):
        now = time.time()
        if cache["token"] is None or now > cache["exp"] - 60:
            tok, exp = _fetch_token(cfg)
            cache.update(token=tok, exp=exp)
        return (cache["token"], cache["exp"], _jwt_sub(cache["token"]), {})

    return oauth_cb


def parse_frb_notice(v: dict) -> dict:
    g = v.get
    return {
        "trigger_time": g("trigger_time"),
        "alert_id": g("id"),
        "alert_type": g("alert_type") or "initial",
        "snr": g("snr"),
        "ra": g("ra"), "dec": g("dec"), "ra_dec_error": g("ra_dec_error"),
        "dm": g("dm"), "dm_error": g("dm_error"),
        "known_source": g("known_source"),
        "association_probability": g("association_probability"),
        "tns_name": g("tns_name") or g("event_name"),
    }


def _ang_arcsec(ra1, dec1, ra2, dec2):
    r1, d1, r2, d2 = map(np.radians, [ra1, dec1, ra2, dec2])
    x = np.clip(np.sin(d1) * np.sin(d2) + np.cos(d1) * np.cos(d2) * np.cos(r1 - r2), -1, 1)
    return float(np.degrees(np.arccos(x)) * 3600.0)


def _dm_close(dm_a, dm_b, tol) -> bool:
    try:
        a, b = float(dm_a), float(dm_b)
    except (TypeError, ValueError):
        return True
    if not (np.isfinite(a) and np.isfinite(b)):
        return True
    return abs(a - b) < tol


def nearest(coords: pd.DataFrame, ra, dec, radius_arcsec):
    try:
        ra_f, dec_f = float(ra), float(dec)
    except (TypeError, ValueError):
        return None, None, None
    if not (np.isfinite(ra_f) and np.isfinite(dec_f)) or coords.empty:
        return None, None, None
    d = np.array([_ang_arcsec(ra_f, dec_f, r.ra, r.dec) for r in coords.itertuples()])
    i = int(d.argmin())
    if d[i] <= radius_arcsec:
        return str(coords.iloc[i]["name"]), float(d[i]), coords.iloc[i]
    return None, None, None


def load_context():
    cat = pd.read_parquet(CAT2, columns=["tns_name", "ra", "dec", "dm_fitb"])
    cat = cat.dropna(subset=["ra", "dec"])
    coords = cat.groupby("tns_name").first().reset_index().rename(columns={"tns_name": "name"})
    v1 = pd.read_parquet(C2V1)
    rows = []
    for r in v1.itertuples():
        rows.append({"name": r.name1, "ra": r.ra1, "dec": r.dec1, "dm": r.dm1})
        rows.append({"name": r.name2, "ra": r.ra2, "dec": r.dec2, "dm": r.dm2})
    watch = (pd.DataFrame(rows).dropna(subset=["ra", "dec"])
             .drop_duplicates("name").reset_index(drop=True))
    pred_sources = set()
    if (PRED_DIR / "REGISTRY.jsonl").exists():
        for ln in (PRED_DIR / "REGISTRY.jsonl").read_text(encoding="utf-8").splitlines():
            if ln.strip():
                pred_sources.add(json.loads(ln)["source"])
    return coords, watch, pred_sources


def handle_frb_notice(v: dict, ctx) -> dict:
    p = parse_frb_notice(v)
    if p["alert_type"] == "retraction":
        return {"action": "RETRACTION", "parsed": p, "bjd_tdb": None,
                "detail": {"of_alert": p["alert_id"]}}
    action, detail = "unmatched", {}
    t_bjd = None
    if p["trigger_time"]:
        try:
            t_bjd = _bjd_from_iso(p["trigger_time"])
        except Exception:  # noqa: BLE001
            detail["bjd_error"] = True
    name = p["known_source"] or p["tns_name"]
    src = None
    if name and name in ctx["pred_sources"]:
        src = name
        detail["match"] = "name"
    else:
        pcoords = ctx["coords"][ctx["coords"]["name"].isin(ctx["pred_sources"])]
        n, sep, row = nearest(pcoords, p["ra"], p["dec"], 300.0)
        if n is not None and _dm_close(p["dm"], row["dm_fitb"], 30.0):
            src = str(n)
            detail.update({"match": "cone", "cone_matched": sep})
    if src and t_bjd is not None:
        r = verify_event(src, t_bjd, PRED_DIR, HIST)
        action = "PREDICTION_" + r["status"]
        detail.update({"pred_source": src, "duplicate": r["duplicate"]})
    wname, wsep, wrow = nearest(ctx["watch"], p["ra"], p["dec"], CONE_WATCH_ARCSEC)
    if wname and _dm_close(p["dm"], wrow["dm"], 40.0):
        action = "UPGRADE_EVIDENCE"
        detail.update({"watch_name": wname, "new_name": name,
                       "alert_id": p["alert_id"],
                       "sep_arcsec": round(wsep, 1), "dm": p["dm"]})
    return {"action": action, "parsed": p, "bjd_tdb": t_bjd, "detail": detail}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--minutes", type=float, default=5.0)
    args = ap.parse_args()
    cfg = json.loads(CFG.read_text(encoding="utf-8"))
    ctx = {"coords": None, "watch": None, "pred_sources": None}
    ctx["coords"], ctx["watch"], ctx["pred_sources"] = load_context()

    seen = set()
    if LEDGER.exists():
        for ln in LEDGER.read_text(encoding="utf-8").splitlines():
            try:
                e = json.loads(ln)
                seen.add((e.get("topic"), e.get("partition"), e.get("offset")))
            except Exception:  # noqa: BLE001
                continue

    from gcn_kafka import Consumer
    consumer = Consumer(config={"oauth_cb": make_oauth_cb(cfg)},
                        client_id=cfg["client_id"], client_secret=cfg["client_secret"])
    try:
        md = consumer.list_topics(timeout=15)
        print(f"[GCN] 连接探针 OK：brokers={len(md.brokers)} topics={len(md.topics)}")
    except Exception as e:  # noqa: BLE001
        print(f"[GCN] 连接探针失败：{e}")
    consumer.subscribe(list(cfg.get("topics", [])))
    print(f"[GCN] 已订阅 {cfg['topics']}；运行 {args.minutes} 分钟...")
    t_end = time.time() + args.minutes * 60
    stats = {"messages": 0, "heartbeat": 0, "frb": 0, "actions": {}}
    LEDGER.parent.mkdir(parents=True, exist_ok=True)
    with LEDGER.open("a", encoding="utf-8") as led:
        while time.time() < t_end:
            for message in consumer.consume(timeout=1):
                if message.error():
                    print(f"[ERR] {message.error()}")
                    continue
                topic = message.topic()
                key = (topic, message.partition(), message.offset())
                stats["messages"] += 1
                if key in seen:
                    continue
                seen.add(key)
                row = {"utc": dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds"),
                       "topic": topic, "partition": message.partition(),
                       "offset": message.offset()}
                if any(t in topic for t in FRB_TOPICS):
                    try:
                        v = json.loads(message.value())
                    except Exception as e:  # noqa: BLE001
                        v = {}
                        row["parse_error"] = str(e)
                    stats["frb"] += 1
                    r = handle_frb_notice(v, ctx)
                    row.update(r)
                    stats["actions"][r["action"]] = stats["actions"].get(r["action"], 0) + 1
                    print(f"[FRB] {topic} action={r['action']} name={r['parsed'].get('tns_name')} "
                          f"detail={r['detail']}")
                else:
                    stats["heartbeat"] += 1
                    row["action"] = "heartbeat"
                led.write(json.dumps(row, ensure_ascii=False) + "\n")
                led.flush()
    consumer.close()
    print(f"[GCN] 结束：消息 {stats['messages']}（FRB {stats['frb']} / heartbeat "
          f"{stats['heartbeat']}）| 动作 {stats['actions']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())