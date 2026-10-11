"""Register a rolled forecast generation with frozen parameters.

Usage: python scripts/register_forecast_roll.py
"""
from __future__ import annotations

import datetime as dt
import json
import sys
from pathlib import Path

from astropy.time import Time

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from frbpp.forecast.predict import make_prediction, register  # noqa: E402

PRED = ROOT / "predictions"


def latest_per_group() -> dict:
    reg = PRED / "REGISTRY.jsonl"
    rows = [json.loads(ln) for ln in
            reg.read_text(encoding="utf-8").splitlines() if ln.strip()]
    groups: dict = {}
    for e in rows:
        p = json.loads((PRED / e["file"]).read_text(encoding="utf-8"))
        gkey = (p["source"], p.get("branch", "scan"))
        if gkey not in groups or e["created_utc"] > groups[gkey][0]["created_utc"]:
            groups[gkey] = (e, p)
    return groups


def main() -> int:
    now = dt.datetime.now(dt.timezone.utc)
    t_now = float(Time(now).tdb.jd)
    tag = now.strftime("%Y%m%d")
    for (src, branch), (_e, p) in sorted(latest_per_group().items()):
        stem = f"PRED_{src}_{branch}_g{tag}"
        pred = make_prediction(
            src, p["P_days"], p["phi0"], p["width_frac"], p["t0_bjd"], t_now,
            n_out=5,
            meta={"branch": branch, "gate": "G8-roll",
                  "n_episodes_train": p.get("n_episodes_train"),
                  "sigma_P": p.get("sigma_P"),
                  "data": p.get("data", "Cat2 (via census)")},
            file_stem=stem)
        st, key, sha = register(pred, PRED)
        st2, _k2, _s2 = register(pred, PRED)
        print(f"[{src} {branch}] {st}（幂等复检 {st2}）key={key} sha={sha[:12]}")
        for a, b in pred["windows"]:
            d0 = Time(a, format="jd", scale="tdb").utc.strftime("%Y-%m-%d")
            d1 = Time(b, format="jd", scale="tdb").utc.strftime("%Y-%m-%d")
            print(f"   窗口: {d0} → {d1}（BJD_TDB {a:.2f}–{b:.2f}）")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())