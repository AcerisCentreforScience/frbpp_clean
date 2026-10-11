"""Period and window hypothesis tests for the 113.69 d pair.

Usage: python scripts/period_113d.py
Output: data/processed/upgrade/pair113_period.json
"""
from __future__ import annotations

import datetime as dt
import json
import sys
from pathlib import Path

import numpy as np
from astropy.time import Time

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from frbpp.catalogs.chime import mjd_to_bjd_tdb                      # noqa: E402
from frbpp.forecast.predict import make_prediction, register         # noqa: E402

OUT = ROOT / "data" / "processed" / "upgrade"
PRED = ROOT / "predictions"
CAT2 = ROOT / "data" / "interim" / "cat2_bursts.parquet"

NAME1, NAME2 = "FRB20200919A", "FRB20210111A"
MJD1, MJD2 = 59111.909114, 59225.596722
DT = MJD2 - MJD1
T_NOW_MJD = 61322.0

BRANCHES = {
    "fundamental": DT,
    "harmonic": DT / 2.0,
    "double": DT * 2.0,
}


def main() -> int:
    import pandas as pd

    r1 = pd.read_parquet(CAT2)
    r1 = r1[r1.tns_name == NAME1].iloc[0]
    ra, dec = float(r1.ra), float(r1.dec)
    t0_bjd = mjd_to_bjd_tdb(MJD1, ra, dec)
    t2_bjd = mjd_to_bjd_tdb(MJD2, ra, dec)
    T_NOW_BJD = mjd_to_bjd_tdb(T_NOW_MJD, ra, dec)

    out = {"built_at": dt.datetime.now().astimezone().isoformat(timespec="seconds"),
           "pair": [NAME1, NAME2], "dt_days": float(DT),
           "t0_bjd_first_burst": float(t0_bjd), "t2_bjd_second_burst": float(t2_bjd),
           "degeneracy_note": "2 事件 → 任何 P=Δt/k 都恒相位自洽；本表不构成显著性证据",
           "branches": {}}

    for tag, P in BRANCHES.items():
        ph1 = float(np.mod((t0_bjd - t0_bjd) / P, 1.0))
        ph2 = float(np.mod((t2_bjd - t0_bjd) / P, 1.0))
        dphi = float(np.mod(ph2 - ph1 + 0.5, 1.0) - 0.5)
        coherent = bool(abs(dphi) < 0.05)
        rec = {"P_days": P, "phase_first": ph1, "phase_second": ph2,
               "phase_diff": dphi, "phase_coherent": coherent,
               "phase_shift_from_dt": DT / P}
        out["branches"][tag] = rec

    W = 0.10
    out["assumed_width_frac"] = W
    out["assumed_width_note"] = "宽度由 2 事件不可约束；取 0.10 作名义值，属假设非测量"

    for tag in ("fundamental", "harmonic"):
        P = BRANCHES[tag]
        pred = make_prediction(
            NAME1, P, 0.0, W, t0_bjd, T_NOW_BJD, n_out=5,
            file_stem=f"PRED_{NAME1}_113d_{tag}",
            meta={
                "gate": "C2-pair113 候选 watch 预注册",
                "eligible_g8": False,
                "eligible_g8_reason": "2 事件 < n_ep≥8；span 不满足 ≥3P；相位检验退化",
                "basis": "C2 时空锥配对间隔 Δt=113.69 d（ΔDM 3.4σ、θ 20.1″）",
                "partner": NAME2,
                "branch": tag,
                "degenerate": True,
                "note": "前瞻试次 = 该源下一次被 GCN 记录的爆发是否落入预测窗；命中=支持成窗，不命中=约束窗口假设",
            })
        status, key, sha = register(pred, PRED)
        wins = [[float(a), float(b)] for a, b in pred["windows"]]
        wins_utc = [[str(Time(a, format="jd", scale="tdb").utc.iso[:10]),
                     str(Time(b, format="jd", scale="tdb").utc.iso[:10])] for a, b in wins]
        out["branches"][tag].update({"register_status": status, "idem_key": key,
                                     "sha256": sha, "file_stem": pred["file_stem"],
                                     "windows_bjd": wins, "windows_utc": wins_utc})
        print(f"[{tag}] P={P:.4f} d 相位差={out['branches'][tag]['phase_diff']:.3f} "
              f"自洽={out['branches'][tag]['phase_coherent']} → {status} key={key}")
        for a, b in wins[:2]:
            print(f"      窗 {a:.2f} → {b:.2f} BJD")

    print(f"[double] P={BRANCHES['double']:.1f} d 相位差="
          f"{out['branches']['double']['phase_diff']:.3f} "
          f"自洽={out['branches']['double']['phase_coherent']}（谐波反相 → 不注册）")

    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "pair113_period.json").write_text(
        json.dumps(out, indent=2, ensure_ascii=False), encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
