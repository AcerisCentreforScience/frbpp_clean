"""Compute and freeze the joint-isolation statistic used by the scatter figure.

Usage: python scripts/build_isolation_data.py
Output: data/processed/upgrade/c2_isolation_figure_v1.json
"""
from __future__ import annotations

import dataclasses
import datetime as dt
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.stdout.reconfigure(encoding="utf-8", errors="replace")
from frbpp.upgrade.crossmatch import build_event_table  # noqa: E402

OUT = ROOT / "data/processed/upgrade/c2_isolation_figure_v1.json"
CAND_PAIR = ("FRB20200511C", "FRB20211128E")
THETA_CUT = 30.0
DDM_CUT = 0.5
CROSS_DAY = 1.0


def main() -> int:
    c1 = pd.read_parquet(ROOT / "data/interim/cat1_bursts.parquet")
    c2 = pd.read_parquet(ROOT / "data/interim/cat2_bursts.parquet")
    ev = build_event_table(c1, c2, event_dedup=True)
    oo = np.array([r == "" for r in ev.src_rep]) & np.isfinite(ev.dm)
    names = np.asarray(ev.name)[oo]
    ra = np.radians(np.asarray(ev.ra)[oo].astype(float))
    dec = np.radians(np.asarray(ev.dec)[oo].astype(float))
    dm = np.asarray(ev.dm)[oo].astype(float)
    mjd = np.asarray(ev.mjd)[oo].astype(float)
    n = len(ra)

    TH = np.full((n, n), np.inf)
    for i0 in range(0, n - 1, 300):
        ii = np.arange(i0, min(i0 + 300, n)); jj = np.arange(i0 + 1, n)
        cs = (np.sin(dec[ii])[:, None] * np.sin(dec[jj])[None, :]
              + np.cos(dec[ii])[:, None] * np.cos(dec[jj])[None, :]
              * np.cos(ra[ii][:, None] - ra[jj][None, :]))
        t = np.degrees(np.arccos(np.clip(cs, -1, 1))) * 3600.0
        TH[np.ix_(ii, jj)] = t
        TH[np.ix_(jj, ii)] = t.T
    np.fill_diagonal(TH, np.inf)

    srt = np.sort(TH, axis=1)
    nn1, nn2 = srt[:, 0], srt[:, 1]
    ki = int(np.where(names == CAND_PAIR[0])[0][0])
    kj = int(np.where(names == CAND_PAIR[1])[0][0])
    assert abs(nn1[ki] - 23.86) < 0.2, f"最近邻与 23.9″ 不符：{nn1[ki]}"

    u = np.triu_indices(n, 1)
    th = TH[u]
    dd = np.abs(dm[u[0]] - dm[u[1]])
    dtt = np.abs(mjd[u[0]] - mjd[u[1]])
    cross = dtt > CROSS_DAY
    box = (th < THETA_CUT) & (dd < DDM_CUT)
    in_box = np.where(box)[0]
    in_box_cross = np.where(box & cross)[0]

    def pair_of(k):
        i, j = int(u[0][k]), int(u[1][k])
        return {"name1": str(names[i]), "name2": str(names[j]),
                "theta_arcsec": round(float(th[k]), 3),
                "ddm": round(float(dd[k]), 4), "dt_days": round(float(dtt[k]), 4)}

    grid = np.geomspace(1.0, 2000.0, 60)
    curve = [{"theta": round(float(g), 3),
              "n_cross": int((cross & (th < g)).sum()),
              "n_cross_dm": int((cross & (th < g) & (dd < DDM_CUT)).sum()),
              "n_all": int((th < g).sum())} for g in grid]

    rec = {
        "built_at": dt.datetime.now().astimezone().isoformat(timespec="seconds"),
        "population": {"n_oneoff_finite_dm": int(n), "n_pairs": int(len(th)),
                       "n_pairs_cross_day": int(cross.sum())},
        "criterion": {"theta_cut_arcsec": THETA_CUT, "ddm_cut": DDM_CUT,
                      "cross_day_d": CROSS_DAY},
        "nearest_neighbour": {
            "member_1nn_arcsec": round(float(nn1[ki]), 2),
            "member_2nn_arcsec": round(float(nn2[ki]), 1),
            "partner_is_nearest": bool(abs(nn1[ki] - TH[ki, kj]) < 1e-6),
            "ratio_2nn_over_pair": round(float(nn2[ki] / nn1[ki]), 1),
            "frac_events_nn_gt_pair": round(float((nn1 > nn1[ki]).mean()), 4),
            "median_1nn_arcsec": round(float(np.median(nn1)), 0)},
        "nn_ecdf": [{"theta": round(float(g), 2),
                     "frac": round(float((nn1 <= g).mean()), 5)}
                    for g in np.geomspace(1.0, 2.0e5, 48)],
        "isolation_box": {
            "n_pairs_in_box_all": int(len(in_box)),
            "n_pairs_in_box_cross_day": int(len(in_box_cross)),
            "members_cross_day": [pair_of(k) for k in in_box_cross]},
        "cross_day_curve": curve,
        "caveat": ("图注里'区内只有两对'仅在跨日域成立；不加 Δt>1 d 时同一框内有 "
                   f"{len(in_box)} 对（多为同夜碎片）。摘要的 4622″ 是次近邻，"
                   "最近邻就是配对伙伴本身。"),
    }
    OUT.write_text(json.dumps(rec, ensure_ascii=False, indent=2), encoding="utf-8")
    print(f"写入 {OUT.name}")
    print(f"  最近邻 {rec['nearest_neighbour']['member_1nn_arcsec']}″（=伙伴）  "
          f"次近邻 {rec['nearest_neighbour']['member_2nn_arcsec']}″  "
          f"比值 {rec['nearest_neighbour']['ratio_2nn_over_pair']}×")
    print(f"  90.1% 检查：{rec['nearest_neighbour']['frac_events_nn_gt_pair']*100:.1f}% 事件最近邻大于本对")
    print(f"  框内：全部 {rec['isolation_box']['n_pairs_in_box_all']} 对 / "
          f"跨日 {rec['isolation_box']['n_pairs_in_box_cross_day']} 对")
    for p in rec["isolation_box"]["members_cross_day"]:
        print(f"     {p['name1']} × {p['name2']}  θ={p['theta_arcsec']}″ "
              f"|ΔDM|={p['ddm']} Δt={p['dt_days']} d")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
