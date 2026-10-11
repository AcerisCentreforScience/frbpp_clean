#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Same-beam and correlated-error audit, including hour-angle statistics.

Usage: python scripts/same_beam_test.py --draws 400 --seed 20261011
Output: data/processed/upgrade/same_beam_test.json
"""
from __future__ import annotations

import argparse
import datetime as dt
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import astropy.units as u
from astropy.time import Time

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.stdout.reconfigure(encoding="utf-8", errors="replace")
from frbpp.upgrade.crossmatch import build_event_table  # noqa: E402

OUT = ROOT / "data/processed/upgrade/same_beam_test.json"
CAT1 = ROOT / "data/interim/cat1_bursts.parquet"
CAT2 = ROOT / "data/interim/cat2_bursts.parquet"
SITE_LON_DEG = -119.6286
SITE_LAT_DEG = 49.3167
PAIR = ("FRB20200511C", "FRB20211128E")
THETA_OBS = 23.86
DDM_OBS = 0.1769
THETA_CUT = 24.0
DDM_CUT = 0.18
HA_SAME_BEAM_H = 0.1
RHO_GRID = (0.0, 0.3, 0.6, 0.9, 0.99)


def sep_arcsec(ra1, dec1, ra2, dec2):
    x = np.sin(dec1) * np.sin(dec2) + np.cos(dec1) * np.cos(dec2) * np.cos(ra1 - ra2)
    return np.degrees(np.arccos(np.clip(x, -1, 1))) * 3600.0


def hour_angle_h(mjd_utc, ra_deg):
    t = Time(mjd_utc, format="mjd", scale="utc")
    lst_h = t.sidereal_time("apparent", longitude=SITE_LON_DEG * u.deg).hour
    if np.ndim(lst_h) == 0:
        lst_h = np.full(np.shape(ra_deg), float(lst_h))
    return (np.asarray(lst_h) - np.asarray(ra_deg) / 15.0 + 12.0) % 24.0 - 12.0


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--draws", type=int, default=2000)
    ap.add_argument("--seed", type=int, default=20261011)
    a = ap.parse_args()
    rng = np.random.default_rng(a.seed)

    c1 = pd.read_parquet(CAT1)
    c2 = pd.read_parquet(CAT2)
    ev = build_event_table(c1, c2, event_dedup=True)
    mjd = np.asarray(getattr(ev, "mjd", None), float) if hasattr(ev, "mjd") else None
    names = np.asarray(ev.name)
    ra = np.radians(np.asarray(ev.ra, float))
    dec = np.radians(np.asarray(ev.dec, float))
    dm = np.asarray(ev.dm, float)
    oo = np.array([r == "" for r in ev.src_rep]) & np.isfinite(dm) & np.isfinite(np.asarray(ev.ra, float))
    idx = np.where(oo)[0]
    assert len(idx) == 3665, f"声明域应为 3665，实得 {len(idx)}"
    ra_d, dec_d, dm_d, nm_d = ra[idx], dec[idx], dm[idx], names[idx]
    t_d = mjd[idx]

    ha = hour_angle_h(t_d, np.degrees(ra_d))
    aha = np.abs(ha)
    c2all = c2[c2.ra.notna() & c2.dec.notna()]
    ha_all = np.abs(hour_angle_h(pd.to_numeric(c2all.mjd_400, errors="coerce").to_numpy(float),
                                 pd.to_numeric(c2all.ra, errors="coerce").to_numpy(float)))
    ha_all = ha_all[np.isfinite(ha_all)]
    b0 = c2[(c2.sub_num == 0) & c2.tns_name.isin(PAIR)].copy()
    h0 = dict(zip(b0.tns_name, hour_angle_h(
        pd.to_numeric(b0.mjd_400, errors="coerce").to_numpy(float),
        pd.to_numeric(b0.ra, errors="coerce").to_numpy(float))))
    pair_dha = float(np.abs(h0[PAIR[0]] - h0[PAIR[1]])) if set(PAIR) <= set(h0) else float("nan")
    rep = c2all["repeater_name"].astype("string").str.strip()
    rok = rep.notna() & ~rep.isin(["", "-9999", "nan", "None"])
    h_rep = hour_angle_h(pd.to_numeric(c2all.mjd_400, errors="coerce").to_numpy(float)[rok.to_numpy()],
                         pd.to_numeric(c2all.ra, errors="coerce").to_numpy(float)[rok.to_numpy()])
    s_rep = rep[rok].to_numpy()
    within = []
    for s in pd.unique(s_rep):
        h = h_rep[s_rep == s]
        if len(h) < 2:
            continue
        i, j = np.triu_indices(len(h), 1)
        dif = np.abs(h[i] - h[j])
        within.append(np.minimum(dif, 24.0 - dif))
    within = np.concatenate(within) if within else np.array([np.nan])
    ha_tdb = np.abs(hour_angle_h(pd.to_numeric(c2all.t_bjd_tdb, errors="coerce").to_numpy(float),
                                 pd.to_numeric(c2all.ra, errors="coerce").to_numpy(float)))
    ha_tdb = ha_tdb[np.isfinite(ha_tdb)]
    ha_stats = {"n": int(len(aha)), "median_abs_ha_h": round(float(np.median(aha)), 4),
                "frac_abs_ha_lt_0p1_h": round(float((aha < 0.1).mean()), 4),
                "frac_abs_ha_lt_0p5_h": round(float((aha < 0.5).mean()), 4),
                "p90_abs_ha_h": round(float(np.percentile(aha, 90)), 4),
                "declared_domain_note": "以上是声明域(3665)口径，同样按 utc",
                "all_cat2_rows": {"n": int(len(ha_all)),
                                  "median_abs_ha_h": round(float(np.median(ha_all)), 4),
                                  "frac_abs_ha_lt_0p1_h": round(float((ha_all < 0.1).mean()), 4),
                                  "note": "§4.3 正文与 Fig 4c 图注的'全目录'数出自此口径"},
                "pair_abs_dha_h": round(pair_dha, 4),
                "same_source_dha": {"median_abs_dha_h": round(float(np.median(within)), 4),
                                    "n_pairs": int(len(within)),
                                    "n_sources": int(pd.unique(s_rep).size),
                                    "definition": "已命名源内部两两 |ΔHA|（折到 ≤12 h）的中位数"},
                "time_column_control": {
                    "median_abs_ha_h_mjd_400_utc": round(float(np.median(ha_all)), 4),
                    "median_abs_ha_h_t_bjd_tdb": round(float(np.median(ha_tdb)), 4),
                    "note": "mjd_400(UTC) 下 |HA| 集中在凌日窗，t_bjd_tdb 下近乎均匀"
                            " ⇒ 赤经反推用的是 mjd_400，不是 BJD_TDB"},
                "definitional_caveat": ("CHIME 为漂扫仪器，Cat2 的赤经由到达时间反推"
                                        "（RA≈该时刻地方恒星时），故 |HA|≈0 与 Δt 落在整数"
                                        "恒星日上属目录构造的结果，不能当作独立证据；本文只把"
                                        "它用作'为什么不能把 565 d 读成物理周期'的守卫。")}

    n = len(ra_d)
    ii, jj = np.triu_indices(n, 1)
    th = sep_arcsec(ra_d[ii], dec_d[ii], ra_d[jj], dec_d[jj])
    ddm = np.abs(dm_d[ii] - dm_d[jj])
    dha = np.abs(hour_angle_h(t_d[ii], np.degrees(ra_d[ii]))
                 - hour_angle_h(t_d[jj], np.degrees(ra_d[jj])))
    dha = (dha + 12.0) % 24.0 - 12.0
    dha = np.abs(dha)
    same = dha < HA_SAME_BEAM_H

    def rate(mask, label):
        k = int(mask.sum())
        p_th = float((th[mask] < THETA_CUT).mean()) if k else float("nan")
        p_both = float(((th[mask] < THETA_CUT) & (ddm[mask] < DDM_CUT)).mean()) if k else float("nan")
        return {"pool": label, "n_pairs": k,
                "P_theta_lt_24": None if k == 0 else round(p_th, 8),
                "P_theta_and_ddm": None if k == 0 else round(p_both, 10),
                "expected_count_theta_and_ddm": None if k == 0 else round(p_both * k, 4),
                "median_theta_arcsec": round(float(np.median(th[mask])), 2) if k else None}

    cross_day = np.abs(t_d[ii] - t_d[jj]) > 1.0
    pools = [rate(np.ones(len(th), bool), "all domain pairs"),
             rate(cross_day, "cross-day pairs (excl. same-night)"),
             rate(same, f"same-beam proxy |dHA|<{HA_SAME_BEAM_H} h"),
             rate(same & cross_day, "same-beam AND cross-day"),
             rate(~same, "different-beam pairs")]

    e1 = c2[(c2.sub_num == 0) & (c2["tns_name"].isin([PAIR[0]]))]
    e2 = c2[(c2.sub_num == 0) & (c2["tns_name"].isin([PAIR[1]]))]
    def sig(row):
        if len(row) == 0:
            return None
        r = row.iloc[0]
        return (float(r.ra_err) * np.cos(np.radians(float(r.dec))) * 3600.0,
                float(r.dec_err) * 3600.0)
    s1, s2 = sig(e1), sig(e2)
    rho_table = []
    if s1 and s2:
        for rho in RHO_GRID:
            sx = np.sqrt(s1[0] ** 2 + s2[0] ** 2 - 2 * rho * s1[0] * s2[0])
            sy = np.sqrt(s1[1] ** 2 + s2[1] ** 2 - 2 * rho * s1[1] * s2[1])
            scomb = np.sqrt(sx ** 2 + sy ** 2) / np.sqrt(2)
            rho_table.append({"rho": rho, "sigma_rel_ra_arcsec": round(float(sx), 1),
                              "sigma_rel_dec_arcsec": round(float(sy), 1),
                              "sigma_rel_comb_arcsec": round(float(scomb), 1),
                              "theta_over_sigma_rel": round(float(THETA_OBS / scomb), 4)})

    look = {}
    for r in c2[c2.sub_num == 0].itertuples():
        if isinstance(r.tns_name, str) and r.tns_name:
            look[r.tns_name] = (float(r.ra_err) if np.isfinite(r.ra_err) else np.nan,
                                float(r.dec_err) if np.isfinite(r.dec_err) else np.nan)
    med_ra = np.nanmedian([v[0] for v in look.values() if np.isfinite(v[0])])
    med_de = np.nanmedian([v[1] for v in look.values() if np.isfinite(v[1])])
    e_ra = np.array([look.get(str(x), (np.nan, np.nan))[0] for x in nm_d])
    e_de = np.array([look.get(str(x), (np.nan, np.nan))[1] for x in nm_d])
    e_ra = np.where(np.isfinite(e_ra), e_ra, med_ra)
    e_de = np.where(np.isfinite(e_de), e_de, med_de)
    sra_rad = np.radians(e_ra) * np.cos(dec_d)
    sde_rad = np.radians(e_de)
    assert np.allclose(dm_d[nm_d == PAIR[0]][0] - dm_d[nm_d == PAIR[1]][0], DDM_OBS, atol=2e-3)

    B = a.draws
    cnt_all = np.empty(B, int)
    cnt_same = np.empty(B, int)
    cnt_cross = np.empty(B, int)
    for b in range(B):
        ra_p = ra_d + rng.normal(0, sra_rad)
        de_p = dec_d + rng.normal(0, sde_rad)
        th_p = sep_arcsec(ra_p[ii], de_p[ii], ra_p[jj], de_p[jj])
        m = (th_p < THETA_CUT) & (ddm < DDM_CUT)
        cnt_all[b] = int(m.sum())
        cnt_same[b] = int((m & same).sum())
        cnt_cross[b] = int((m & cross_day).sum())
    conv = {"draws": B, "seed": a.seed,
            "all_pairs": {"median": float(np.median(cnt_all)), "mean": round(float(cnt_all.mean()), 4),
                          "p95": float(np.percentile(cnt_all, 95)), "max": int(cnt_all.max()),
                          "P_ge_1": round(float((cnt_all >= 1).mean()), 5)},
            "same_beam": {"median": float(np.median(cnt_same)), "mean": round(float(cnt_same.mean()), 4),
                          "p95": float(np.percentile(cnt_same, 95)), "max": int(cnt_same.max()),
                          "P_ge_1": round(float((cnt_same >= 1).mean()), 5)},
            "cross_day": {"median": float(np.median(cnt_cross)), "mean": round(float(cnt_cross.mean()), 5),
                          "p95": float(np.percentile(cnt_cross, 95)), "max": int(cnt_cross.max()),
                          "P_ge_1": round(float((cnt_cross >= 1).mean()), 5)},
            "note": ("承重的限制是**跨日**而不是同束：小 θ 背景对里绝大多数是同夜碎片；"
                     "且本目录里 θ 小的对必然同凌日窗（不同束的对 θ 极大），"
                     "所以'同束'过滤对小 θ 池几乎是空操作——不能拿它当独立条件。")}

    rec = {"built_at": dt.datetime.now().astimezone().isoformat(timespec="seconds"),
           "purpose": "审稿 A3：23.9″ 是否为同束/相关系统误差的自然产物",
           "data_gaps": ["Cat2 无 beam/beamformer/zenith-angle 列，束身份用 |ΔHA|<0.1 h 近似",
                         "Cat2 只有 ra_err/dec_err，无误差相关系数 ⇒ 相关性只能做 ρ 敏感性扫描，不是真协方差",
                         "θ 的联合分布用球面角距；抖动按各事件 1σ 高斯，逐轴独立"],
           "observed_pair": {"names": list(PAIR), "theta_arcsec": THETA_OBS, "ddm": DDM_OBS,
                             "theta_cut": THETA_CUT, "ddm_cut": DDM_CUT},
           "ha_stats": ha_stats,
           "pools": pools,
           "rho_sensitivity": rho_table,
           "error_convolved_counts": conv,
           "verdict": None}
    allp = next(p for p in pools if p["pool"] == "cross-day pairs (excl. same-night)")
    samep = next(p for p in pools if p["pool"] == "same-beam AND cross-day")
    if samep["P_theta_and_ddm"] and allp["P_theta_and_ddm"]:
        lift = samep["P_theta_and_ddm"] / allp["P_theta_and_ddm"]
        rec["verdict"] = {
            "same_beam_lift_vs_all": round(float(lift), 3),
            "pools_compared": {"numerator": samep["pool"], "denominator": allp["pool"]},
            "reading": ("同为跨日池的前提下，再加同束条件后 θ<24″∧|ΔDM|<0.18 的条件率倍数。"
                        "接近 1 ⇒ 小 θ 不是束造成的；远大于 1 ⇒ 审稿人的担心成立。")}
    OUT.write_text(json.dumps(rec, ensure_ascii=False, indent=2), encoding="utf-8")

    print(f"写入 {OUT.name}")
    print(f"|HA| 中位 {ha_stats['median_abs_ha_h']} h，<0.1 h 占 {ha_stats['frac_abs_ha_lt_0p1_h']*100:.1f}%"
          f"  （稿里现值 0.049 h / 78.9%）")
    for p in pools:
        print(f"  {p['pool']:44s} 对数 {p['n_pairs']:>10,}  P(θ<24)={p['P_theta_lt_24']}  "
              f"P(θ&ΔDM)={p['P_theta_and_ddm']}  期望 {p['expected_count_theta_and_ddm']}")
    for r in rho_table:
        print(f"  ρ={r['rho']:<5} σ_rel={r['sigma_rel_comb_arcsec']:>7}\"  θ/σ_rel={r['theta_over_sigma_rel']}")
    print(f"  卷积后：全池 P(≥1)={conv['all_pairs']['P_ge_1']}  同束池 P(≥1)={conv['same_beam']['P_ge_1']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
