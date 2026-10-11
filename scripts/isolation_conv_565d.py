"""Localisation-error convolution Monte Carlo for the joint spatial isolation of the primary pair.

Usage: python scripts/isolation_conv_565d.py
Output: data/processed/upgrade/pair565_isolation.json
"""
from __future__ import annotations

import datetime as dt
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from frbpp.upgrade.crossmatch import (build_event_table, geometry_pairs,  # noqa: E402
                                      same_event_mask)

OUT = ROOT / "data" / "processed" / "upgrade"
CAT2 = ROOT / "data" / "interim" / "cat2_bursts.parquet"
CAT1 = ROOT / "data" / "interim" / "cat1_bursts.parquet"

N_MC = 1000
SEED = 20261009
THETA_CONFIG = 600.0
DM_CONFIG = 5.0
THETA_PRIMARY = 23.86
DM_PRIMARY = 0.5
DM_PREFILTER = 20.0


def main() -> int:
    rng = np.random.default_rng(SEED)
    c1 = pd.read_parquet(CAT1)
    c2 = pd.read_parquet(CAT2)
    ev = build_event_table(c1, c2)

    c2p = c2[c2.sub_num == 0].set_index("tns_name")
    c2p = c2p[~c2p.index.duplicated(keep="first")]
    c1p = (c1[c1.is_burst_primary] if "is_burst_primary" in c1 else c1).set_index("tns_name")
    c1p = c1p[~c1p.index.duplicated(keep="first")]
    n_ev = len(ev.ra)
    names = np.asarray(ev.name, dtype=object)
    cat = np.asarray(ev.catalog)
    oneoff = np.array([not s for s in ev.src_rep])
    keep = oneoff & np.isfinite(ev.ra) & np.isfinite(ev.dm)

    sig_ra = np.full(n_ev, np.nan)
    sig_dec = np.full(n_ev, np.nan)
    sig_dm = np.full(n_ev, np.nan)
    n_from_c2 = n_from_c1 = n_missing = 0
    med_ra = np.nanmedian(np.concatenate([c2.ra_err.dropna().to_numpy(float), c1.ra_err.dropna().to_numpy(float)]))
    med_dec = np.nanmedian(np.concatenate([c2.dec_err.dropna().to_numpy(float), c1.dec_err.dropna().to_numpy(float)]))
    med_dm = np.nanmedian(c2.dm_fitb_err.dropna().to_numpy(float))
    for i in range(n_ev):
        if not keep[i]:
            continue
        row = None
        if names[i] in c2p.index:
            row = c2p.loc[names[i]]
            src = "c2"
        elif names[i] in c1p.index:
            row = c1p.loc[names[i]]
            src = "c1"
        re_ = float(row.ra_err) if row is not None and pd.notna(row.ra_err) else np.nan
        de_ = float(row.dec_err) if row is not None and pd.notna(row.dec_err) else np.nan
        de_dm = float(row.dm_fitb_err) if row is not None and pd.notna(row.dm_fitb_err) else np.nan
        if np.isnan(re_):
            re_, de_, de_dm = med_ra, med_dec, med_dm
            n_missing += 1
        elif src == "c2":
            n_from_c2 += 1
        else:
            n_from_c1 += 1
        sig_ra[i], sig_dec[i], sig_dm[i] = re_, de_, de_dm if np.isfinite(de_dm) else med_dm

    idx = np.where(keep)[0]
    ra = np.radians(ev.ra[idx]); dec = np.radians(ev.dec[idx]); dm = np.asarray(ev.dm, float)[idx]
    sra = sig_ra[idx]; sdec = sig_dec[idx]; sdm = sig_dm[idx]
    print(f"one-off 事件（有限 DM+坐标）{idx.size}；误差回联：c2={n_from_c2} c1={n_from_c1} 缺→中位={n_missing}；"
          f"σ_ra 中位={np.median(sra)*3600:.0f}″ σ_dec 中位={np.median(sdec)*3600:.0f}″")

    gi, gj = geometry_pairs(ev, 2000.0)
    keep_pair = ~same_event_mask(ev, gi, gj)
    gi, gj = gi[keep_pair], gj[keep_pair]
    oo_pair = np.array([keep[a] and keep[b] for a, b in zip(gi, gj)])
    gi, gj = gi[oo_pair], gj[oo_pair]
    pos = {g: k for k, g in enumerate(idx)}
    ii = np.array([pos[a] for a in gi]); jj = np.array([pos[b] for b in gj])
    ddm0 = np.abs(dm[ii] - dm[jj])
    m = ddm0 < DM_PREFILTER
    ii, jj = ii[m], jj[m]
    ddm0 = ddm0[m]
    print(f"规范配对表（θ<2000″、去重、双侧 one-off、|ΔDM|<20）：{ii.size} 对")

    def sep_arcsec(ra1, dec1, ra2, dec2):
        x = np.sin(dec1) * np.sin(dec2) + np.cos(dec1) * np.cos(dec2) * np.cos(ra1 - ra2)
        return np.degrees(np.arccos(np.clip(x, -1, 1))) * 3600.0

    th0 = sep_arcsec(ra[ii], dec[ii], ra[jj], dec[jj])
    ddm0v = ddm0
    cfg0 = int(np.sum((th0 < THETA_CONFIG) & (ddm0v < DM_CONFIG)))
    prim0 = int(np.sum((th0 < THETA_PRIMARY) & (ddm0v < DM_PRIMARY)))
    minth5_0 = float(th0[ddm0v < DM_CONFIG].min()) if np.any(ddm0v < DM_CONFIG) else None
    minth05_0 = float(th0[ddm0v < DM_PRIMARY].min()) if np.any(ddm0v < DM_PRIMARY) else None
    mindm_0 = float(ddm0v[th0 < THETA_CONFIG].min()) if np.any(th0 < THETA_CONFIG) else None

    cnt_cfg = np.empty(N_MC); cnt_prim = np.empty(N_MC)
    minth5 = np.full(N_MC, np.inf); minth05 = np.full(N_MC, np.inf); mindm = np.full(N_MC, np.inf)
    for r in range(N_MC):
        ra_p = ra + rng.normal(0, sra)
        dec_p = dec + rng.normal(0, sdec)
        dm_p = dm + rng.normal(0, sdm)
        th = sep_arcsec(ra_p[ii], dec_p[ii], ra_p[jj], dec_p[jj])
        dd = np.abs(dm_p[ii] - dm_p[jj])
        cfg = (th < THETA_CONFIG) & (dd < DM_CONFIG)
        cnt_cfg[r] = int(cfg.sum())
        prim = (th < THETA_PRIMARY) & (dd < DM_PRIMARY)
        cnt_prim[r] = int(prim.sum())
        if np.any(dd < DM_CONFIG):
            minth5[r] = th[dd < DM_CONFIG].min()
        if np.any(dd < DM_PRIMARY):
            minth05[r] = th[dd < DM_PRIMARY].min()
        if np.any(th < THETA_CONFIG):
            mindm[r] = dd[th < THETA_CONFIG].min()

    def stats(x):
        return {"mean": float(np.mean(x)), "p95": float(np.quantile(x, .95)),
                "max": float(np.max(x)), "median": float(np.median(x))}

    p_cfg = float(np.mean(cnt_cfg >= cfg0))
    p_prim = float(np.mean(cnt_prim >= prim0))
    p_minth5 = float(np.mean(minth5 <= minth5_0)) if minth5_0 else None
    p_minth05 = float(np.mean(minth05 <= minth05_0)) if minth05_0 else None
    p_mindm = float(np.mean(mindm <= mindm_0)) if mindm_0 else None
    p_prim_atleast1 = float(np.mean(cnt_prim >= 1))

    out = {
        "built_at": dt.datetime.now().astimezone().isoformat(timespec="seconds"),
        "config": {"n_mc": N_MC, "seed": SEED,
                   "error_model": "每事件独立 1σ 高斯（ra_err/dec_err/dm_fitb_err 按名回联；缺→人口中位）",
                   "assumption": "ra_err/dec_err 按 1σ 处理；两暴同束→系统误差相关 ⇒ 独立误差假设对主候选是保守的",
                   "n_oneoff": int(idx.size), "err_source": {"c2": n_from_c2, "c1": n_from_c1, "median_filled": n_missing},
                   "sigma_ra_median_arcsec": float(np.median(sra) * 3600),
                   "sigma_dec_median_arcsec": float(np.median(sdec) * 3600),
                   "prefilter_pairs": int(ii.size)},
        "observed": {"n_(600,5)": cfg0, "n_(24,0.5)": prim0,
                     "min_theta_(dm<5)_arcsec": minth5_0, "min_theta_(dm<0.5)_arcsec": minth05_0,
                     "min_|dm|(theta<600)": mindm_0},
        "mc": {"n_(600,5)": stats(cnt_cfg), "n_(24,0.5)": stats(cnt_prim),
               "min_theta_(dm<5)": stats(minth5), "min_theta_(dm<0.5)": stats(minth05),
               "min_|dm|(theta<600)": stats(mindm)},
        "p_values": {"P(mc_n(600,5) >= 44)": p_cfg,
                     "P(mc_n(24,0.5) >= 1)": p_prim_atleast1,
                     "P(mc_min_theta(dm<0.5) <= 23.9arcsec)": p_minth05,
                     "P(mc_min_theta(dm<5) <= obs)": p_minth5,
                     "P(mc_min_|dm|(theta<600) <= obs)": p_mindm},
        "samples": {
            "min_theta05_arcsec": [round(float(x), 1) if np.isfinite(x) else 20000.0 for x in minth05],
            "n_cfg": [int(x) for x in cnt_cfg],
            "n_prim": [int(x) for x in cnt_prim]},
        "interpretation": {
            "isolation": ("误差卷积后，机会符合背景中出现主候选级对 (θ<24″ & |ΔDM|<0.5) 的期望次数与概率见 mc/p_values；"
                          "该数即联合空间孤立度的严格量化"),
            "note": ("独立误差假设保守：两暴同束、系统误差高度相关 → 其相对 23.9″ 的真实精度优于 σ_comb≈1189″；"
                     "MC 只会低估孤立度"),
        }}
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "pair565_isolation.json").write_text(
        json.dumps(out, indent=2, ensure_ascii=False), encoding="utf-8")

    print(f"\n观测: n(600,5)={cfg0}  n(24,0.5)={prim0}  minθ(dm<0.5)={minth05_0 and round(minth05_0,1)}″  minDM(θ<600)={mindm_0 and round(mindm_0,3)}")
    print(f"MC  n(600,5): 均值 {np.mean(cnt_cfg):.2f} p95 {np.quantile(cnt_cfg,.95):.0f} max {cnt_cfg.max():.0f} → P(≥obs)={p_cfg}")
    print(f"MC  n(24,0.5): 均值 {np.mean(cnt_prim):.3f} p95 {np.quantile(cnt_prim,.95):.0f} max {cnt_prim.max():.0f} → P(≥1)={p_prim_atleast1}")
    print(f"MC  minθ(dm<0.5): 中位 {np.median(minth05[np.isfinite(minth05)]):.0f}″ → P(≤23.9″)={p_minth05}")
    print(f"MC  minDM(θ<600): 中位 {np.median(mindm[np.isfinite(mindm)]):.2f} → P(≤obs)={p_mindm}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
