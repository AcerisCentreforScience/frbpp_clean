"""Stage 7 acceptance gate G7.

Usage: python scripts/gate_g7.py
Output: data/processed/analysis/g7_report.json
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

from frbpp.analysis.dm import fit_dm_models                       # noqa: E402
from frbpp.upgrade.crossmatch import clean_rep                    # noqa: E402

RAW2 = ROOT / "data" / "interim" / "cat2_bursts.parquet"
RAW1 = ROOT / "data" / "interim" / "cat1_bursts.parquet"
SCAN = ROOT / "data" / "processed" / "census" / "a_layer_scan.parquet"
OUT = ROOT / "data" / "processed" / "analysis"

N_REP = 20
P_SYN = 16.35


def _num(s):
    return pd.to_numeric(s, errors="coerce").to_numpy(dtype=float)


def build_source_dm(cat2_df):
    df = cat2_df.copy()
    rep = clean_rep(df.repeater_name)
    src = np.where(np.array([r != "" for r in rep]), rep, df.tns_name.astype(str))
    mjd, dm, err = _num(df.mjd_400), _num(df.dm_fitb), _num(df.dm_fitb_err)
    sub = _num(df.sub_num)
    ok = np.isfinite(mjd) & np.isfinite(dm) & np.isfinite(err) & (err > 0)
    d = pd.DataFrame({"src": src[ok], "mjd": mjd[ok], "dm": dm[ok], "err": err[ok],
                      "pri": np.nan_to_num(sub)[ok] == 0})
    d["b"] = (d.mjd / 1e-6).round().astype("int64")
    d = (d.sort_values(["src", "b", "pri"], ascending=[True, True, False])
          .drop_duplicates(["src", "b"], keep="first"))
    out = {s: (g.mjd.to_numpy(), g.dm.to_numpy(), g.err.to_numpy())
           for s, g in d.groupby("src")}
    return out, int(np.sum(~ok))


def delta_cat(cat1_df, cat2_df):
    a = pd.DataFrame({"tns": cat1_df.tns_name.astype(str),
                      "mjd": _num(cat1_df.mjd_400), "dm": _num(cat1_df.dm_fitb),
                      "err": _num(cat1_df.dm_fitb_err)})
    a = a[cat1_df.is_burst_primary.to_numpy()]
    b = pd.DataFrame({"tns": cat2_df.tns_name.astype(str),
                      "mjd": _num(cat2_df.mjd_400), "dm": _num(cat2_df.dm_fitb),
                      "err": _num(cat2_df.dm_fitb_err)})
    m = a.merge(b, on="tns", suffixes=("_1", "_2"))
    m = m[np.isfinite(m.mjd_1) & np.isfinite(m.mjd_2) & (np.abs(m.mjd_1 - m.mjd_2) < 0.01)]
    m = (m.assign(dabs=np.abs(m.mjd_1 - m.mjd_2))
          .sort_values("dabs").drop_duplicates("tns", keep="first"))
    ddm = (m.dm_2 - m.dm_1).to_numpy(dtype=float)
    ec = np.sqrt(m.err_1.to_numpy(dtype=float) ** 2 + m.err_2.to_numpy(dtype=float) ** 2)
    ok = np.isfinite(ddm) & np.isfinite(ec)
    ddm, ec = ddm[ok], ec[ok]
    w = 1.0 / ec
    wmed = float(np.average(ddm, weights=w)) if ddm.size else np.nan
    return {"n_pairs": int(ddm.size), "wmean_ddm": wmed,
            "median_ddm": float(np.median(ddm)) if ddm.size else np.nan,
            "mad_ddm": float(1.4826 * np.median(np.abs(ddm - np.median(ddm)))) if ddm.size else np.nan,
            "median_combined_err": float(np.median(ec)) if ddm.size else np.nan}


def synthetic_sinusoid(rng, n=40, span=400.0, A=0.04, jitter=0.02, P=P_SYN):
    t = np.sort(rng.uniform(0, span, size=n))
    phi = rng.uniform(0, 2 * np.pi)
    dm = 350.0 + A * np.sin(2 * np.pi * t / P + phi) + rng.normal(0, jitter, n)
    formal = np.full(n, jitter / 3.0)
    return t, dm, formal


def era_offset_control(rng, P=P_SYN, n=40, delta_d=1.0 / 24, sigma_d=0.3 / 24,
                       band=(10, 30)):
    k = np.arange(n)
    era = ((k >= band[0]) & (k < band[1])).astype(float)
    t = k * P + era * delta_d + rng.normal(0, sigma_d, k.size)
    w = 1.0 / sigma_d
    wv = np.full(k.size, w)
    g1, cov1 = np.polyfit(k, t, 2, w=wv, cov=True)
    z1 = abs(g1[0]) / np.sqrt(cov1[0, 0])
    X = np.column_stack([np.ones_like(k, dtype=float), k, k.astype(float) ** 2, era])
    Xw = X * w
    beta, *_ = np.linalg.lstsq(Xw, t * w, rcond=None)
    cov2 = np.linalg.inv(Xw.T @ Xw)
    z2 = abs(beta[2]) / np.sqrt(cov2[2, 2])
    d_hat = beta[3]
    return {"fake_z": float(z1), "fake_reproduced": bool(z1 >= 3.0),
            "absorbed_c_z": float(z2), "delta_hat_days": float(d_hat),
            "absorbed_ok": bool(z2 <= 3.0 and abs(d_hat - delta_d) <= 0.3 * delta_d)}


def main() -> int:
    rng = np.random.default_rng(20261013)
    out = {"built_at": dt.datetime.now().astimezone().isoformat(timespec="seconds"),
           "criteria": "docs/TASKS.md G7 冻结口径 v1(+v1.1)"}

    rec, amp_ratio = 0, []
    for _ in range(N_REP):
        t, dm, formal = synthetic_sinusoid(rng)
        res = fit_dm_models(t, dm, formal, p_burst=P_SYN, rng=rng, n_perm=300)
        s = res["sinusoid"]
        if s["significant"] and 0.5 * 0.04 < s["amp"] < 1.5 * 0.04:
            rec += 1
        if np.isfinite(res["sigma_jit"]) and res["sigma_jit"] > 0:
            amp_ratio.append(res["sigma_jit"] / float(np.median(formal)))
    ok1 = rec >= 16
    print(f"[G7① 正弦回收] {rec}/20（判 ≥16）| σ_jit/formal 中位 = "
          f"{np.median(amp_ratio):.2f}× → {'OK' if ok1 else 'FAIL'}")
    out["g7_1"] = {"recovered": rec, "n": N_REP, "inflation_median": float(np.median(amp_ratio)),
                   "ok": bool(ok1)}

    fake, absorbed, det = 0, 0, []
    for _ in range(N_REP):
        r = era_offset_control(rng)
        det.append(r)
        fake += r["fake_reproduced"]
        absorbed += r["absorbed_ok"]
    ok2 = fake >= 16 and absorbed >= 16
    print(f"[G7② δ 对照] 假 dP/dt 复现 {fake}/20 | 吸收后恢复 {absorbed}/20（各判 ≥16）"
          f" → {'OK' if ok2 else 'FAIL'}")
    out["g7_2"] = {"fake_reproduced": fake, "absorbed_ok": absorbed, "n": N_REP,
                   "detail": det, "ok": bool(ok2)}

    c1 = pd.read_parquet(RAW1)
    c2 = pd.read_parquet(RAW2)
    dc = delta_cat(c1, c2)
    thresh = 0.5 * dc["median_combined_err"] if np.isfinite(dc["median_combined_err"]) else np.nan
    dc["exceeds_threshold"] = bool(np.isfinite(dc["wmean_ddm"]) and abs(dc["wmean_ddm"]) > thresh)
    print(f"[G7④ δ_cat] 同暴对 {dc['n_pairs']}：加权均值 ΔDM = {dc['wmean_ddm']:.4f} "
          f"(median 误差 {dc['median_combined_err']:.4f}) | 超阈 {dc['exceeds_threshold']}")
    out["delta_cat"] = dc

    src_dm, n_drop = build_source_dm(c2)
    scan = pd.read_parquet(SCAN)
    real, missing, bare, trend_list = [], 0, 0, []
    required = ["sigma_jit", "inflation", "floor_applied", "trend_significant",
                "dbic_M0_M1", "mc_p_trend", "models", "err_scale"]
    for row in scan.itertuples():
        if row.src not in src_dm:
            missing += 1
            continue
        t, dm, err = src_dm[row.src]
        if t.size < 4:
            real.append({"src": row.src, "n": int(t.size), "note": "样本过少（<4 暴）"})
            continue
        res = fit_dm_models(t, dm, err, p_burst=float(row.P_hat), rng=rng, n_perm=300)
        res["src"] = row.src
        res["P_scan"] = float(row.P_hat)
        for k in required:
            if k not in res:
                missing += 1
        if res.get("trend_significant"):
            trio = (np.isfinite(res.get("dbic_M0_M1", np.nan))
                    and np.isfinite(res.get("mc_p_trend", np.nan))
                    and res.get("err_scale") == "jitter-floored")
            bare += int(not trio)
            trend_list.append({"src": row.src, "n": res["n"],
                               "slope_per_day": res["models"]["M1"]["coef_high2low"][0],
                               "dbic_M0_M1": res.get("dbic_M0_M1"),
                               "mc_p": res.get("mc_p_trend"),
                               "sigma_jit": res.get("sigma_jit")})
        real.append(res)
    ok3 = (missing == 0 and bare == 0 and len(real) == len(scan))
    pc = next((r for r in real if r.get("src") == "FRB20180916B"), None)
    print(f"[G7③ 真实输出] {len(real)}/{len(scan)} 源（缺 {missing}、裸趋势 {bare}）→ "
          f"{'OK' if ok3 else 'FAIL'}")
    for tr in sorted(trend_list, key=lambda x: -x["dbic_M0_M1"])[:6]:
        print(f"   [趋势] {tr['src']:14s} n={tr['n']:3d} slope={tr['slope_per_day']:+.5f}/d "
              f"dBIC={tr['dbic_M0_M1']:.1f} p={tr['mc_p']:.4f} σ_jit={tr['sigma_jit']:.3f}")
    if pc:
        tr = pc.get("trend_significant")
        print(f"   180916: n={pc['n']} σ_jit={pc.get('sigma_jit'):.4f} "
              f"infl={pc.get('inflation'):.2f} 趋势显著={tr} "
              f"(dBIC={pc.get('dbic_M0_M1')}, p={pc.get('mc_p_trend')})")
    out["real_sources"] = real
    out["g7_3"] = {"n_sources": len(real), "missing": int(missing), "bare_trend": int(bare),
                   "dropped_bursts_nonfinite": int(n_drop), "ok": bool(ok3)}

    out["rm_status"] = {"status": "deferred-v2",
                        "reason": "Cat1/Cat2 无 RM 列；外部 RM 时间序列未接入（孪生要求保留在账）"}
    ok4 = True
    print(f"[G7④ RM] 状态: deferred-v2（无数据；不阻塞）")

    verdict = bool(ok1 and ok2 and ok3 and ok4)
    out["verdict"] = "PASS" if verdict else "FAIL"
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "g7_report.json").write_text(json.dumps(out, indent=2, ensure_ascii=False,
                                                   default=str), encoding="utf-8")
    print(f"G7: {out['verdict']}  →  data/processed/analysis/g7_report.json")
    return 0 if verdict else 1


if __name__ == "__main__":
    raise SystemExit(main())