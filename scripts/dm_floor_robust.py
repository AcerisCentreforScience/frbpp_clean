#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Bootstrap uncertainty and window sensitivity for the DM error floor; appendix material.

Usage: python scripts/dm_floor_robust.py [--boot N] [--seed S]
Output: data/processed/upgrade/c2_dm_floor_bins_v1.csv, c2_dm_floor_robust_v1.json
"""
from __future__ import annotations

import argparse
import datetime as dt
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw" / "chime" / "cat2_25" / "chimefrbcat2.csv"
RANKED = ROOT / "data/processed/upgrade/c2_ranked_crossday_v14.parquet"
OUTD = ROOT / "data/processed/upgrade"
FRACS = (0.15, 0.25, 0.40)
MIN_N = 60
PAIR = ("FRB20200511C", "FRB20211128E")
DM_TARGET = 1834.952
DDM = 0.17688324101300168
BINS = [0, 100, 200, 300, 400, 500, 700, 900, 1100, 1400, 1800, 2500, 4000]


def mad_scale(d: np.ndarray) -> float:
    return float(1.4826 * np.median(np.abs(d - np.median(d))))


def local_window(dm_ref: np.ndarray, target: float, frac: float):
    f = frac
    for _ in range(6):
        m = np.abs(dm_ref - target) <= f * target
        if m.sum() >= MIN_N:
            return m, f
        f *= 1.5
    m = np.abs(dm_ref - target) <= f * target
    return m, f


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--boot", type=int, default=2000)
    ap.add_argument("--seed", type=int, default=20261011)
    a = ap.parse_args()
    rng = np.random.default_rng(a.seed)

    raw = pd.read_csv(RAW, low_memory=False)
    b = pd.to_numeric(raw.bonsai_dm, errors="coerce").to_numpy(float)
    f = pd.to_numeric(raw.dm_fitb, errors="coerce").to_numpy(float)
    ok = np.isfinite(b) & np.isfinite(f)
    dm_ref, D_ref = f[ok], (f - b)[ok]
    print(f"可比事件 n={ok.sum()}（bonsai_dm 非空 {np.isfinite(b).sum()}，"
          f"dm_fitb 非空 {np.isfinite(f).sum()}，两者同时 {ok.sum()}）")

    e_all = pd.to_numeric(raw.dm_fitb_err, errors="coerce").to_numpy(float)
    e_ref = e_all[ok]
    rows = []
    rom, mor = [], []
    for lo, hi in zip(BINS[:-1], BINS[1:]):
        m = (dm_ref >= lo) & (dm_ref < hi)
        n = int(m.sum())
        if n == 0:
            continue
        d = D_ref[m]
        s = mad_scale(d)
        em = e_ref[m]
        k = np.isfinite(em) & (em > 0)
        med_err = float(np.median(em[k])) if k.sum() else float("nan")
        ratio_med = float(np.median(np.abs(d)) / med_err) if k.sum() else float("nan")
        med_ratio = float(np.median(np.abs(d[k]) / em[k])) if k.sum() else float("nan")
        rom.append(ratio_med); mor.append(med_ratio)
        rows.append({"dm_lo": lo, "dm_hi": hi, "dm_center": round((lo + hi) / 2, 1),
                     "n": n, "median_abs_D": round(float(np.median(np.abs(d))), 4),
                     "MAD": round(float(np.median(np.abs(d - np.median(d)))), 4),
                     "s_floor": round(s, 4), "sigma_comb_sqrt2_s": round(np.sqrt(2) * s, 4),
                     "median_formal_err": round(med_err, 4),
                     "ratio_of_medians": round(ratio_med, 1),
                     "median_of_ratios": round(med_ratio, 1)})
    bins = pd.DataFrame(rows)
    bins.to_csv(OUTD / "c2_dm_floor_bins_v1.csv", index=False)
    print(f"分箱表 {len(bins)} 行 -> c2_dm_floor_bins_v1.csv")
    ke = np.isfinite(e_ref) & (e_ref > 0)
    optimism = {"n_with_err": int(ke.sum()),
                "frac_absD_gt_formal_err": round(float((np.abs(D_ref[ke]) > e_ref[ke]).mean()), 4),
                "ratio_of_medians_range": [round(float(np.nanmin(rom)), 1), round(float(np.nanmax(rom)), 1)],
                "median_of_ratios_range": [round(float(np.nanmin(mor)), 1), round(float(np.nanmax(mor)), 1)],
                "overall_median_of_ratios": round(float(np.median(np.abs(D_ref[ke]) / e_ref[ke])), 1),
                "note": ("正文原写'每个 DM 分箱里形式误差都偏小 30--45 倍'，分箱实测不支持：最"
                         "高约 47 倍（700--900），最低的箱（200--300，n=974）只有约 8 倍。故此处"
                         "并列给出两种定义（中位数之比、逐事件比的中位数）与全表分位数。")}

    out = {"built_at": dt.datetime.now().astimezone().isoformat(timespec="seconds"),
           "n_comparable": int(ok.sum()), "reference_pair": {"names": list(PAIR),
                                                            "dm_mean": DM_TARGET, "ddm": DDM},
           "bins_csv": "c2_dm_floor_bins_v1.csv", "fracs": list(FRACS),
           "formal_error_optimism": optimism,
           "min_local_n": MIN_N, "logic_convention": "s=1.4826·MAD(D) 局部窗；σ_comb=√2·s 为保守约定"}
    for fr in FRACS:
        m, used = local_window(dm_ref, DM_TARGET, fr)
        s_hat = mad_scale(D_ref[m])
        idx = np.where(m)[0]
        bs = np.empty(a.boot)
        for k in range(a.boot):
            pick = rng.choice(idx, size=idx.size, replace=True)
            bs[k] = mad_scale(D_ref[pick])
        out[f"window_frac_{fr}"] = {
            "used_frac": round(used, 4), "n_local": int(m.sum()),
            "s_hat": round(s_hat, 4), "sigma_comb_hat": round(np.sqrt(2) * s_hat, 4),
            "s_boot_ci95": [round(float(np.percentile(bs, 2.5)), 4),
                            round(float(np.percentile(bs, 97.5)), 4)],
            "ratio_loose_DMs": round(DDM / s_hat, 4),
            "ratio_loose_ci95": [round(float(DDM / np.percentile(bs, 97.5)), 4),
                                 round(float(DDM / np.percentile(bs, 2.5)), 4)],
            "ratio_conservative_DMsqrt2": round(DDM / (np.sqrt(2) * s_hat), 4),
            "ratio_conservative_ci95": [round(float(DDM / (np.sqrt(2) * np.percentile(bs, 97.5))), 4),
                                        round(float(DDM / (np.sqrt(2) * np.percentile(bs, 2.5))), 4)],
            "boot_draws": a.boot, "seed": a.seed}
        print(f"  窗 ±{fr}: n={m.sum():4d} s={s_hat:.4f}  "
              f"CI95 [{out[f'window_frac_{fr}']['s_boot_ci95'][0]}, "
              f"{out[f'window_frac_{fr}']['s_boot_ci95'][1]}]  "
              f"|ΔDM|/(√2 s)={DDM/(np.sqrt(2)*s_hat):.4f} "
              f"CI95 {out[f'window_frac_{fr}']['ratio_conservative_ci95']}")

    rk = pd.read_parquet(RANKED)
    stab = {}
    for conv in ("r_formal", "rA", "rB"):
        for fr in (15, 25, 40):
            col = conv if conv == "r_formal" else f"{conv}_F{fr}"
            if col not in rk.columns:
                continue
            order = rk.sort_values(col)
            stab[f"{conv}_F{fr}"] = [{"rank": i + 1, "pair": f"{r.name1} x {r.name2}",
                                      "ratio": round(float(getattr(r, col)), 4),
                                      "dt_days": round(float(r.dt_days), 2)}
                                     for i, r in enumerate(order.head(3).itertuples())]
    out["top3_by_window_and_convention"] = stab
    ref_key = "rB_F25"
    first_rank_stable = all(v[0]["pair"] == stab[ref_key][0]["pair"] for v in stab.values())
    top3_sets = {k: {e["pair"] for e in v} for k, v in stab.items()}
    formal_top3 = {e["pair"] for e in stab.get("r_formal_F25", [])}
    out["stability"] = {
        "rank1_identical_in_all_cells": bool(first_rank_stable),
        "top3_formal_vs_conservative": {"formal": sorted(formal_top3),
                                        "conservative": sorted(top3_sets.get("rB_F25", set()))},
        "top3_changes_with_sigma_convention": bool(formal_top3 != top3_sets.get("rB_F25")),
        "top3_identical_within_conservative_conv":
            bool(top3_sets.get("rB_F15") == top3_sets.get("rB_F25") == top3_sets.get("rB_F40")),
        "top3_identical_within_loose_conv":
            bool(top3_sets.get("rA_F15") == top3_sets.get("rA_F25") == top3_sets.get("rA_F40")),
        "note": ("只有第一名在所有 6 个格子里不变；前三名的**组成**随 σ 约定变化"
                 "（形式误差下第二名的一个成员报了被高估的逐事件误差），"
                 "所以正文只主张第一名，其余按无序集合处理。")}

    out["catalogue_consistency"] = {
        "file": str(RAW.relative_to(ROOT)),
        "both_columns_present": bool(np.isfinite(b).any() and np.isfinite(f).any()),
        "n_bonsai": int(np.isfinite(b).sum()), "n_fitb": int(np.isfinite(f).sum()),
        "median_D": round(float(np.median(D_ref)), 4),
        "frac_D_gt_formal_err_median": None,
        "checkable": "同一发布文件内两列同时存在、可比样本与 D 的分布",
        "not_checkable": ("跨 release 的 bonsai_dm/dm_fitb 定义是否变过，无法从单个文件核实；"
                          "本文所有地板统计只基于这一份文件，换 release 需重算")}
    err = pd.to_numeric(raw.dm_fitb_err, errors="coerce").to_numpy(float)[ok]
    both = np.isfinite(err) & np.isfinite(D_ref)
    if both.sum():
        out["catalogue_consistency"]["frac_D_gt_formal_err_median"] = round(
            float((np.abs(D_ref[both]) > np.median(err[both])).mean()), 4)

    (OUTD / "c2_dm_floor_robust_v1.json").write_text(
        json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")
    print(f"写入 c2_dm_floor_robust_v1.json")
    print(f"  稳定性：第一名全格子一致={out['stability']['rank1_identical_in_all_cells']} "
          f"保守约定前三一致={out['stability']['top3_identical_within_conservative_conv']} "
          f"宽松约定前三一致={out['stability']['top3_identical_within_loose_conv']}  "
          f"形式误差前三与地板不同={out['stability']['top3_changes_with_sigma_convention']}")
    return 0


if __name__ == "__main__":
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    raise SystemExit(main())
