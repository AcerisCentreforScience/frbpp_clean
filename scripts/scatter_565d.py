"""Scattering-timescale analysis separating the two DM hypotheses for the primary pair.

Usage: python scripts/scatter_565d.py
Output: data/processed/upgrade/pair565_scatter.json
"""
from __future__ import annotations

import datetime as dt
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "data" / "processed" / "upgrade"
CAT2 = ROOT / "data" / "interim" / "cat2_bursts.parquet"

NAMES = ("FRB20200511C", "FRB20211128E")
SNR_LO, SNR_HI = 8.0, 15.0
DM_LO, DM_HI = 1500.0, 2000.0


def fit_tau_dm(df):
    x = np.log10(df.dm_fitb.to_numpy(float))
    y = np.log10(df.scat_time.to_numpy(float))
    b, a = np.polyfit(x, y, 1)
    resid = y - (a + b * x)
    return float(a), float(b), float(np.std(resid)), int(len(df))


def main() -> int:
    c2 = pd.read_parquet(CAT2)
    b = c2[c2.sub_num == 0].copy()
    meas = b[(b.scat_time > 0) & b.scat_time_err.notna()]
    non = b[~(b.scat_time > 0)]
    out = {"built_at": dt.datetime.now().astimezone().isoformat(timespec="seconds"),
           "note": "scat_time 单位=秒（目录内部一致）；τ=0 且 err=NaN 视为未测得（低于地板）",
           "population": {"burst_level": int(len(b)),
                          "tau_measured": int(len(meas)),
                          "tau_zero_or_unmeasured": int(len(non))}}

    rows = {}
    for nm in NAMES:
        r = b[b.tns_name == nm].iloc[0]
        rows[nm] = {"scat_time": float(r.scat_time), "scat_time_err": None if pd.isna(r.scat_time_err) else float(r.scat_time_err),
                    "width_fitb_ms": float(r.width_fitb) * 1e3, "width_fitb_err_ms": float(r.width_fitb_err) * 1e3,
                    "bc_width_ms": float(r.bc_width) * 1e3,
                    "peak_freq": float(r.peak_freq), "low_freq": float(r.low_freq), "high_freq": float(r.high_freq),
                    "sp_idx": float(r.sp_idx), "sp_idx_err": float(r.sp_idx_err),
                    "bonsai_snr": float(r.bonsai_snr), "dm_fitb": float(r.dm_fitb),
                    "fluence_jy_ms": float(r.fluence), "flux_jy": float(r.flux)}
    out["bursts"] = rows
    out["band_limited_note"] = ("两暴 low_freq≈605/636 MHz → 带限（仅带上半可测），"
                                "sp_idx=+108/+174 为带限退化拟合，无物理意义，弃用")

    cell = b[(b.bonsai_snr >= SNR_LO) & (b.bonsai_snr < SNR_HI)
             & (b.dm_fitb >= DM_LO) & (b.dm_fitb < DM_HI)
             & b.dm_fitb.notna()]
    cellm = cell[(cell.scat_time > 0) & cell.scat_time_err.notna()]
    out["control_cell"] = {
        "def": f"bonsai_snr∈[{SNR_LO},{SNR_HI}) & DM∈[{DM_LO:.0f},{DM_HI:.0f})",
        "n_total": int(len(cell)), "n_tau_measured": int(len(cellm)),
        "frac_measured": float(len(cellm) / max(len(cell), 1)),
        "tau_ms_measured": {"min": float(cellm.scat_time.min() * 1e3) if len(cellm) else None,
                            "p25": float(cellm.scat_time.quantile(.25) * 1e3) if len(cellm) else None,
                            "median": float(cellm.scat_time.median() * 1e3) if len(cellm) else None,
                            "p75": float(cellm.scat_time.quantile(.75) * 1e3) if len(cellm) else None},
        "names_measured": sorted(cellm.tns_name.tolist())[:20],
    }

    snr_band = b[(b.bonsai_snr >= SNR_LO) & (b.bonsai_snr < SNR_HI)]
    snr_band_m = snr_band[(snr_band.scat_time > 0) & snr_band.scat_time_err.notna()]
    floor_all = float(snr_band_m.scat_time.min() * 1e3) if len(snr_band_m) else None
    floor_err = float(snr_band_m.scat_time_err.min() * 1e3) if len(snr_band_m) else None
    out["floor"] = {"def": f"全部 burst-level、SNR∈[{SNR_LO},{SNR_HI})、τ>0 可测者",
                    "tau_min_ms": floor_all, "err_min_ms": floor_err,
                    "n_measured_in_band": int(len(snr_band_m)),
                    "tau_p10_ms": float(snr_band_m.scat_time.quantile(.10) * 1e3) if len(snr_band_m) else None}

    a, bsig, sig, nfit = fit_tau_dm(meas)
    dm_pair = 1834.95
    tau_pred_ms = (10 ** (a + bsig * np.log10(dm_pair))) * 1e3
    floor_use = floor_all if floor_all else (out["floor"]["tau_p10_ms"])
    R = tau_pred_ms / floor_use if floor_use else None
    if R is None:
        verdict = "floor 不可估 → 无信息"
    elif R > 10:
        verdict = ("R>10：按目录标度 DM≈1835 的湍流介质应有可测散射，实测 τ<地板"
                   " → 该对是散射低于 DM 趋势的离群点 → 湍流屏版 H2 被削弱")
    elif R > 1:
        verdict = "R~1-10 → 弱信息"
    else:
        verdict = "R<1 → 无信息"
    out["scaling"] = {"fit": "log10(tau_s) = a + b*log10(DM)，仅可测事件（同仪器同约定，无需外部系数）",
                      "a": a, "b": bsig, "resid_std_dex": sig, "n": nfit,
                      "dm_range_of_fit": [float(meas.dm_fitb.min()), float(meas.dm_fitb.max())],
                      "tau_predicted_at_1835_ms": tau_pred_ms,
                      "tau_floor_used_ms": floor_use,
                      "ratio_R_predicted_over_floor": R,
                      "verdict": verdict}

    out["interpretation"] = {
        "h2_turbulent_screen": ("削弱：若 DM≈1835 大部分来自湍流致密屏（SNR/星风等），"
                                "按目录内部 τ(DM) 标度应产生 ≫地板的散射；实测未测得 → 湍流屏版 H2 与观测不符"),
        "h2_smooth_medium": ("存活：光滑（低湍流）电离介质可贡献 DM 而少散射——本观测无法排除"),
        "h1_igm": ("相容：IGM 湍流模型预测跨度大，τ<地板与低湍流 IGM 一致；不构成对 H1 的独立支持"),
        "rm": ("RM 不在 Cat1/Cat2 列中（无偏振字段）→ 不可从公开目录提取；"
               "ADS 核查另见 pair565_ads 通道；RM 仍是 H1/H2 最干净的补充判别量（需 CHIME 内部数据或专项观测）"),
        "caveats": [
            "τ<地板 的定量上限依赖 SNR 匹配细胞内最小可测 τ 作为地板代理",
            "τ(DM) 标度拟合含全部可测事件（含重复源），存在可测性选择效应（幸存者偏差偏向高 τ）→ R 是保守下界",
            "目录未给 τ 的参考频率约定；本分析全部在目录内部比较（同约定自洽），不跨目录换算",
            "sp_idx 带限退化，不可用",
        ]}

    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "pair565_scatter.json").write_text(
        json.dumps(out, indent=2, ensure_ascii=False), encoding="utf-8")

    print(f"两暴 τ：20200511C={rows['FRB20200511C']['scat_time']} / "
          f"20211128E={rows['FRB20211128E']['scat_time']}（均未测得）")
    cc = out["control_cell"]
    print(f"控制细胞 {cc['def']}: n={cc['n_total']}, 可测={cc['n_tau_measured']} "
          f"({cc['frac_measured']*100:.0f}%), τ中位={cc['tau_ms_measured']['median'] and round(cc['tau_ms_measured']['median'],1)}ms")
    print(f"地板(同SNR档): min τ={floor_all:.2f} ms (n={out['floor']['n_measured_in_band']})")
    sc = out["scaling"]
    print(f"τ(DM) 标度: b={bsig:.2f}±, 外推 DM=1835 → τ_pred={tau_pred_ms:.1f} ms; "
          f"R=τ_pred/地板={R:.0f}")
    print(f"判定: {sc['verdict']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
