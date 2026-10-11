"""Exposure analysis over the silent gap of the 113.69 d pair.

Usage: python scripts/silence_113d.py
Output: data/processed/upgrade/pair113_silence.json
"""
from __future__ import annotations

import datetime as dt
import json
import sys
from pathlib import Path

import numpy as np
from astropy import units as u
from astropy.time import Time
from astropy_healpix import HEALPix

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from frbpp.exposure.model import transit_grid   # noqa: E402

OUT = ROOT / "data" / "processed" / "upgrade"
CAT2 = ROOT / "data" / "interim" / "cat2_bursts.parquet"
EXPO_DIR = ROOT / "data" / "raw" / "chime" / "05_cat1"

NAME1, NAME2 = "FRB20200919A", "FRB20210111A"
MAP_T0, MAP_T1 = 58358.0, 58666.0
T_NOW_MJD = float(Time("2026-10-09T00:00:00", scale="utc").mjd)
SAMPLE_S = 4.0

MJD1, MJD2 = 59111.909114, 59225.596722
DEAD_T = T_NOW_MJD - MJD2

COVERAGE = [(58324.0, 60202.0),
            (61288.0, T_NOW_MJD)]
BLIND = (60202.0, 61288.0)


def _transits(ra, dec, t0, t1):
    up = transit_grid(ra, dec, t0, t1)
    lo = transit_grid((ra + 180.0) % 360.0, dec, t0, t1)
    return int(up.size), int(lo.size)


def _transits_in_coverage(ra, dec, t0, t1, coverage):
    n = 0
    for a, b in coverage:
        lo, hi = max(a, t0), min(b, t1)
        if hi > lo:
            n += sum(_transits(ra, dec, lo, hi))
    return n


def main() -> int:
    import pandas as pd

    cat = pd.read_parquet(CAT2)
    t = cat.mjd_400.to_numpy(dtype=float)
    span = (float(np.nanmin(t)), float(np.nanmax(t)))

    r1 = cat[cat.tns_name == NAME1].iloc[0]
    r2 = cat[cat.tns_name == NAME2].iloc[0]
    ra = 0.5 * (float(r1.ra) + float(r2.ra))
    dec = 0.5 * (float(r1.dec) + float(r2.dec))
    ra1, dec1 = float(r1.ra), float(r1.dec)

    hpix = HEALPix(nside=4096, order="ring")
    ipix = int(hpix.lonlat_to_healpix(ra1 * u.deg, dec1 * u.deg))
    texp_map_s = 0.0
    per_map = {}
    for tag in ("L", "U"):
        fp = EXPO_DIR / (f"exposure_int_20180828_20190702_transit_{tag}"
                         f"_beam_FWHM-600_res_4s_0.86_arcmin.npz")
        arr = np.load(fp, allow_pickle=True)["exposure"]
        s = float(arr[ipix]) * SAMPLE_S
        per_map[tag] = s
        texp_map_s += s

    n_map = sum(_transits(ra1, dec1, MAP_T0, MAP_T1))
    n_gap = sum(_transits(ra1, dec1, MJD1, MJD2))
    n_full = sum(_transits(ra1, dec1, span[0], span[1]))
    n_after2 = sum(_transits(ra1, dec1, MJD2, T_NOW_MJD))
    n_after2_cov = _transits_in_coverage(ra1, dec1, MJD2, T_NOW_MJD, COVERAGE)
    n_blind = sum(_transits(ra1, dec1, BLIND[0], BLIND[1]))

    per_transit_s = texp_map_s / n_map

    def teff_h(n):
        return per_transit_s * n / 3600.0

    exp_up_h = 0.5 * (float(r1.exp_up) + float(r2.exp_up))
    frac_gap = (MJD2 - MJD1) / (span[1] - span[0])
    teff_gap_B_h = exp_up_h * frac_gap

    teff_full_A_h = teff_h(n_full)
    teff_gap_A_h = teff_h(n_gap)

    n_det = 2
    lam_per_h = n_det / teff_full_A_h
    exp_in_gap = lam_per_h * teff_gap_A_h
    p_zero_in_gap = float(np.exp(-exp_in_gap))
    lam_upper95 = 3.0 / teff_gap_A_h
    exp_after_second = lam_per_h * teff_h(n_after2)
    exp_after_second_cov = lam_per_h * teff_h(n_after2_cov)
    teff_blind_h = teff_h(n_blind)
    exp_blind = lam_per_h * teff_blind_h

    expo_ratio_gap = teff_gap_A_h / teff_full_A_h
    sparse = bool(teff_gap_A_h < 1.0)
    verdict = ("exposure_sufficient_zero_burst" if not sparse
               else "exposure_sparse_silence_artifact")

    out = {
        "built_at": dt.datetime.now().astimezone().isoformat(timespec="seconds"),
        "pair": [NAME1, NAME2],
        "position_mean": {"ra": ra, "dec": dec},
        "position_errors_arcsec": {
            "ra_err_1": float(r1.ra_err) * 3600.0, "dec_err_1": float(r1.dec_err) * 3600.0,
            "ra_err_2": float(r2.ra_err) * 3600.0, "dec_err_2": float(r2.dec_err) * 3600.0,
            "note": "形式定位误差 ~0.25° ≈ 810-1065″，远大于 180″ 锥；绝对定位不可信"},
        "epochs": {"mjd_burst1": MJD1, "mjd_burst2": MJD2,
                   "dt_days": MJD2 - MJD1,
                   "cat2_span_mjd": span, "to_now_days": DEAD_T},
        "exposure_map": {"epoch_mjd": [MAP_T0, MAP_T1], "ipix": ipix,
                         "L_s": per_map["L"], "U_s": per_map["U"],
                         "T_eff_map_s": texp_map_s, "T_eff_map_h": texp_map_s / 3600.0},
        "transit_counts": {"map_epoch": n_map, "gap": n_gap,
                           "full_cat2_span": n_full, "after_burst2_to_now": n_after2,
                           "after_burst2_in_coverage": n_after2_cov},
        "coverage": {"segments_mjd": COVERAGE, "blind_mjd": list(BLIND),
                     "note": "Cat2 止于 2023-09-15、GCN 归档自 2026-09-05 → 中间约 3 年为盲区"},
        "per_transit_eff_s": per_transit_s,
        "T_eff_gap_h": {"methodA_map_extrap": teff_gap_A_h, "methodB_cat2_expup": teff_gap_B_h},
        "T_eff_full_h": teff_full_A_h,
        "exposure_ratio_gap_over_full": expo_ratio_gap,
        "rate_constraint": {
            "lambda_per_hour_stationary": lam_per_h,
            "expected_bursts_in_gap_if_stationary": exp_in_gap,
            "P_zero_given_stationary": p_zero_in_gap,
            "lambda_upper95_per_hour": lam_upper95,
            "expected_after_second_burst_naive": exp_after_second,
            "expected_after_second_burst_coverage_aware": exp_after_second_cov,
            "P_zero_after_second_coverage_aware": float(np.exp(-exp_after_second_cov)),
            "blind_zone": {"mjd": list(BLIND), "n_transits": n_blind,
                           "T_eff_h": teff_blind_h,
                           "expected_bursts_if_observing": exp_blind,
                           "P_at_least_one": 1.0 - float(np.exp(-exp_blind)),
                           "note": "盲区=Cat2 止(2023-09-15)→GCN 归档始(2026-09-05)；"
                                   "此期望假设 CHIME 全程正常观测且逐中天效率纪元无关"},
        },
        "verdict": verdict,
        "caveats": [
            "曝光图为 2018-08→2019-07 纪元，用逐中天有效时长外推到 2020-2021，"
            "假设逐中天效率纪元无关（CHIME 灵敏度/RFI 环境可能已变）→ 系统学未消除",
            "事件覆盖非连续：Cat2 止 2023-09-15、GCN 归档自 2026-09-05 → 中间约 3 年盲区，"
            "『第二暴后至今』的期望暴数必须用覆盖感知口径（本脚本已给）",
            "仅 2 个事件 → 本分析给『活动率约束』，不能单独作为『严格成窗』证据",
            "两次探测的定位高度相关（同仪器同位置），20.1″ 一致不改善绝对定位",
        ],
    }
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "pair113_silence.json").write_text(
        json.dumps(out, indent=2, ensure_ascii=False), encoding="utf-8")

    print(f"位置(均值) RA={ra:.5f} Dec={dec:.5f} | ipix={ipix}")
    print(f"曝光图纪元 T_eff={texp_map_s/3600:.2f} h（L={per_map['L']:.0f}s U={per_map['U']:.0f}s）"
          f" 逐中天={per_transit_s:.1f} s")
    print(f"中天数: 图纪元={n_map} 间隔={n_gap} 全Cat2={n_full} "
          f"第二暴后(墙钟)={n_after2} 第二暴后(覆盖内)={n_after2_cov}")
    print(f"间隔内 T_eff: A(图外推)={teff_gap_A_h:.2f} h | B(Cat2 exp_up 折算)={teff_gap_B_h:.2f} h")
    print(f"全纪元 T_eff(A)={teff_full_A_h:.2f} h ; 间隔占比={expo_ratio_gap:.3f}")
    print(f"稳态率 λ={lam_per_h:.4f} /h ; 间隔期望={exp_in_gap:.2f} ; P(0)={p_zero_in_gap:.3f}")
    print(f"零暴 95% 上限 λ<{lam_upper95:.4f} /h")
    print(f"第二暴后期望: 朴素(含盲区)={exp_after_second:.2f} | "
          f"覆盖感知={exp_after_second_cov:.2f} (P0={np.exp(-exp_after_second_cov):.3f})")
    print(f"盲区(2023-09-15→2026-09-05): 中天数={n_blind} T_eff={teff_blind_h:.1f} h → "
          f"若在观测则期望 {exp_blind:.2f} 次 (P(≥1)={1 - np.exp(-exp_blind):.3f})")
    print(f"判定: {verdict}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
