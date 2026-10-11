"""Stage 2 acceptance gate G2.

Usage: python scripts/gate_g2.py
"""
from __future__ import annotations

import datetime as dt
import json
import subprocess
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from astropy import units as u                      # noqa: E402
from astropy.time import Time                       # noqa: E402
from astropy_healpix import HEALPix                 # noqa: E402

from frbpp.catalogs import chime                    # noqa: E402
from frbpp.exposure import model as em              # noqa: E402

RAW = ROOT / "data" / "raw" / "chime"
R3 = RAW / "03_180916" / "R3_CHIME_exposure.npz"
TOA = RAW / "03_180916" / "R3_CHIME_TOA_20Feb4.txt"
LMAP = RAW / "05_cat1" / ("exposure_int_20180828_20190702_transit_L_beam_"
                          "FWHM-600_res_4s_0.86_arcmin.npz")
UMAP = RAW / "05_cat1" / ("exposure_int_20180828_20190702_transit_U_beam_"
                          "FWHM-600_res_4s_0.86_arcmin.npz")

P_DAYS = 16.35
RA_180916, DEC_180916 = 29.503, 65.717
N_BINS = 20
EXPECTED_CLASS = {"in_window": 17, "gap_or_edge": 9, "beyond_range": 12}


def main() -> int:
    out = {"built_at": dt.datetime.now().astimezone().isoformat(timespec="seconds"),
           "criteria": "docs/TASKS.md G2 冻结口径 v1", "P_days": P_DAYS,
           "radec": [RA_180916, DEC_180916]}

    mjd, meta = em.load_r3_exposure_series(R3)
    wins = em.series_to_windows(mjd)
    toa_df, _ = chime.load_r3_180916_toa(TOA)
    toa = toa_df["mjd_native"].to_numpy()
    t0 = meta["t_start_mjd_utc"]

    def classify(x):
        for a, b in wins:
            if a <= x <= b:
                return "in_window"
        return "gap_or_edge" if mjd.min() <= x <= mjd.max() else "beyond_range"

    cls = np.array([classify(x) for x in toa])
    counts = {k: int((cls == k).sum()) for k in EXPECTED_CLASS}
    g2a_ok = counts == EXPECTED_CLASS
    log = em.obslog_from_series(mjd, RA_180916, DEC_180916)
    eps = em.ExposureModel(log).epsilon(toa_df["t_bjd_tdb"].to_numpy())
    eps_ok = bool(np.all(eps[cls == "in_window"] > 0)) and \
        bool(np.all(eps[cls != "in_window"] == 0))
    print(f"[G2a] 分类 {counts}（期望 {EXPECTED_CLASS}）→ {'OK' if g2a_ok else 'FAIL'}")
    print(f"[G2a] 窗内 ε>0 且 其余 ε=0: {eps_ok}")
    out["G2a"] = {"counts": counts, "expected": EXPECTED_CLASS,
                  "class_ok": g2a_ok, "epsilon_guard_ok": eps_ok}

    map_end_off = float(Time("2019-07-02T00:00:00", scale="utc").mjd - t0)
    sel = (mjd - t0) <= map_end_off
    t_series_s = float(sel.sum()) * meta["sample_stride_median_s"]
    hpix = HEALPix(nside=4096, order="ring")
    ipix = int(hpix.lonlat_to_healpix(RA_180916 * u.deg, DEC_180916 * u.deg))
    t_map_s, maps = 0.0, []
    for name, p in (("L", LMAP), ("U", UMAP)):
        if p.exists():
            arr = np.asarray(np.load(p, allow_pickle=True)["exposure"], dtype=float)
            if ipix < arr.size:
                t_map_s += float(arr[ipix]) * 4.0
                maps.append(name)
            del arr
    ratio = float(t_map_s / t_series_s) if t_series_s > 0 else float("nan")
    g2b_ok = bool(0.5 <= ratio <= 2.0)
    print(f"[G2b] T_eff 序列={t_series_s / 3600:.1f} h | 图({'/'.join(maps)})={t_map_s / 3600:.1f} h "
          f"| 比值={ratio:.3f} → {'OK' if g2b_ok else 'FAIL'}（判据 [0.5, 2]）")
    out["G2b"] = {"t_eff_series_hours": t_series_s / 3600.0,
                  "t_eff_map_hours": t_map_s / 3600.0, "maps": maps,
                  "ratio": ratio, "ok": g2b_ok,
                  "map_end_offset_days": map_end_off, "ipix": ipix}

    ph_w = np.mod((mjd - toa[0]) / P_DAYS, 1.0)
    cov, _ = np.histogram(ph_w, bins=N_BINS, range=(0.0, 1.0))
    ph_t = np.mod((toa - toa[0]) / P_DAYS, 1.0)
    bin_t = np.clip((ph_t * N_BINS).astype(int), 0, N_BINS - 1)
    zero_bins = int((cov == 0).sum())
    inwin_zero = int(sum(1 for k, b in zip(cls, bin_t) if k == "in_window" and cov[b] == 0))
    g2c_ok = inwin_zero == 0
    print(f"[G2c] 覆盖箱 min/max = {cov.min()}/{cov.max()}；空箱 {zero_bins}；"
          f"窗内 TOA 落空箱 {inwin_zero} → {'OK' if g2c_ok else 'FAIL'}")
    out["G2c"] = {"coverage_hist": cov.tolist(), "zero_bins": zero_bins,
                  "in_window_toa_in_zero_bins": inwin_zero, "ok": g2c_ok,
                  "note": "覆盖曲线(=曝光在相位上的分布)供 G3 复用"}

    verdict = bool(g2a_ok and eps_ok and g2b_ok and g2c_ok)
    out["verdict"] = "PASS" if verdict else "FAIL"
    try:
        sha = subprocess.run(["git", "-C", str(ROOT), "rev-parse", "HEAD"],
                             capture_output=True, text=True, timeout=10)
        out["git"] = sha.stdout.strip() if sha.returncode == 0 else "no-git"
    except Exception:  # noqa: BLE001
        out["git"] = "no-git"

    dst = ROOT / "data" / "processed" / "exposure"
    dst.mkdir(parents=True, exist_ok=True)
    (dst / "g2_report.json").write_text(json.dumps(out, indent=2, ensure_ascii=False),
                                        encoding="utf-8")
    print(f"G2: {out['verdict']}  →  data/processed/exposure/g2_report.json")
    return 0 if verdict else 1


if __name__ == "__main__":
    raise SystemExit(main())