"""Stage 4 acceptance gate G4.

Usage: python scripts/gate_g4.py
"""
from __future__ import annotations

import datetime as dt
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from frbpp.catalogs import chime                              # noqa: E402
from frbpp.exposure import model as em                        # noqa: E402
from frbpp.pointprocess.periodogram import (                  # noqa: E402
    StepWindowScanner, obs_grid_from_obslog)
from frbpp.validation.fpr import fpr_battery_n1, fpr_battery_n2     # noqa: E402
from frbpp.validation.inject_recover import completeness_surface     # noqa: E402

RAW = ROOT / "data" / "raw" / "chime"
OUT = ROOT / "data" / "processed" / "validation"

N_EVENTS = 17
ALPHA_N1, B_N1, N_NULL = 0.01, 400, 400
B_1B, N_NULL_B = 300, 300
ALPHA_N2, B_N2, N_SCR = 0.05, 50, 24
N_REP, N_REP_REF = 20, 40
TOL = 0.01
CONTRAST = 50.0
DT_DAYS = 60.0 / 86400.0


def main() -> int:
    rng_n1 = np.random.default_rng(20261009)
    rng_n1b = np.random.default_rng(20261010)
    rng_n2 = np.random.default_rng(20261011)
    rng_comp = np.random.default_rng(20261012)
    rng_gauss = np.random.default_rng(20261013)
    out = {"built_at": dt.datetime.now().astimezone().isoformat(timespec="seconds"),
           "criteria": "docs/TASKS.md G4 冻结口径 v1.2"}

    mjd, _ = em.load_r3_exposure_series(RAW / "03_180916/R3_CHIME_exposure.npz")
    wins = em.series_to_windows(mjd)
    log = em.obslog_from_series(mjd, 29.503, 65.717)
    grid_t, grid_w = obs_grid_from_obslog(log, dt_days=DT_DAYS)
    toa_df, _ = chime.load_r3_180916_toa(RAW / "03_180916/R3_CHIME_TOA_20Feb4.txt")
    toa = toa_df["mjd_native"].to_numpy()
    events = toa[[any(a <= x <= b for a, b in wins) for x in toa]]
    t0 = float(events[0])
    scanner = StepWindowScanner.build(grid_t, grid_w, t0,
                                      np.linspace(5.0, 35.0, 1200),
                                      np.arange(24) / 24.0,
                                      np.array([0.05, 0.10, 0.15, 0.20, 0.25]))
    print(f"网格 {grid_t.size} 点; 事件基准 N={N_EVENTS}")

    from scipy import stats

    b1 = fpr_battery_n1(scanner, grid_t, grid_w, N_EVENTS, B_N1, N_NULL,
                        ALPHA_N1, rng_n1)
    p_theory = 3.0 / (N_NULL + 1)
    exp01 = B_N1 * p_theory
    sig01 = float(np.sqrt(B_N1 * p_theory * (1 - p_theory)))
    thr01 = int(np.ceil(exp01 + 2 * sig01))
    k05 = int(np.sum(b1["p_values"] < 0.05))
    exp05 = B_N1 * 0.05
    sig05 = float(np.sqrt(B_N1 * 0.05 * 0.95))
    thr05 = int(np.ceil(exp05 + 2 * sig05))
    ks = stats.kstest(b1["p_values"], "uniform")
    ok_n1_main = bool(b1["n_false"] <= thr01 and k05 <= thr05)
    print(f"[FPR-N1] 假阳性 {b1['n_false']}/{b1['n_real']} = {b1['fpr']:.3%} "
          f"(理论期望 {exp01:.1f}±{sig01:.1f}, 判据 ≤{thr01}) | "
          f"p<0.05: {k05}/{b1['n_real']} (期望 {exp05:.0f}±{sig05:.1f}, 判据 ≤{thr05}) "
          f"| KS(报告): p={ks.pvalue:.3f} → {'OK' if ok_n1_main else 'FAIL'}")
    out["fpr_n1"] = {k: (v.tolist() if isinstance(v, np.ndarray) else v)
                     for k, v in b1.items()}
    out["fpr_n1"].update({"expected01": exp01, "sigma01": sig01, "threshold01": thr01,
                          "k05": k05, "expected05": exp05, "threshold05": thr05,
                          "ks_pvalue": float(ks.pvalue), "ok": ok_n1_main})

    b1b = fpr_battery_n1(scanner, grid_t, grid_w, N_EVENTS, B_1B, N_NULL_B,
                         ALPHA_N1, rng_n1b)
    p_theory_b = 3.0 / (N_NULL_B + 1)
    exp01b = B_1B * p_theory_b
    sig01b = float(np.sqrt(B_1B * p_theory_b * (1 - p_theory_b)))
    thr01b = int(np.ceil(exp01b + 2 * sig01b))
    ok_n1_rep = bool(b1b["n_false"] <= thr01b)
    print(f"[FPR-N1 复核] 假阳性 {b1b['n_false']}/{b1b['n_real']} = {b1b['fpr']:.3%} "
          f"(期望 {exp01b:.1f}±{sig01b:.1f}, 判据 ≤{thr01b}) → {'OK' if ok_n1_rep else 'FAIL'}")
    out["fpr_n1_replication"] = {"n_false": b1b["n_false"], "n_real": b1b["n_real"],
                                 "fpr": b1b["fpr"], "expected": exp01b,
                                 "threshold": thr01b, "ok": ok_n1_rep}

    b2 = fpr_battery_n2(scanner, grid_t, grid_w, N_EVENTS, B_N2, N_SCR,
                        ALPHA_N2, rng_n2)
    g4_n2_ok = bool(b2["fpr"] <= 0.10)
    print(f"[FPR-N2] 假阳性 {b2['n_false']}/{b2['n_real']} = {b2['fpr']:.3%} "
          f"(α={ALPHA_N2}; 零分布=加扰后完整扫描 max) → {'OK' if g4_n2_ok else 'FAIL'}")
    out["fpr_n2"] = {k: (v.tolist() if isinstance(v, np.ndarray) else v)
                     for k, v in b2.items()}
    out["fpr_n2"]["ok"] = g4_n2_ok

    ref_configs = [{"P": 16.35, "width": 0.25, "n_ep": n} for n in (10, 20, 40)]
    other_configs = ([{"P": P, "width": w, "n_ep": n}
                      for P in (8.0, 30.0) for w in (0.10, 0.25) for n in (10, 20, 40)]
                     + [{"P": 16.35, "width": 0.10, "n_ep": n} for n in (10, 20, 40)])
    rows = completeness_surface(scanner, grid_t, grid_w, t0, ref_configs, N_REP_REF,
                                CONTRAST, TOL, rng_comp)
    rows += completeness_surface(scanner, grid_t, grid_w, t0, other_configs, N_REP,
                                 CONTRAST, TOL, rng_comp)
    print("[完备度] P      w     n_ep  恢复率  (95% CI)")
    for r in rows:
        print(f"          {r['P']:>6.2f} {r['width']:>5.2f} {r['n_ep']:>5d}  "
              f"{r['recovery']:>5.0%}  [{r['ci95'][0]:.0%},{r['ci95'][1]:.0%}]")
    out["completeness"] = rows

    g_rows = completeness_surface(scanner, grid_t, grid_w, t0,
                                  [{"P": 16.35, "width": 0.25, "n_ep": 20,
                                    "gaussian_sigma": 0.125}], N_REP,
                                  CONTRAST, TOL, rng_gauss)
    print(f"[高斯窗] σ=0.125, P=16.35, n_ep=20: 恢复率 {g_rows[0]['recovery']:.0%}")
    out["gaussian_window"] = g_rows[0]

    ref = [r for r in rows if r["P"] == 16.35 and r["width"] == 0.25 and r["n_ep"] >= 20]
    ref_ok = all(r["recovery"] >= 0.8 for r in ref)
    verdict = bool(ok_n1_main and ok_n1_rep and g4_n2_ok and ref_ok)
    out["verdict"] = "PASS" if verdict else "FAIL"
    out["config"] = {"n_events": N_EVENTS, "alpha_n1": ALPHA_N1, "b_n1": B_N1,
                     "b_n1_replication": B_1B, "alpha_n2": ALPHA_N2, "b_n2": B_N2,
                     "n_surrogate_n2": N_SCR, "n_rep": N_REP, "n_rep_ref": N_REP_REF,
                     "tol": TOL, "contrast": CONTRAST,
                     "reference_cells": "P=16.35, w=0.25, n_ep>=20 (v1.2 校准)"}
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "g4_report.json").write_text(json.dumps(out, indent=2, ensure_ascii=False),
                                        encoding="utf-8")
    print(f"G4: {out['verdict']}  →  data/processed/validation/g4_report.json")
    return 0 if verdict else 1


if __name__ == "__main__":
    raise SystemExit(main())