"""Stage 3 acceptance gate G3.

Usage: python scripts/gate_g3.py
"""
from __future__ import annotations

import datetime as dt
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from frbpp.catalogs import chime                       # noqa: E402
from frbpp.exposure import model as em                 # noqa: E402
from frbpp.pointprocess import null as ppnull          # noqa: E402
from frbpp.pointprocess.decluster import decluster     # noqa: E402
from frbpp.pointprocess.periodogram import (           # noqa: E402
    StepWindowScanner, kuiper_statistic, obs_grid_from_obslog)

RAW = ROOT / "data" / "raw" / "chime"
R3 = RAW / "03_180916" / "R3_CHIME_exposure.npz"
TOA = RAW / "03_180916" / "R3_CHIME_TOA_20Feb4.txt"
OUT = ROOT / "data" / "processed" / "pointprocess"

P_TRUE = 16.35
P_TOL_FRAC = 0.01
P_LO, P_HI, N_PERIOD = 5.0, 35.0, 1200
N_PHI0, WIDTHS = 24, np.array([0.05, 0.10, 0.15, 0.20, 0.25])
N_SURR_N1 = 200
N_SYNTH = 10
N_UNIFORM = 20
N_SCRAMBLE = 200
DT_DAYS = 60.0 / 86400.0
T_CLUSTER_D = 1.0
CONTRAST = 50.0


def main() -> int:
    rng = np.random.default_rng(20261008)
    out = {"built_at": dt.datetime.now().astimezone().isoformat(timespec="seconds"),
           "criteria": "docs/TASKS.md G3 冻结口径 v1", "P_true": P_TRUE}

    mjd, meta = em.load_r3_exposure_series(R3)
    wins = em.series_to_windows(mjd)
    log = em.obslog_from_series(mjd, 29.503, 65.717)
    grid_t, grid_w = obs_grid_from_obslog(log, dt_days=DT_DAYS)
    toa_df, _ = chime.load_r3_180916_toa(TOA)
    toa = toa_df["mjd_native"].to_numpy()
    t_bjd = toa_df["t_bjd_tdb"].to_numpy()
    in_win = np.array([any(a <= x <= b for a, b in wins) for x in toa])
    events = t_bjd[in_win]
    t0 = float(events[0])
    print(f"观测网格: {grid_t.size} 点, T_tot={grid_w.sum():.1f} d; "
          f"事件集: 窗内 {events.size}/{toa.size}")

    phi0_grid = np.arange(N_PHI0) / N_PHI0
    p_grid = np.linspace(P_LO, P_HI, N_PERIOD)
    scanner = StepWindowScanner.build(grid_t, grid_w, t0, p_grid, phi0_grid, WIDTHS)

    lam, best = scanner.scan(events)
    i_hat = int(np.argmax(lam))
    p_hat = float(p_grid[i_hat])
    lam_max = float(lam[i_hat])
    _, w_hat = phi0_grid[best[i_hat][0]], WIDTHS[best[i_hat][1]]
    kui = kuiper_statistic(events, grid_t, grid_w, t0, p_hat)
    print(f"[观测] P̂={p_hat:.4f} d  Λ_max={lam_max:.2f}  φ0={phi0_grid[best[i_hat][0]]:.3f} "
          f"w={w_hat:.2f}  (Kuiper={kui:.3f})")
    out["observed"] = {"P_hat": p_hat, "lambda_max": lam_max,
                       "phi0": float(phi0_grid[best[i_hat][0]]), "width": float(w_hat),
                       "kuiper": kui, "n_events": int(events.size),
                       "dev_from_true_frac": float(abs(p_hat - P_TRUE) / P_TRUE)}

    n1_lams = ppnull.max_stat_null(scanner, grid_t, grid_w, events.size, N_SURR_N1, rng)
    p_n1 = ppnull.empirical_p(lam_max, n1_lams)
    q99 = float(np.quantile(n1_lams, 0.99))
    g3_2_ok = abs(p_hat - P_TRUE) / P_TRUE <= P_TOL_FRAC and p_n1 < 0.01
    print(f"[G3②] N1: max零分布 中位={np.median(n1_lams):.2f} q99={q99:.2f} "
          f"| p_N1={p_n1:.4f} | 偏差={abs(p_hat - P_TRUE) / P_TRUE:.3%} "
          f"→ {'OK' if g3_2_ok else 'FAIL'}")
    out["G3_2"] = {"p_n1": p_n1, "null_median": float(np.median(n1_lams)),
                   "null_q99": q99, "n_surrogate": N_SURR_N1, "ok": bool(g3_2_ok)}

    ep, ep_meta = decluster(events, t_cluster_days=T_CLUSTER_D)
    scr = np.empty(N_SCRAMBLE)
    for i in range(N_SCRAMBLE):
        shift = ppnull.scramble_episodes(ep, p_hat, rng)
        s, _ = scanner.stat_at(shift, p_hat)
        scr[i] = s
    p_n2 = ppnull.empirical_p(lam_max, scr)
    print(f"[G3②-N2] episodes={ep_meta['n_episodes']} (t_cluster={T_CLUSTER_D}d); "
          f"p_N2={p_n2:.4f} (中位={np.median(scr):.2f})")
    out["G3_2_N2"] = {"p_n2": p_n2, "n_episodes": ep_meta["n_episodes"],
                      "null_median": float(np.median(scr)),
                      "consistent_with_n1": bool((p_n2 < 0.05) == (p_n1 < 0.05))}

    hits, peaks = 0, []
    for _ in range(N_SYNTH):
        phi0_true = rng.uniform(0.0, 1.0)
        sig = ppnull.sample_periodic_signal(grid_t, grid_w, t0, P_TRUE, phi0_true,
                                            0.25, CONTRAST, 40, rng)
        lam_s, _ = scanner.scan(sig)
        pk = float(p_grid[int(np.argmax(lam_s))])
        peaks.append(pk)
        if abs(pk - P_TRUE) / P_TRUE <= P_TOL_FRAC:
            hits += 1
    g3_1a_ok = hits >= 8
    print(f"[G3①-a] 合成恢复(contrast={CONTRAST}) {hits}/{N_SYNTH} 落 ±1% → "
          f"{'OK' if g3_1a_ok else 'FAIL'} | peaks={np.round(peaks, 3).tolist()}")
    out["G3_1a"] = {"hits": hits, "n": N_SYNTH, "peaks": peaks, "contrast": CONTRAST,
                    "ok": bool(g3_1a_ok)}
    sens, sens_hits = [], 0
    for _ in range(N_SYNTH):
        phi0_true = rng.uniform(0.0, 1.0)
        sig = ppnull.sample_periodic_signal(grid_t, grid_w, t0, P_TRUE, phi0_true,
                                            0.25, 10.0, 40, rng)
        pk = float(p_grid[int(np.argmax(scanner.scan(sig)[0]))])
        sens.append(pk)
        sens_hits += abs(pk - P_TRUE) / P_TRUE <= P_TOL_FRAC
    out["G3_1a_contrast10"] = {"hits": int(sens_hits), "n": N_SYNTH, "peaks": sens}
    print(f"[G3①-a] 敏感性 contrast=10: {int(sens_hits)}/{N_SYNTH}")

    exceed = 0
    for _ in range(N_UNIFORM):
        ev = ppnull.sample_uniform_given_exposure(grid_t, grid_w, events.size, rng)
        lam_u, _ = scanner.scan(ev)
        if float(lam_u.max()) > q99:
            exceed += 1
    g3_1b_ok = exceed <= 1
    print(f"[G3①-b] 纯曝光零假报: {exceed}/{N_UNIFORM} 超 q99 → {'OK' if g3_1b_ok else 'FAIL'}")
    out["G3_1b"] = {"exceed_q99": exceed, "n": N_UNIFORM, "ok": bool(g3_1b_ok)}

    verdict = bool(g3_1a_ok and g3_1b_ok and g3_2_ok)
    out["verdict"] = "PASS" if verdict else "FAIL"
    out["engine"] = {"model": "episode-level NHPP step activity window",
                     "p_grid": [P_LO, P_HI, N_PERIOD], "n_phi0": N_PHI0,
                     "widths": WIDTHS.tolist(), "dt_days": DT_DAYS,
                     "t_cluster_days": T_CLUSTER_D}
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "g3_report.json").write_text(json.dumps(out, indent=2, ensure_ascii=False),
                                        encoding="utf-8")
    np.save(OUT / "g3_periodogram.npy", np.vstack([p_grid, lam]))
    print(f"G3: {out['verdict']}  →  data/processed/pointprocess/g3_report.json")
    return 0 if verdict else 1


if __name__ == "__main__":
    raise SystemExit(main())