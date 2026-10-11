"""Stage 5 acceptance gate G5: all-source census.

Usage: python scripts/gate_g5.py [--workers N]
"""
from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import json
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from frbpp.census.routing import build_event_frame, route_sources     # noqa: E402
from frbpp.exposure.model import obslog_from_teff                     # noqa: E402
from frbpp.pointprocess.fdr import bh_fdr                              # noqa: E402

RAW2 = ROOT / "data" / "interim" / "cat2_bursts.parquet"
LMAP = ROOT / "data" / "raw" / "chime" / "05_cat1" / (
    "exposure_int_20180828_20190702_transit_L_beam_FWHM-600_res_4s_0.86_arcmin.npz")
UMAP = ROOT / "data" / "raw" / "chime" / "05_cat1" / (
    "exposure_int_20180828_20190702_transit_U_beam_FWHM-600_res_4s_0.86_arcmin.npz")
OUT = ROOT / "data" / "processed" / "census"

P_LO, P_HI, N_P = 1.0, 200.0, 1500
PHI0 = np.arange(24) / 24.0
WIDTHS = np.array([0.05, 0.10, 0.15, 0.20, 0.25])
DT_DAYS = 120.0 / 86400.0
N_SURR = 300
MASTER_SEED = "frbpp-census-v1-20261008"
P_TRUE_180916 = 16.35


def _scan_source(task):
    src, ra, dec, ep, t_start, t_end, sens = task
    from frbpp.contracts import ObservationLog
    from frbpp.pointprocess.null import empirical_p, max_stat_null
    from frbpp.pointprocess.periodogram import (StepWindowScanner, kuiper_statistic,
                                                obs_grid_from_obslog)

    ep = np.asarray(ep, dtype=float)
    log = ObservationLog(t_start=np.asarray(t_start), t_end=np.asarray(t_end),
                         sens=np.asarray(sens), tgt_ra=float(ra), tgt_dec=float(dec),
                         beam_ra=np.full(len(t_start), np.nan),
                         beam_dec=np.full(len(t_start), np.nan))
    grid_t, grid_w = obs_grid_from_obslog(log, dt_days=DT_DAYS)
    scanner = StepWindowScanner.build(grid_t, grid_w, float(ep[0]),
                                      np.linspace(P_LO, P_HI, N_P), PHI0, WIDTHS)
    lam, best = scanner.scan(ep)
    ihat = int(np.argmax(lam))
    p_hat = float(scanner.p_grid[ihat])
    lam_max = float(lam[ihat])
    i_phi, i_w = int(best[ihat][0]), int(best[ihat][1])
    w_hat = float(WIDTHS[i_w])
    ph = np.mod((ep - float(ep[0])) / p_hat, 1.0)
    occ = float((np.mod(ph - PHI0[i_phi], 1.0) < w_hat).mean())

    seed = int(hashlib.sha256((src + "|" + MASTER_SEED).encode()).hexdigest()[:16], 16)
    rng = np.random.default_rng(seed)
    nulls = max_stat_null(scanner, grid_t, grid_w, ep.size, N_SURR, rng)
    p_global = empirical_p(lam_max, nulls)

    k = ep.size // 2
    lam1 = scanner.scan(ep[:k])[0]
    lam2 = scanner.scan(ep[k:])[0]
    p1 = float(scanner.p_grid[int(np.argmax(lam1))])
    p2 = float(scanner.p_grid[int(np.argmax(lam2))])
    step = float(scanner.p_grid[1] - scanner.p_grid[0])
    p_at_bound = bool(p_hat <= P_LO + 1.5 * step or p_hat >= P_HI - 1.5 * step)
    width_at_bound = bool(w_hat >= float(WIDTHS.max()) - 1e-9)
    seg_ok = bool(abs(p1 - p_hat) <= 0.01 * p_hat and abs(p2 - p_hat) <= 0.01 * p_hat)
    kui = kuiper_statistic(ep, grid_t, grid_w, float(ep[0]), p_hat)
    guards_clean = bool((not p_at_bound) and (not width_at_bound) and seg_ok and occ >= 0.5)
    j2 = int(np.argmin(np.abs(scanner.p_grid - 2.0 * p_hat)))
    jh = int(np.argmin(np.abs(scanner.p_grid - 0.5 * p_hat)))
    harmonic_dlam = float(min(abs(lam_max - lam[j2]), abs(lam_max - lam[jh])))
    harmonic_ambiguous = bool(harmonic_dlam < 3.0)
    return {"src": src, "P_hat": p_hat, "lambda_max": lam_max, "w_hat": w_hat,
            "occupancy": occ, "p_global": p_global,
            "null_median": float(np.median(nulls)), "null_q99": float(np.quantile(nulls, 0.99)),
            "half_P1": p1, "half_P2": p2, "seg_reproducible": seg_ok,
            "p_at_bound": p_at_bound, "width_at_bound": width_at_bound,
            "guards_clean": guards_clean, "kuiper": float(kui),
            "harmonic_dlam": harmonic_dlam, "harmonic_ambiguous": harmonic_ambiguous,
            "n_grid": int(grid_t.size), "n_episodes": int(ep.size)}


def _label(q, guards_clean):
    if q < 0.05 and guards_clean:
        return "疑似(T1)"
    if q < 0.10 or (q < 0.05 and not guards_clean):
        return "疑似(T2)"
    return "无结论"


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=0)
    args = ap.parse_args()

    raw = pd.read_parquet(RAW2)
    events = build_event_frame(raw)
    sources, episodes = route_sources(events)
    counts = sources.category.value_counts().to_dict()
    print(f"事件 {len(events)} 行 | 源 {len(sources)} | 路由 {counts}")

    a_src = sources[sources.category == "A"].sort_values("n_episodes", ascending=False)
    print("A 层源:", ", ".join(a_src.src.head(12)), f"…（共 {len(a_src)}）")

    from astropy import units as u
    from astropy_healpix import HEALPix
    hpix = HEALPix(nside=4096, order="ring")
    from frbpp.upgrade.crossmatch import clean_rep
    rep_names = set(r for r in clean_rep(raw.repeater_name) if r)
    a_src = a_src.assign(is_rep=a_src.src.isin(rep_names))
    t_eff = {}
    for name, path in (("L", LMAP), ("U", UMAP)):
        arr = np.asarray(np.load(path, allow_pickle=True)["exposure"], dtype=float)
        for r in a_src.itertuples():
            ipix = int(hpix.lonlat_to_healpix(float(r.ra) * u.deg, float(r.dec) * u.deg))
            if ipix < arr.size:
                t_eff[r.src] = t_eff.get(r.src, 0.0) + float(arr[ipix]) * 4.0
        del arr
    print("T_eff 提取完成（L+U）")

    from frbpp.catalogs.chime import mjd_to_bjd_tdb   # noqa: E402
    tasks = []
    for r in a_src.itertuples():
        ep = episodes[r.src]
        ep_bjd = mjd_to_bjd_tdb(ep, float(r.ra), float(r.dec))
        log, _meta = obslog_from_teff(float(r.ra), float(r.dec),
                                      float(ep.min()), float(ep.max()),
                                      t_eff[r.src])
        tasks.append((r.src, float(r.ra), float(r.dec), ep_bjd,
                      log.t_start, log.t_end, log.sens))

    n_workers = args.workers if args.workers > 0 else max(1, min(6, (os.cpu_count() or 4) - 1))
    print(f"启动 {n_workers} 进程 × {len(tasks)} 源 …")
    import multiprocessing as mp
    ctx = mp.get_context("spawn")
    with ctx.Pool(n_workers) as pool:
        results = pool.map(_scan_source, tasks)
    res = pd.DataFrame(results)
    q = bh_fdr(res.p_global.to_numpy())
    res["q_value"] = q
    res["label"] = [_label(qi, gi) for qi, gi in zip(res.q_value, res.guards_clean)]
    OUT.mkdir(parents=True, exist_ok=True)
    res.to_parquet(OUT / "a_layer_scan.parquet", index=False)
    print("\nA 层结果（按 q 排序）:")
    show = res.sort_values(["q_value", "lambda_max"], ascending=[True, False])
    for r in show.head(12).itertuples():
        print(f"  {r.src:16s} P̂={r.P_hat:8.2f} d Λ={r.lambda_max:7.1f} p={r.p_global:.4f} "
              f"q={r.q_value:.4f} w={r.w_hat:.2f} occ={r.occupancy:.2f} {r.label}")

    census = sources.merge(res.drop(columns=["n_episodes"]), on="src", how="left")
    census["label"] = census["label"].fillna("不适用(路由)")
    census["sigma_inflation"] = np.nan

    a = census[census.category == "A"]
    checks = {}
    checks["routing_no_violation"] = bool(
        (census.loc[census.category != "A", "label"] == "不适用(路由)").all())
    checks["no_bare_p"] = bool(
        a[["p_global", "q_value", "guards_clean", "n_episodes", "label"]].notna().all().all())
    claims = a[a.label != "无结论"]
    checks["suspected_only"] = bool(
        all(str(x).startswith("疑似") for x in claims.label))
    pc = census[census.src == "FRB20180916B"].iloc[0]
    d_1 = abs(pc.P_hat - P_TRUE_180916) / P_TRUE_180916
    d_2 = abs(2 * pc.P_hat - P_TRUE_180916) / P_TRUE_180916
    d_h = abs(pc.P_hat / 2 - P_TRUE_180916) / P_TRUE_180916
    p_floor = 1.0 / (1.0 + N_SURR)
    pc_ok = bool(pc.p_global <= p_floor + 1e-12 and min(d_1, d_2, d_h) <= 0.01)
    checks["positive_control_180916"] = pc_ok
    print(f"\n[G5] 180916: P̂={pc.P_hat:.3f} d p={pc.p_global:.4f}(floor={p_floor:.4f}) "
          f"q={pc.q_value:.4f}(报告) label={pc.label} 谐波歧义={pc.harmonic_ambiguous}"
          f"(ΔΛ={pc.harmonic_dlam:.2f}) 偏差(1x/2x/half)={d_1:.3%}/{d_2:.3%}/{d_h:.3%} "
          f"→ {'OK' if pc_ok else 'FAIL'}")
    print(f"[G5] 标签分布（A 层）: {a.label.value_counts().to_dict()}")

    verdict = all(checks.values())
    out = {"built_at": dt.datetime.now().astimezone().isoformat(timespec="seconds"),
           "criteria": "docs/TASKS.md G5 冻结口径 v1",
           "routing_counts": {k: int(v) for k, v in counts.items()},
           "checks": checks, "verdict": "PASS" if verdict else "FAIL",
           "a_layer": res.to_dict(orient="records"),
           "config": {"P_range": [P_LO, P_HI], "n_p": N_P, "n_phi0": len(PHI0),
                      "widths": WIDTHS.tolist(), "dt_days": DT_DAYS,
                      "n_surrogate": N_SURR, "master_seed": MASTER_SEED}}
    OUT.mkdir(parents=True, exist_ok=True)
    census.to_parquet(OUT / "census_v1.parquet", index=False)
    (OUT / "g5_report.json").write_text(json.dumps(out, indent=2, ensure_ascii=False),
                                        encoding="utf-8")
    print(f"G5: {out['verdict']}  →  data/processed/census/{{census_v1.parquet, g5_report.json}}")
    return 0 if verdict else 1


if __name__ == "__main__":
    raise SystemExit(main())