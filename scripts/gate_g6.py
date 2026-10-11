"""Stage 6 acceptance gate G6.

Usage: python scripts/gate_g6.py
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

from frbpp.census.routing import build_event_frame, route_sources      # noqa: E402
from frbpp.exposure.model import obslog_from_teff                      # noqa: E402
from frbpp.models.ladder import ladder_fit                             # noqa: E402
from frbpp.pointprocess.alias import alias_test                        # noqa: E402
from frbpp.pointprocess.null import sample_periodic_signal             # noqa: E402
from frbpp.pointprocess.periodogram import obs_grid_from_obslog        # noqa: E402

RAW2 = ROOT / "data" / "interim" / "cat2_bursts.parquet"
SCAN = ROOT / "data" / "processed" / "census" / "a_layer_scan.parquet"
LMAP = ROOT / "data" / "raw" / "chime" / "05_cat1" / (
    "exposure_int_20180828_20190702_transit_L_beam_FWHM-600_res_4s_0.86_arcmin.npz")
UMAP = ROOT / "data" / "raw" / "chime" / "05_cat1" / (
    "exposure_int_20180828_20190702_transit_U_beam_FWHM-600_res_4s_0.86_arcmin.npz")
OUT = ROOT / "data" / "processed" / "models"

KNOWN_299 = ["FRB20220912A", "FRB20201124A", "FRB20190303A", "FRB20190117A",
             "FRB20200223B", "FRB20220529A", "FRB20190208A"]
P_SYN = 16.35
N_REP = 10
DT_DAYS = 120.0 / 86400.0


def _sample_weighted(grid_t, grid_w, lam, n, rng):
    p = np.asarray(grid_w) * np.asarray(lam)
    p = np.clip(p, 0, None)
    p = p / p.sum()
    idx = rng.choice(len(grid_t), size=int(n), replace=True, p=p)
    return np.sort(np.asarray(grid_t)[idx])


def sample_sine(grid_t, grid_w, t0, P, phi, A, n, rng):
    th = 2 * np.pi * np.mod((np.asarray(grid_t) - t0) / P, 1.0)
    return _sample_weighted(grid_t, grid_w, 1 + A * np.sin(th + phi), n, rng)


def sample_trend_window(grid_t, grid_w, t0, P, phi0, w, b, n, rng):
    g = np.asarray(grid_t)
    tau = 2.0 * (g - (g.min() + g.max()) / 2.0) / max(g.max() - g.min(), 1e-9)
    Wm = np.mod(np.mod((g - t0) / P, 1.0) - phi0, 1.0) < w
    lam = Wm * (1.0 + b * tau) + 0.05
    return _sample_weighted(grid_t, grid_w, lam, n, rng)


def main() -> int:
    rng = np.random.default_rng(20261011)
    out = {"built_at": dt.datetime.now().astimezone().isoformat(timespec="seconds"),
           "criteria": "docs/TASKS.md G6 冻结口径 v1"}

    raw = pd.read_parquet(RAW2)
    events = build_event_frame(raw)
    sources, episodes = route_sources(events)
    scan = pd.read_parquet(SCAN).set_index("src")
    a_src = sources[sources.category == "A"].copy()

    from astropy import units as u
    from astropy_healpix import HEALPix
    hpix = HEALPix(nside=4096, order="ring")
    t_eff = {}
    for _, path in (("L", LMAP), ("U", UMAP)):
        arr = np.asarray(np.load(path, allow_pickle=True)["exposure"], dtype=float)
        for r in a_src.itertuples():
            ipix = int(hpix.lonlat_to_healpix(float(r.ra) * u.deg, float(r.dec) * u.deg))
            if ipix < arr.size:
                t_eff[r.src] = t_eff.get(r.src, 0.0) + float(arr[ipix]) * 4.0
        del arr

    grids, starts = {}, {}
    for r in a_src.itertuples():
        ep = episodes[r.src]
        log, _ = obslog_from_teff(float(r.ra), float(r.dec), float(ep.min()),
                                  float(ep.max()), t_eff[r.src])
        grids[r.src] = obs_grid_from_obslog(log, dt_days=DT_DAYS)
        starts[r.src] = log.t_start

    alias_rows = {}
    for r in a_src.itertuples():
        p_hat = float(scan.loc[r.src, "P_hat"])
        alias_rows[r.src] = {**alias_test(grids[r.src][0], grids[r.src][1], p_hat),
                             "P_hat": p_hat}
    hit_map = {s: bool(alias_rows[s]["flagged"]) for s in KNOWN_299}
    hit299 = bool(alias_rows["FRB20190303A"]["flagged"])
    pc_alias = alias_rows["FRB20180916B"]["flagged"]
    bstarts = np.arange(0.0, 800.0, 0.99727)
    beat = alias_test(bstarts, np.ones_like(bstarts), p_hat=2.9918)
    ok4 = bool(hit299 and (not pc_alias) and beat["flagged"])
    print(f"[G6④ 别名] 20190303A 命中={hit299} | 180916 flagged={pc_alias} | "
          f"合成节拍 rank={beat['rank']:.4f} → {'OK' if ok4 else 'FAIL'}")
    print(f"   （原 2.99 簇逐源报告: { {k: v for k, v in hit_map.items()} }）")
    out["alias"] = {"rows": alias_rows, "hit_20190303A": hit299,
                    "known299_map": hit_map,
                    "pc_180916_flagged": bool(pc_alias), "synthetic_beat": beat, "ok": ok4}

    gt, gw = grids["FRB20180916B"]
    t0g = float(gt.min())

    def run_batch(sampler, n_runs=N_REP):
        finals, deg_flags = [], []
        for _ in range(n_runs):
            ev = sampler()
            res = ladder_fit(ev, float(ev[0]), P_SYN, gt, gw)
            finals.append(res["final"])
            deg_flags.append(bool(res["degenerate"]))
        return finals, deg_flags

    f1, _ = run_batch(lambda: sample_periodic_signal(gt, gw, float(gt.min()), P_SYN,
                                                     float(rng.uniform()), 0.15, 50.0, 40, rng))
    ok1 = sum(1 for x in f1 if x == "M1") >= 8
    print(f"[G6① 方波→M1] {f1} → {'OK' if ok1 else 'FAIL'}")

    f2, _ = run_batch(lambda: sample_trend_window(gt, gw, float(gt.min()), P_SYN,
                                                  float(rng.uniform()), 0.2, 0.8, 150, rng))
    ok2 = sum(1 for x in f2 if x == "M3") >= 8
    print(f"[G6② 趋势+周期→M3] {f2} → {'OK' if ok2 else 'FAIL'}")

    f3, dflags = run_batch(
        lambda: sample_periodic_signal(gt, gw, float(gt.min()), P_SYN,
                                       float(rng.uniform()), 0.45, 50.0, 40, rng),
        n_runs=20)
    correct3 = sum(1 for f_, d_ in zip(f3, dflags) if d_ or f_ == "M1")
    overclaim3 = sum(1 for f_, d_ in zip(f3, dflags) if (f_ == "M2" and not d_))
    ok3 = bool(correct3 >= 18 and overclaim3 == 0)
    print(f"[G6③ 简并/真形] correct {correct3}/20（简并标记或 M1）| 高置信错判 {overclaim3} "
          f"→ {'OK' if ok3 else 'FAIL'}")
    out["synthetic"] = {"step_to_M1": f1, "ok1": bool(ok1),
                        "trend_to_M3": f2, "ok2": bool(ok2),
                        "degenerate_correct": int(correct3), "overclaim": int(overclaim3),
                        "ok3": bool(ok3)}

    t2 = [s for s in scan.index if scan.loc[s, "label"] == "疑似(T2)"]
    from frbpp.catalogs.chime import mjd_to_bjd_tdb as _m2b
    refit = {}
    for s in t2:
        rrow = a_src[a_src.src == s].iloc[0]
        ep = _m2b(episodes[s], float(rrow.ra), float(rrow.dec))
        res = ladder_fit(ep, float(ep[0]), float(scan.loc[s, "P_hat"]), *grids[s])
        refit[s] = {"final": res["final"], "supported": res["supported"],
                    "degenerate": res["degenerate"],
                    "M1_width": res["table"]["M1"]["det"].get("width"),
                    "alias": alias_rows[s]["flagged"]}
        print(f"  [复核] {s:14s} P̂={scan.loc[s, 'P_hat']:8.2f} → {res['final']} "
              f"(支持集 {','.join(res['supported']) or '无'})"
              f"{' 简并' if res['degenerate'] else ''}"
              f"{' 别名' if alias_rows[s]['flagged'] else ''} w*={refit[s]['M1_width']}")
    out["refit"] = refit

    checks = {"step_to_M1": bool(ok1), "trend_to_M3": bool(ok2),
              "degenerate": bool(ok3), "alias": bool(ok4)}
    verdict = all(checks.values())
    out["checks"] = checks
    out["verdict"] = "PASS" if verdict else "FAIL"
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "g6_report.json").write_text(json.dumps(out, indent=2, ensure_ascii=False,
                                                   default=str), encoding="utf-8")
    print(f"G6: {out['verdict']}  →  data/processed/models/g6_report.json")
    return 0 if verdict else 1


if __name__ == "__main__":
    raise SystemExit(main())