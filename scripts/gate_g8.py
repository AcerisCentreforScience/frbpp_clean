"""Stage 8 acceptance gate G8.

Usage: python scripts/gate_g8.py
Output: data/processed/forecast/g8_report.json
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

from frbpp.catalogs.chime import mjd_to_bjd_tdb                    # noqa: E402
from frbpp.census.routing import build_event_frame                 # noqa: E402
from frbpp.exposure.model import obslog_from_teff                  # noqa: E402
from frbpp.forecast.predict import (make_prediction, register,     # noqa: E402
                                    verify_event, window_occurs)
from frbpp.pointprocess.alias import alias_test                    # noqa: E402
from frbpp.pointprocess.decluster import decluster                 # noqa: E402
from frbpp.pointprocess.periodogram import (StepWindowScanner,     # noqa: E402
                                            obs_grid_from_obslog)
from frbpp.upgrade.crossmatch import clean_rep                     # noqa: E402

RAW2 = ROOT / "data" / "interim" / "cat2_bursts.parquet"
SCAN = ROOT / "data" / "processed" / "census" / "a_layer_scan.parquet"
LMAP = ROOT / "data" / "raw" / "chime" / "05_cat1" / (
    "exposure_int_20180828_20190702_transit_L_beam_FWHM-600_res_4s_0.86_arcmin.npz")
UMAP = ROOT / "data" / "raw" / "chime" / "05_cat1" / (
    "exposure_int_20180828_20190702_transit_U_beam_FWHM-600_res_4s_0.86_arcmin.npz")
PRED = ROOT / "predictions"
OUT = ROOT / "data" / "processed" / "forecast"
DT = 120.0 / 86400.0
P_GRID = np.linspace(1.0, 200.0, 1500)
PHI0 = np.arange(24) / 24.0
WIDTHS = np.array([0.05, 0.10, 0.15, 0.20, 0.25])


def load_context():
    raw = pd.read_parquet(RAW2)
    evf = build_event_frame(raw)
    scan = pd.read_parquet(SCAN).set_index("src")
    rep = clean_rep(raw.repeater_name)
    src_arr = np.where(np.array([r != "" for r in rep]), rep, raw.tns_name.astype(str))
    mjd = pd.to_numeric(raw.mjd_400, errors="coerce").to_numpy(dtype=float)
    ok = np.isfinite(mjd)
    df = pd.DataFrame({"src": src_arr[ok], "mjd": mjd[ok]})
    bursts = {s: np.sort(g.mjd.to_numpy()) for s, g in df.groupby("src")}
    coords = {s: (float(g.ra_n.iloc[0]), float(g.dec_n.iloc[0]))
              for s, g in evf.groupby("src")}
    return scan, bursts, coords


def build_teff(sources):
    from astropy import units as u
    from astropy_healpix import HEALPix
    hpix = HEALPix(nside=4096, order="ring")
    t_eff = {}
    for path in (LMAP, UMAP):
        arr = np.asarray(np.load(path, allow_pickle=True)["exposure"], dtype=float)
        for s, (ra, dec) in sources.items():
            ip = int(hpix.lonlat_to_healpix(float(ra) * u.deg, float(dec) * u.deg))
            if ip < arr.size:
                t_eff[s] = t_eff.get(s, 0.0) + float(arr[ip]) * 4.0
        del arr
    return t_eff


def fine_refine(grid_t, grid_w, ep_bjd, P0, frac=0.02, n=800):
    lo, hi = P0 * (1 - frac), P0 * (1 + frac)
    sc = StepWindowScanner.build(grid_t, grid_w, float(ep_bjd[0]),
                                 np.linspace(lo, hi, n), PHI0, WIDTHS)
    lam, best = sc.scan(ep_bjd)
    j = int(np.argmax(lam))
    return {"P": float(sc.p_grid[j]), "phi0": float(PHI0[int(best[j][0])]),
            "w": float(WIDTHS[int(best[j][1])]), "sigma_P": float((hi - lo) / n)}


def occasional_hits(t_bjd, t0, P, phi0, w, sigP):
    t = np.asarray(t_bjd, dtype=float)
    tol = np.maximum(0.25, 0.5 * w * P) + 2.0 * np.abs(t - t0) * sigP / P
    tol_ph = tol / P
    ph = np.mod((t - t0) / P, 1.0)
    hit = np.mod(ph - phi0 + tol_ph, 1.0) < (w + 2 * tol_ph)
    p_eff = np.minimum(w + 2 * tol_ph, 1.0)
    return hit, p_eff


def main() -> int:
    rng = np.random.default_rng(20261014)
    out = {"built_at": dt.datetime.now().astimezone().isoformat(timespec="seconds"),
          "criteria": "docs/TASKS.md G8 冻结口径 v1(+v1.1)"}

    scan, bursts, coords = load_context()

    qual = []
    for s in scan.index:
        lab = scan.loc[s, "label"]
        if lab not in ("疑似(T1)", "疑似(T2)") or s not in bursts or s not in coords:
            continue
        P_hat = float(scan.loc[s, "P_hat"])
        span = float(bursts[s].max() - bursts[s].min())
        if span < 3 * P_hat or int(scan.loc[s, "n_episodes"]) < 8:
            continue
        qual.append({"src": s, "P_hat": P_hat, "span": span,
                     "n_ep": int(scan.loc[s, "n_episodes"]),
                     "harmonic": bool(scan.loc[s, "harmonic_ambiguous"])})
    ctrl_srcs = [s for s in scan.index
                 if scan.loc[s, "label"] == "无结论" and s in bursts and s in coords
                 and float(bursts[s].max() - bursts[s].min()) >=
                 3 * float(scan.loc[s, "P_hat"])
                 and int(scan.loc[s, "n_episodes"]) >= 8][:2]
    print(f"合格预测源: {[q['src'] for q in qual]} | 阴性对照: {ctrl_srcs}")

    need = {q["src"]: coords[q["src"]] for q in qual}
    need.update({s: coords[s] for s in ctrl_srcs})
    t_eff = build_teff(need)

    reg_rows, ok1 = [], True
    for q in qual:
        s = q["src"]
        ra, dec = coords[s]
        ep_mjd = decluster(bursts[s], 1.0)[0]
        ep_bjd = mjd_to_bjd_tdb(ep_mjd, ra, dec)
        log, _ = obslog_from_teff(ra, dec, float(ep_mjd.min()), float(ep_mjd.max()),
                                  t_eff[s])
        grid_t, grid_w = obs_grid_from_obslog(log, dt_days=DT)
        if alias_test(grid_t, grid_w, q["P_hat"])["flagged"]:
            print(f"  [跳过] {s}：别名命中")
            continue
        branches = [(q["P_hat"], "scan")]
        if q["harmonic"]:
            branches.append((2 * q["P_hat"], "fundamental"))
        for P0, tag in branches:
            ref = fine_refine(grid_t, grid_w, ep_bjd, P0, frac=0.01, n=600)
            P, phi0, w = ref["P"], ref["phi0"], ref["w"]
            stem = f"PRED_{s}_{tag}" if q["harmonic"] else f"PRED_{s}"
            pred = make_prediction(s, P, phi0, w, float(ep_bjd[0]), float(ep_bjd.max()),
                                   meta={"branch": tag, "gate": "G8",
                                         "n_episodes_train": int(ep_bjd.size),
                                         "sigma_P": ref["sigma_P"],
                                         "data": "Cat2 (via census)"},
                                   file_stem=stem)
            st, _key, sha = register(pred, PRED)
            st2, _k2, _s2 = register(pred, PRED)
            ok1 &= (st in ("registered", "exists")) and (st2 == "exists")
            reg_rows.append({"src": s, "branch": tag, "P": P, "width": w, "phi0": phi0,
                             "status": st, "idem_ok": st2 == "exists",
                             "sha256": sha[:16], "windows": pred["windows"]})
            print(f"  [注册] {stem}: P={P:.4f}±{ref['sigma_P']:.4f} w={w:.2f} "
                  f"φ0={phi0:.3f} {st}（幂等 {st2}）SHA={sha[:12]}")
    out["predictions"] = reg_rows
    out["g8_1"] = {"n_registered": len(reg_rows), "ok": bool(ok1)}
    print(f"[G8①] 注册 {len(reg_rows)} 份；幂等全过={bool(ok1)}")

    def blind(src, n_perm=500):
        t_all = bursts[src]
        cut = float(np.quantile(t_all, 0.70))
        train, hold = t_all[t_all <= cut], t_all[t_all > cut]
        if train.size < 5 or hold.size < 3:
            return None
        ra, dec = coords[src]
        ep_mjd = decluster(train, 1.0)[0]
        ep_tr = mjd_to_bjd_tdb(ep_mjd, ra, dec)
        log, _ = obslog_from_teff(ra, dec, float(ep_mjd.min()), float(ep_mjd.max()),
                                  t_eff[src])
        grid_t, grid_w = obs_grid_from_obslog(log, dt_days=DT)
        sc = StepWindowScanner.build(grid_t, grid_w, float(ep_tr[0]),
                                     P_GRID, PHI0, WIDTHS)
        lam, _ = sc.scan(ep_tr)
        P_coarse = float(P_GRID[int(np.argmax(lam))])
        hold_bjd = mjd_to_bjd_tdb(hold, ra, dec)
        res = {}
        for fac, tag in [(1.0, "scan"), (2.0, "fundamental")]:
            ref = fine_refine(grid_t, grid_w, ep_tr, P_coarse * fac, frac=0.02, n=800)
            hit, p_eff = occasional_hits(hold_bjd, float(ep_tr[0]), ref["P"],
                                         ref["phi0"], ref["w"], ref["sigma_P"])
            n_hit = int(hit.sum())
            E = float(p_eff.sum())
            V = float((p_eff * (1 - p_eff)).sum())
            counts = np.empty(n_perm, dtype=int)
            for i in range(n_perm):
                shift = rng.uniform(0.0, ref["P"])
                h2, _ = occasional_hits(hold_bjd + shift, float(ep_tr[0]), ref["P"],
                                        ref["phi0"], ref["w"], ref["sigma_P"])
                counts[i] = int(h2.sum())
            p_mc = (1 + int(np.sum(counts >= n_hit))) / (n_perm + 1)
            z_mc = (n_hit - float(counts.mean())) / max(float(counts.std()), 1e-9)
            res[tag] = {"P_train": ref["P"], "P_coarse": P_coarse, "w": ref["w"],
                        "phi0": ref["phi0"], "sigma_P": ref["sigma_P"],
                        "n_hold": int(hold.size), "hits": n_hit,
                        "rate": n_hit / hold.size, "E": E, "sqrtV": float(np.sqrt(V)),
                        "p_mc": p_mc, "z_mc": float(z_mc)}
        return res

    blind_rows = {q["src"]: blind(q["src"]) for q in qual}
    blind_rows = {k: v for k, v in blind_rows.items() if v}
    pc = blind_rows.get("FRB20180916B", {})
    pc_ok = False
    for tag, v in pc.items():
        crit = bool(v["rate"] >= 0.8)
        v["criterion_ok"] = crit
        pc_ok = pc_ok or crit
        print(f"  [盲搜] FRB20180916B {tag:11s} P_tr={v['P_train']:.3f}±{v['sigma_P']:.4f} "
              f"回收 {v['hits']}/{v['n_hold']} = {v['rate']:.2f}（E={v['E']:.1f}±{v['sqrtV']:.1f}, "
              f"p_mc={v['p_mc']:.4f}≈覆盖率下限, z_mc={v['z_mc']:+.2f}）→ {'OK' if crit else 'FAIL'}")
    ctrl_rows, ctrl_ok = {}, True
    for s in ctrl_srcs:
        r = blind(s)
        if not r:
            continue
        v = r["scan"]
        good = bool(abs(v["z_mc"]) <= 3.0)
        ctrl_ok &= good
        ctrl_rows[s] = {**v, "ok": good}
        print(f"  [阴性对照] {s}: 命中 {v['hits']}/{v['n_hold']}（E={v['E']:.1f}, "
              f"z_mc={v['z_mc']:+.2f}）→ {'OK' if good else 'FAIL'}")
    ok2 = bool(pc_ok and ctrl_ok)
    out["blind"] = blind_rows
    out["negative_controls"] = ctrl_rows
    out["g8_2"] = {"pc_ok": bool(pc_ok), "ctrl_ok": bool(ctrl_ok), "ok": ok2}
    print(f"[G8②] 180916 主判据 {'OK' if pc_ok else 'FAIL'} | "
          f"阴性对照 {'OK' if ctrl_ok else 'FAIL'}")

    dry_hist = PRED / "HISTORY_DRYRUN.jsonl"
    anchor = next((r for r in reg_rows if r["src"] == "FRB20180916B"), reg_rows[0])
    w0 = anchor["windows"][0]
    t_hit = 0.5 * (w0[0] + w0[1])
    r1 = verify_event(anchor["src"], t_hit, PRED, dry_hist)
    t_miss, r2 = None, None
    for off in range(3, 40):
        cand = w0[0] - float(off)
        st = "MISS"
        for fp in PRED.glob("PRED_FRB20180916B*.json"):
            p = json.loads(fp.read_text(encoding="utf-8"))
            for a, b in p["windows"]:
                if a - p["tolerance_days"] <= cand <= b + p["tolerance_days"]:
                    st = "HIT"
        if st == "MISS":
            t_miss = cand
            break
    if t_miss is not None:
        r2 = verify_event(anchor["src"], t_miss, PRED, dry_hist)
    r3 = verify_event("FRB000000X", t_hit, PRED, dry_hist)
    r4 = verify_event(anchor["src"], t_hit, PRED, dry_hist)
    ok3 = bool(r1["status"] == "HIT" and r2 is not None and r2["status"] == "MISS"
               and r3["status"] == "NO_MATCH" and r4["duplicate"] is True)
    print(f"[G8③ dry-run] HIT={r1['status']} MISS={r2['status']} "
          f"NO_MATCH={r3['status']} 幂等={r4['duplicate']} → {'OK' if ok3 else 'FAIL'}")
    out["g8_3"] = {"hit": r1["status"], "miss": r2["status"], "no_match": r3["status"],
                   "duplicate": bool(r4["duplicate"]), "ok": ok3}

    verdict = bool(ok1 and ok2 and ok3 and len(reg_rows) > 0)
    out["verdict"] = "PASS" if verdict else "FAIL"
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "g8_report.json").write_text(json.dumps(out, indent=2, ensure_ascii=False,
                                                   default=str), encoding="utf-8")
    print(f"G8: {out['verdict']}  →  data/processed/forecast/g8_report.json")
    return 0 if verdict else 1


if __name__ == "__main__":
    raise SystemExit(main())