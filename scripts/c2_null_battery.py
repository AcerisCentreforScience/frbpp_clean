"""Null battery for candidate significance: permutation models, multiplicity, injection recovery, physical consistency.

Usage: python scripts/c2_null_battery.py [--sections ...] [--B 1000] [--seed 20261009]
Output: data/processed/upgrade/c2_null_battery.json
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
sys.path.insert(0, str(ROOT / "src"))

from frbpp.upgrade.crossmatch import (                       # noqa: E402
    LOGIC_VERSION, EventTable, build_event_table, candidate_groups, clean_rep,
    cone_pairs, dm_shuffle_counts, geometry_pairs, repeater_recall, same_event_mask)

RAW1 = ROOT / "data" / "interim" / "cat1_bursts.parquet"
RAW2 = ROOT / "data" / "interim" / "cat2_bursts.parquet"
RAW2CSV = ROOT / "data" / "raw" / "chime" / "cat2_25" / "chimefrbcat2.csv"
BVJSON = ROOT / "data" / "raw" / "blinkverse" / "blinkverse_sources.json"
CAND1 = ROOT / "data" / "processed" / "upgrade" / "c2_candidates_v1.parquet"
OUT = ROOT / "data" / "processed" / "upgrade"

THETAS = (60.0, 180.0, 600.0)
DMWS = (1.0, 5.0, 20.0)
N_SHUFFLE = 200
DEC_BINS = 8
N_EFF_TRIALS = 9


def haversine_arcsec(ra1, dec1, ra2, dec2):
    r1, d1, r2, d2 = map(np.radians, [ra1, dec1, ra2, dec2])
    x = np.clip(np.sin(d1) * np.sin(d2) + np.cos(d1) * np.cos(d2) * np.cos(r1 - r2), -1, 1)
    return float(np.degrees(np.arccos(x)) * 3600.0)


def build_with_aux(c1: pd.DataFrame, c2: pd.DataFrame, event_dedup: bool = True):
    ev = build_event_table(c1, c2, event_dedup=event_dedup)
    num = lambda s: pd.to_numeric(s, errors="coerce").to_numpy(dtype=float)
    ra = np.concatenate([num(c1.ra), num(c2.ra)])
    dec = np.concatenate([num(c1.dec), num(c2.dec)])
    ok = np.isfinite(ra) & np.isfinite(dec)
    mjd = np.concatenate([num(c1.mjd_400), num(c2.mjd_400)])
    ok = ok & np.isfinite(mjd)
    if event_dedup:
        from frbpp.upgrade.crossmatch import _event_layer_keep, _str
        name = np.concatenate([_str(c1.tns_name), _str(c2.tns_name)])
        cat = np.concatenate([np.full(len(c1), "cat1"), np.full(len(c2), "cat2")])
        ok = ok & _event_layer_keep(name, cat, mjd, ok)
    dmerr = np.concatenate([num(c1.dm_fitb_err), num(c2.dm_fitb_err)])[ok]
    snr = np.concatenate([num(c1.bonsai_snr), num(c2.bonsai_snr)])[ok]
    return ev, dmerr, snr


def subset(ev: EventTable, m: np.ndarray) -> EventTable:
    return EventTable(*(a[m] for a in (ev.ra, ev.dec, ev.mjd, ev.dm, ev.name,
                                       ev.catalog, ev.src_rep, ev.multi_component, ev.excluded)))


def grid_search(ev: EventTable, rng, n_shuffle=N_SHUFFLE, thetas=THETAS, dmws=DMWS,
                dedup_bg: bool = False):
    cfgs = []
    for th in thetas:
        gi, gj = geometry_pairs(ev, th)
        keep = ~same_event_mask(ev, gi, gj)
        for dw in dmws:
            ddm = np.abs(ev.dm[gi] - ev.dm[gj])
            ok = keep & np.isfinite(ddm) & (ddm < dw)
            obs = int(ok.sum())
            bgi, bgj = (gi[keep], gj[keep]) if dedup_bg else (gi, gj)
            bg = dm_shuffle_counts(ev, bgi, bgj, dw, n_shuffle, rng)
            mean = float(bg.mean())
            thr = mean + 3.0 * np.sqrt(mean)
            i, j = gi[ok], gj[ok]
            cfgs.append({"theta": th, "dm_window": dw, "obs": obs, "bg_mean": mean,
                         "thr": thr, "sig": bool(obs > thr), "i": i, "j": j,
                         "n_geom": int(gi.size), "n_dedup": int(keep.sum()),
                         "ngroups": len(candidate_groups(ev, i, j))})
    return cfgs


def union_candidates(cfgs, ev: EventTable):
    best = {}
    for c in cfgs:
        if not c["sig"]:
            continue
        for a, b in zip(c["i"].tolist(), c["j"].tolist()):
            key = tuple(sorted([str(ev.name[a]), str(ev.name[b])]))
            if key not in best or c["theta"] < best[key]["theta"]:
                best[key] = {"theta": c["theta"], "a": a, "b": b}
    return best


def cand_stats(best, ev, dmerr):
    n = len(best)
    if n == 0:
        return {"n_cand": 0, "max_ddm_sigma": 0.0, "min_theta": float("nan"),
                "n_sig_cfg": 0, "sum_groups": 0}
    ms, ths = [], []
    for v in best.values():
        a, b = v["a"], v["b"]
        e = float(np.hypot(dmerr[a], dmerr[b]))
        ms.append(abs(ev.dm[a] - ev.dm[b]) / e if e > 0 else np.inf)
        ths.append(haversine_arcsec(ev.ra[a], ev.dec[a], ev.ra[b], ev.dec[b]))
    return {"n_cand": n, "max_ddm_sigma": float(np.max(ms)), "min_theta": float(np.min(ths))}


def one_rep(evp, dmerr_p, rng, n_shuffle=N_SHUFFLE):
    cfgs = grid_search(evp, rng, n_shuffle)
    best = union_candidates(cfgs, evp)
    st = cand_stats(best, evp, dmerr_p)
    st["n_sig_cfg"] = int(sum(c["sig"] for c in cfgs))
    st["sum_groups"] = int(sum(c["ngroups"] for c in cfgs))
    return st, best, cfgs


def perm_ra_strat(ev, rng, nbins=DEC_BINS):
    ra = ev.ra.copy()
    order = np.argsort(ev.dec)
    for b in np.array_split(order, nbins):
        ra[b] = rng.permutation(ev.ra[b])
    return EventTable(ra, ev.dec, ev.mjd, ev.dm, ev.name, ev.catalog, ev.src_rep,
                      ev.multi_component, ev.excluded)


def perm_dm_shuffle(ev, dmerr, rng):
    idx = rng.permutation(ev.dm.size)
    return (EventTable(ev.ra, ev.dec, ev.mjd, ev.dm[idx], ev.name, ev.catalog,
                       ev.src_rep, ev.multi_component, ev.excluded), dmerr[idx])


def perm_ra_dm(ev, dmerr, rng):
    evp = perm_ra_strat(ev, rng)
    idx = rng.permutation(ev.dm.size)
    return (EventTable(evp.ra, evp.dec, evp.mjd, ev.dm[idx], ev.name, ev.catalog,
                       ev.src_rep, evp.multi_component, ev.excluded), dmerr[idx])


def perm_full(ev, dmerr, rng):
    n = ev.dm.size
    ra = rng.uniform(0.0, 360.0, size=n)
    dec = ev.dec[rng.permutation(n)]
    dm = ev.dm[rng.permutation(n)]
    mjd = ev.mjd[rng.permutation(n)]
    return (EventTable(ra, dec, mjd, dm, ev.name, ev.catalog, ev.src_rep,
                       ev.multi_component, ev.excluded), dmerr[rng.permutation(n)])


def run_permutation(kind, ev_oo, dmerr_oo, real, B, seed):
    rng = np.random.default_rng(seed)
    rec = {k: np.empty(B) for k in
           ("n_cand", "n_sig_cfg", "sum_groups", "max_ddm_sigma", "min_theta")}
    for b in range(B):
        if kind == "ra":
            evp, dp = perm_ra_strat(ev_oo, rng), dmerr_oo
        elif kind == "dm":
            evp, dp = perm_dm_shuffle(ev_oo, dmerr_oo, rng)
        elif kind == "ra_dm":
            evp, dp = perm_ra_dm(ev_oo, dmerr_oo, rng)
        else:
            evp, dp = perm_full(ev_oo, dmerr_oo, rng)
        st, _, _ = one_rep(evp, dp, rng)
        for k in rec:
            rec[k][b] = st[k]
    n = rec["n_cand"]
    real_n = real["n_cand"]
    ge = int(np.sum(n >= real_n))
    single_p = (1 + ge) / (B + 1)
    return {
        "n_reps": B,
        "cand_mean": float(n.mean()), "cand_p95": float(np.quantile(n, 0.95)),
        "cand_max": float(n.max()), "cand_min": float(n.min()),
        "cand_median": float(np.median(n)),
        "real_cand": int(real_n),
        "real_percentile": float(np.mean(n < real_n)),
        "n_ge_real": ge, "single_p": float(single_p),
        "global_p_bonferroni": float(min(1.0, single_p * N_EFF_TRIALS)),
        "n_sig_cfg_mean": float(rec["n_sig_cfg"].mean()),
        "sum_groups_mean": float(rec["sum_groups"].mean()),
        "max_ddm_sigma_mean": float(rec["max_ddm_sigma"].mean()),
        "max_ddm_sigma_p95": float(np.quantile(rec["max_ddm_sigma"], 0.95)),
        "min_theta_mean": (float(np.nanmean(rec["min_theta"]))
                           if np.isfinite(rec["min_theta"]).any() else float("nan")),
        "hist": np.histogram(n, bins=np.arange(0, max(30, int(n.max()) + 2)))[0].tolist(),
        "hist_edges": np.arange(0, max(30, int(n.max()) + 2)).tolist(),
    }


def known_repeater_analysis(ev, dmerr):
    srcs = sorted({r for r in ev.src_rep if r})
    rows, miss = [], {}
    ddm_sigma_all, theta_all = [], []
    for s in srcs:
        idx = np.where(ev.src_rep == s)[0]
        n = idx.size
        pairs = []
        for a in range(n):
            for b in range(a + 1, n):
                ia, ib = idx[a], idx[b]
                if same_event_mask(ev, np.array([ia]), np.array([ib]))[0]:
                    continue
                th = haversine_arcsec(ev.ra[ia], ev.dec[ia], ev.ra[ib], ev.dec[ib])
                ddm = abs(ev.dm[ia] - ev.dm[ib])
                e = float(np.hypot(dmerr[ia], dmerr[ib]))
                pairs.append((th, ddm, ddm / e if e > 0 else np.inf))
        if pairs:
            theta_all += [p[0] for p in pairs]
            ddm_sigma_all += [p[2] for p in pairs if p[0] < 600.0]
        within = [p for p in pairs if p[0] < 600.0]
        rec = any(p[1] < 5.0 for p in within)
        if rec:
            reason = "已召回"
        elif n < 2:
            reason = "目录事件<2（无法配对）"
        elif not within:
            reason = "角距>θ（逐暴定位散布）"
        elif not any(p[1] < 5.0 for p in within):
            reason = "DM变化>容差(±5)"
        else:
            reason = "其他"
        miss[reason] = miss.get(reason, 0) + 1
        bp = min(pairs, key=lambda p: p[0]) if pairs else None
        rows.append({"src": s, "n_events": int(n), "recovered": bool(rec),
                     "min_theta_arcsec": float(bp[0]) if bp else float("nan"),
                     "best_pair_ddm_sigma": float(bp[2]) if bp else float("nan"),
                     "best_pair_ddm": float(bp[1]) if bp else float("nan"),
                     "reasons": reason})
    ddm_sigma_all = np.array([x for x in ddm_sigma_all if np.isfinite(x)])
    theta_all = np.array(theta_all)
    mint = np.array([r["min_theta_arcsec"] for r in rows if np.isfinite(r["min_theta_arcsec"])])
    rec_sig = np.array([r["best_pair_ddm_sigma"] for r in rows
                        if r["recovered"] and np.isfinite(r["best_pair_ddm_sigma"])])
    return {
        "n_known": len(srcs), "n_recovered": int(sum(r["recovered"] for r in rows)),
        "miss_reasons": miss,
        "per_source_min_theta": {
            "n": int(mint.size),
            "lt30_frac": float(np.mean(mint < 30)) if mint.size else None,
            "lt600_frac": float(np.mean(mint < 600)) if mint.size else None,
            "median_arcsec": float(np.median(mint)) if mint.size else None},
        "recovered_best_pair_ddm_sigma": {
            "n": int(rec_sig.size),
            "median": float(np.median(rec_sig)) if rec_sig.size else None,
            "lt1_frac": float(np.mean(rec_sig < 1)) if rec_sig.size else None},
        "ddm_sigma_within600": {
            "n": int(ddm_sigma_all.size),
            "median": float(np.median(ddm_sigma_all)) if ddm_sigma_all.size else None,
            "lt1_frac": float(np.mean(ddm_sigma_all < 1)) if ddm_sigma_all.size else None,
            "p90": float(np.quantile(ddm_sigma_all, 0.9)) if ddm_sigma_all.size else None},
        "theta_all": {"median_arcsec": float(np.median(theta_all)) if theta_all.size else None,
                      "lt30_frac": float(np.mean(theta_all < 30)) if theta_all.size else None,
                      "lt600_frac": float(np.mean(theta_all < 600)) if theta_all.size else None},
        "rows": rows,
    }


def independence(pairs):
    nodes, deg = set(), {}
    parent = {}

    def find(x):
        parent.setdefault(x, x)
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for a, b in pairs:
        nodes |= {a, b}
        deg[a] = deg.get(a, 0) + 1
        deg[b] = deg.get(b, 0) + 1
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[ra] = rb
    comps = {}
    for x in nodes:
        comps.setdefault(find(x), set()).add(x)
    sizes = sorted((len(v) for v in comps.values()), reverse=True)
    return {"n_pairs": len(pairs), "n_unique_events": len(nodes),
            "max_event_degree": max(deg.values()) if deg else 0,
            "top_event": max(deg, key=deg.get) if deg else None,
            "n_components": len(comps), "components": [sorted(c) for c in comps.values()],
            "component_sizes": sizes}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--B", type=int, default=1000)
    ap.add_argument("--sections", type=str, default="all")
    ap.add_argument("--seed", type=int, default=20261009)
    a = ap.parse_args()
    secs = ({"perm", "grid", "strong", "known", "indep", "blink", "inject", "ring"}
            if a.sections == "all" else set(a.sections.split(",")))

    c1, c2 = pd.read_parquet(RAW1), pd.read_parquet(RAW2)
    ev, dmerr, snr = build_with_aux(c1, c2)
    m_oo = np.array([r == "" for r in ev.src_rep])
    ev_oo, dmerr_oo = subset(ev, m_oo), dmerr[m_oo]
    print(f"事件表 {len(ev)}（one-off 子空间 {len(ev_oo)}，DM err 有限 "
          f"{int(np.isfinite(dmerr_oo).sum())}）")

    res = {"built_at": dt.datetime.now().astimezone().isoformat(timespec="seconds"),
           "criteria": "用户 2026-10-09 检验协议；管线 = 事件层去重 v1.3（5014 / 3742 / 3665）",
           "logic_version": LOGIC_VERSION,
           "config": {"thetas": list(THETAS), "dm_windows": list(DMWS),
                      "n_shuffle": N_SHUFFLE, "dec_bins": DEC_BINS,
                      "n_eff_trials": N_EFF_TRIALS}}

    rng0 = np.random.default_rng(a.seed)
    real, real_best, real_cfgs = one_rep(ev_oo, dmerr_oo, rng0)
    real["n_sig_cfg"] = int(sum(c["sig"] for c in real_cfgs))
    real["sum_groups"] = int(sum(c["ngroups"] for c in real_cfgs))
    print(f"[真实] one-off 候选 {real['n_cand']} 对 | 显著档 {real['n_sig_cfg']} | "
          f"maxΔDM/σ {real['max_ddm_sigma']:.2f} | minθ {real['min_theta']:.1f}\"")
    res["real_oneoff"] = {**real, "real_cfgs": [
        {k: c[k] for k in ("theta", "dm_window", "obs", "bg_mean", "thr", "sig", "ngroups")}
        for c in real_cfgs]}

    if "perm" in secs:
        print(f"\n===== 一、随机零假设置换（B={a.B}）=====")
        res["perm"] = {}
        _off = {"ra": 1, "dm": 2, "ra_dm": 3, "full": 4}
        for kind, name in (("ra", "RA置换(Dec分层)"), ("dm", "DM置换(保分布洗牌)"),
                           ("ra_dm", "RA+DM置换"), ("full", "完全随机化")):
            r = run_permutation(kind, ev_oo, dmerr_oo, real, a.B, a.seed + _off[kind])
            res["perm"][kind] = {"name": name, **r}
            print(f"  {name:18s}: 候选均值 {r['cand_mean']:.2f} | 95分位 {r['cand_p95']:.0f} "
                  f"| max {r['cand_max']:.0f} | 真实{r['real_cand']}分位 {r['real_percentile']*100:.1f}% "
                  f"| 单次p {r['single_p']:.4f} | 全局p(Bonf×9) {r['global_p_bonferroni']:.4f}")

    if "grid" in secs:
        print(f"\n===== 二、搜索空间与多重比较 =====")

        def space_table(evx, cfgs, label):
            rows = []
            for c in cfgs:
                i, j = c["i"], c["j"]
                events = set()
                for aa, bb in zip(i.tolist(), j.tolist()):
                    events |= {str(evx.name[aa]), str(evx.name[bb])}
                hit28 = 0
                for v in real_best.values():
                    na, nb = str(ev_oo.name[v["a"]]), str(ev_oo.name[v["b"]])
                    if {na, nb} <= events:
                        hit28 += 1
                rows.append({"domain": label, "theta": c["theta"], "dm_window": c["dm_window"],
                             "obs": c["obs"], "bg_mean": round(c["bg_mean"], 2),
                             "thr": round(c["thr"], 2), "sig": c["sig"], "ngroups": c["ngroups"],
                             "n_unique_events": len(events), "n_of_28": hit28,
                             "n_geom": c["n_geom"], "n_dedup": c["n_dedup"]})
            return rows

        rngf = np.random.default_rng(a.seed + 1)
        cfgs_full = grid_search(ev, rngf)
        rows_full = space_table(ev, cfgs_full, "full")
        rows_oo = space_table(ev_oo, real_cfgs, "one_off")
        rngc = np.random.default_rng(a.seed + 2)
        cfgs_corr = grid_search(ev_oo, rngc, dedup_bg=True)
        rows_corr = space_table(ev_oo, cfgs_corr, "one_off_corrected_bg")
        res["grid"] = {"full": rows_full, "one_off": rows_oo, "one_off_corrected_bg": rows_corr}
        for g in rows_full + rows_oo:
            print(f"  [{g['domain']:8s}] θ={g['theta']:5.0f}\" ±{g['dm_window']:4.0f}: "
                  f"对 {g['obs']:5d} (背景 {g['bg_mean']:8.2f}) 显著 {str(g['sig']):5s} "
                  f"组 {g['ngroups']:3d} 唯一事件 {g['n_unique_events']:4d} 含28对 {g['n_of_28']:2d}")
        n_sig_frozen = sum(g["sig"] for g in rows_oo)
        n_sig_corr = sum(g["sig"] for g in rows_corr)
        print(f"  → 冻结口径显著档 {n_sig_frozen}/9；一致性修正背景后显著档 {n_sig_corr}/9")
        res["grid_summary"] = {"n_sig_frozen": n_sig_frozen, "n_sig_corrected_bg": n_sig_corr,
                               "corrected_consequence":
                               "背景去重一致性修正后阈值下降、更多档位显著 → 候选并集膨胀，"
                               "冻结口径 28 对为偏保守子集"}

    if "strong" in secs:
        print(f"\n===== 三、最强候选物理一致性 =====")
        c2r = pd.read_csv(RAW2CSV, low_memory=False)
        sub = c2r[c2r.tns_name.isin(["FRB20200511C", "FRB20211128E"])]
        rows = []
        for r in sub.itertuples():
            rows.append({"name": r.tns_name, "ra": r.ra, "dec": r.dec,
                         "ra_err": r.ra_err, "dec_err": r.dec_err,
                         "dm": r.dm_fitb, "dm_err": r.dm_fitb_err,
                         "snr_bonsai": r.bonsai_snr, "snr_fitb": r.snr_fitb,
                         "fluence": r.fluence, "fluence_err": r.fluence_err,
                         "width": r.width_fitb, "mjd": r.mjd_400})
        A, B_ = rows[0], rows[1]
        ddm = abs(A["dm"] - B_["dm"])
        e_off = float(np.hypot(A["dm_err"], B_["dm_err"]))
        sys_scan = {}
        for s in (0.5, 1.0):
            e = float(np.sqrt(A["dm_err"] ** 2 + B_["dm_err"] ** 2 + 2 * s ** 2))
            sys_scan[s] = {"sigma_comb": e, "ddm_sigma": ddm / e}
        th = haversine_arcsec(A["ra"], A["dec"], B_["ra"], B_["dec"])
        s_ra = float(np.hypot(A["ra_err"], B_["ra_err"])) * 3600
        s_dec = float(np.hypot(A["dec_err"], B_["dec_err"])) * 3600
        res["strong"] = {"events": rows, "ddm": ddm,
                         "sigma_official": e_off, "ddm_sigma_official": ddm / e_off,
                         "sigma_sys_scan": sys_scan, "theta_arcsec": th,
                         "sigma_ra_arcsec": s_ra, "sigma_dec_arcsec": s_dec,
                         "theta_over_sigma_1d": th / np.hypot(s_ra, s_dec),
                         "beam_column": "Cat2 导出无波束列（ra_dec_notes/exp_notes 均空）",
                         "external_localization": "无（CHIME 单站，无 FAST/Parkes/ASKAP/MeerKAT 亚角秒定位）"}
        print(f"  ΔDM={ddm:.4f}; σ官方={e_off:.4f} → {ddm/e_off:.2f}σ; "
              f"σ_sys=0.5 → {sys_scan[0.5]['ddm_sigma']:.2f}σ; "
              f"σ_sys=1.0 → {sys_scan[1.0]['ddm_sigma']:.2f}σ")
        print(f"  θ={th:.1f}\"; σ_ra≈{s_ra:.0f}\" σ_dec≈{s_dec:.0f}\" → θ/σ_1d="
              f"{th/np.hypot(s_ra,s_dec):.3f}")

    if "known" in secs:
        print(f"\n===== 四、已知重复源对照 =====")
        kr = known_repeater_analysis(ev, dmerr)
        res["known"] = kr
        print(f"  已知重复源 {kr['n_known']}；召回 {kr['n_recovered']}；漏检原因 {kr['miss_reasons']}")
        print(f"  [每源最近对] min θ：中位 {kr['per_source_min_theta']['median_arcsec']:.0f}\" "
              f"<30\" 占比 {kr['per_source_min_theta']['lt30_frac']*100:.0f}% "
              f"<600\" 占比 {kr['per_source_min_theta']['lt600_frac']*100:.0f}%")
        print(f"  [已召回源最近对] ΔDM/σ 中位 "
              f"{kr['recovered_best_pair_ddm_sigma']['median']:.2f} <1σ 占比 "
              f"{kr['recovered_best_pair_ddm_sigma']['lt1_frac']*100:.0f}%")
        print(f"  [源内 600\" 内全部对] ΔDM/σ n={kr['ddm_sigma_within600']['n']} 中位 "
              f"{kr['ddm_sigma_within600']['median']:.2f}")

    v1 = pd.read_parquet(CAND1)
    pairs = list(zip(v1.name1.astype(str), v1.name2.astype(str)))
    if "indep" in secs:
        print(f"\n===== 五、候选对独立性 =====")
        ind = independence(pairs)
        res["independence"] = ind
        print(f"  28 对 → 唯一事件 {ind['n_unique_events']} | 最大出现次数 "
              f"{ind['max_event_degree']}（{ind['top_event']}）| 连通分量 "
              f"{ind['n_components']} 大小 {ind['component_sizes']}")

    if "blink" in secs:
        print(f"\n===== 六、Blinkverse 零重复 N_exp =====")
        bv = json.loads(BVJSON.read_text(encoding="utf-8"))["data"]
        by = {}
        for x in bv:
            c = x.get("content", {})
            by[str(c.get("source", "")).upper()] = c
        st = {n: (by.get(n, {}) or {}).get("type", "未收录")
              for n in ("FRB20200511C", "FRB20211128E")}
        t_eff_h = 179.73
        lam = 0.0111
        eff = 0.261
        res["blinkverse"] = {
            "status": st,
            "fluence_jy_ms": 3.4056, "width_ms": 2.885, "dm": 1835.04,
            "assumptions": {"T_eff_h_cat2_epoch": t_eff_h, "lambda_per_hour": lam,
                            "fluence_efficiency": eff,
                            "note": "Blinkverse 为聚合库（4020 源），自身不监测；"
                                    "N_exp = T_eff × λ × eff × f_beam，f_beam=1（单波束在源中心）"},
            "N_exp_chime_cat2_epoch": t_eff_h * lam * eff,
            "N_exp_per_100h_any_monitor": 100.0 * lam * eff,
        }
        print(f"  Blinkverse 状态 {st}")
        print(f"  N_exp(CHIME Cat2 纪元 T_eff={t_eff_h}h) = {t_eff_h*lam*eff:.2f}")

    if "inject" in secs:
        print(f"\n===== 七、注入-回收与假阳性率 =====")
        res["inject"] = injection_battery(ev_oo, dmerr_oo, real, a.B, a.seed)
        print(f"  随机 one-off 平均候选 {res['inject']['random_oneoff']['mean']:.2f} "
              f"| 95分位 {res['inject']['random_oneoff']['p95']:.0f}")
        print(f"  真重复注入回收率（θ=600,ΔDMdrift=1）："
              f"{res['inject']['true_repeater']['recovery']:.2f}")

    if "ring" in secs:
        print(f"\n===== 八、环形平移/单簇宽窗 p_floor（180916）=====")
        res["ringshift"] = ring_shift_180916(a.seed)
        print(f"  {res['ringshift']['summary']}")

    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "c2_null_battery.json").write_text(
        json.dumps(res, indent=2, ensure_ascii=False, default=str), encoding="utf-8")
    print(f"\n→  {OUT / 'c2_null_battery.json'}")
    return 0


def injection_battery(ev_oo, dmerr_oo, real, B, seed):
    rng = np.random.default_rng(seed + 7)
    fin = np.where(np.isfinite(ev_oo.dm))[0]
    nf = fin.size
    decf, mjdf, dmf, def_ = ev_oo.dec[fin], ev_oo.mjd[fin], ev_oo.dm[fin], dmerr_oo[fin]

    def synth_catalog(m, rng):
        di = rng.integers(0, nf, m)
        mi = rng.integers(0, nf, m)
        dmi = rng.integers(0, nf, m)
        ra = rng.uniform(0.0, 360.0, size=m)
        names = np.array([f"SYN{k:05d}" for k in range(m)])
        evt = EventTable(ra, decf[di], mjdf[mi], dmf[dmi], names,
                         np.array(["cat2"] * m), np.array([""] * m),
                         np.zeros(m, bool), np.zeros(m, bool))
        return evt, def_[dmi]

    fp = np.empty(B)
    for b in range(B):
        evs, ds = synth_catalog(nf, rng)
        st, _, _ = one_rep(evs, ds, rng)
        fp[b] = st["n_cand"]

    dec0 = float(np.nanmedian(ev_oo.dec))
    mjd0 = float(np.nanmedian(ev_oo.mjd))
    dm0 = float(np.nanmedian(ev_oo.dm))

    def simulate_source(k, sig_pos_arcsec, drift, span_days, th_test, dw_test, rng):
        ra0 = 180.0
        ra = ra0 + sig_pos_arcsec / 3600.0 * rng.standard_normal(k) / np.cos(np.radians(dec0))
        dec = dec0 + sig_pos_arcsec / 3600.0 * rng.standard_normal(k)
        dm = dm0 + rng.uniform(-drift, drift, size=k)
        mjd = mjd0 + rng.uniform(0, span_days, size=k)
        evs = EventTable(ra, dec, mjd, dm, np.array([f"S{k}_{i}" for i in range(k)]),
                         np.array(["cat2"] * k), np.array(["SYN"] * k),
                         np.zeros(k, bool), np.zeros(k, bool))
        i, j = cone_pairs(evs, th_test, dw_test)
        return i.size > 0

    grid = []
    for th_test in (60.0, 180.0, 600.0):
        for dw_test in (1.0, 5.0):
            for drift in (0.0, 1.0, 2.0, 5.0):
                Nsim = 300
                ok = sum(simulate_source(3, 30.0, drift, 100.0, th_test, dw_test, rng)
                         for _ in range(Nsim))
                grid.append({"theta": th_test, "dm_window": dw_test, "dm_drift": drift,
                             "recovery": ok / Nsim, "n_sim": Nsim})
    base = [g for g in grid if g["theta"] == 600.0 and g["dm_window"] == 5.0
            and g["dm_drift"] == 1.0][0]

    cut = float(np.quantile(ev_oo.mjd, 0.70))
    tr, ho = ev_oo.mjd <= cut, ev_oo.mjd > cut
    ev_tr, ev_ho = subset(ev_oo, tr), subset(ev_oo, ho)
    rng2 = np.random.default_rng(seed + 11)
    ctr = grid_search(ev_tr, rng2)
    bg_map = {(c["theta"], c["dm_window"]): c["bg_mean"] for c in ctr}
    thr_map = {(c["theta"], c["dm_window"]): c["thr"] for c in ctr}
    c_ho = grid_search(ev_ho, np.random.default_rng(seed + 12), n_shuffle=50)
    ho_rec = [{"theta": c["theta"], "dm_window": c["dm_window"], "obs": c["obs"],
               "train_thr": round(thr_map[(c["theta"], c["dm_window"])], 2),
               "train_bg": round(bg_map[(c["theta"], c["dm_window"])], 2),
               "sig_holdout": c["sig"]} for c in c_ho]
    extra = np.empty(200)
    for b in range(200):
        evs, ds = synth_catalog(500, rng2)
        merged = EventTable(
            np.concatenate([ev_ho.ra, evs.ra]), np.concatenate([ev_ho.dec, evs.dec]),
            np.concatenate([ev_ho.mjd, evs.mjd]), np.concatenate([ev_ho.dm, evs.dm]),
            np.concatenate([ev_ho.name, evs.name]), np.concatenate([ev_ho.catalog, evs.catalog]),
            np.concatenate([ev_ho.src_rep, evs.src_rep]),
            np.concatenate([ev_ho.multi_component, evs.multi_component]),
            np.concatenate([ev_ho.excluded, evs.excluded]))
        cfgs = grid_search(merged, rng2, n_shuffle=30)
        best = union_candidates(cfgs, merged)
        extra[b] = sum(1 for k in best if k[0].startswith("SYN") or k[1].startswith("SYN"))
    return {"random_oneoff": {"mean": float(fp.mean()), "p95": float(np.quantile(fp, 0.95)),
                              "max": float(fp.max()), "n_reps": B},
            "true_repeater": {"recovery": base["recovery"], "grid": grid},
            "blind_holdout": {"n_train": int(tr.sum()), "n_hold": int(ho.sum()),
                              "holdout_sig_cfgs": int(sum(1 for c in c_ho if c["sig"])),
                              "per_config": ho_rec,
                              "neg_ctrl_synth_involved_pairs_mean": float(extra.mean()),
                              "neg_ctrl_synth_involved_pairs_max": float(extra.max())}}


def ring_shift_180916(seed, n_perm=1000):
    raw = pd.read_parquet(RAW2)
    rep = clean_rep(raw.repeater_name)
    m = np.array([r == "FRB20180916B" for r in rep])
    t = np.sort(pd.to_numeric(raw.mjd_400[m], errors="coerce").dropna().to_numpy(float))
    if t.size < 5:
        return {"summary": "无 FRB20180916B 事件", "n": int(t.size)}
    P, w = 16.33128518622877, 0.25
    cut = float(np.quantile(t, 0.70))
    train, hold = t[t <= cut], t[t > cut]
    t0 = float(train[0])
    tol_days = max(0.25, 0.5 * w * P)
    tol_ph = tol_days / P
    ph_tr = np.mod((train - t0) / P, 1.0)
    best = max((int(np.sum(np.mod(ph_tr - p0 + tol_ph, 1.0) < (w + 2 * tol_ph))), p0)
               for p0 in np.arange(0, 1, 0.01))
    phi0 = float(best[1])
    ph_ho = np.mod((hold - t0) / P, 1.0)
    hit = np.mod(ph_ho - phi0 + tol_ph, 1.0) < (w + 2 * tol_ph)
    n_hit = int(hit.sum())
    p_eff = min(1.0, w + 2 * tol_ph)
    rng = np.random.default_rng(seed + 21)
    counts = np.empty(n_perm, dtype=int)
    for i in range(n_perm):
        shift = rng.uniform(0.0, P)
        ph2 = np.mod((hold + shift - t0) / P, 1.0)
        counts[i] = int(np.sum(np.mod(ph2 - phi0 + tol_ph, 1.0) < (w + 2 * tol_ph)))
    p_mc = (1 + int(np.sum(counts >= n_hit))) / (n_perm + 1)
    return {"src": "FRB20180916B", "P": P, "w": w, "tol0_days": tol_days, "drift_days": 0.0,
            "p_eff": p_eff, "n_total": int(t.size), "n_train": int(train.size),
            "n_hold": int(hold.size), "hits": n_hit, "rate": n_hit / hold.size,
            "phi0_window": phi0,
            "ringshift_mean_hits": float(counts.mean()),
            "ringshift_p95_hits": float(np.quantile(counts, 0.95)),
            "p_mc": p_mc, "g8_frozen": {"n_hold": 43, "hits": 43, "E": 22.19,
                                        "p_mc": 0.2495, "z_mc": 1.23},
            "summary": (f"命中 {n_hit}/{hold.size}; w={w}; p_eff=w+2·tol/P={p_eff:.3f}; "
                        f"ring-shift 均值 {counts.mean():.1f} 命中; p_mc={p_mc:.3f} "
                        f"(G8 冻结 p_mc=0.2495≈w) → p 下限≈覆盖率，不可能 <0.05")}


if __name__ == "__main__":
    raise SystemExit(main())
