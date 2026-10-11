"""Known-repeater control statistics with explicit pool definitions.

Usage: python scripts/known_repeater_control.py
Output: data/processed/upgrade/known_repeater_control{.json,_rows.csv}
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from dm_sigma_floor import local_scale          # noqa: E402
from rank_c2_candidates import clean_rep        # noqa: E402

OUT = ROOT / "data" / "processed" / "upgrade"
REF_FRAC = 0.25
CAND = dict(name="FRB20200511C x FRB20211128E", dm1=1835.0404, dm2=1834.8635,
            e1=0.327129, e2=0.378592, theta=23.9)


def burst_pool_cat2():
    raw = pd.read_csv(ROOT / "data/raw/chime/cat2_25/chimefrbcat2.csv", low_memory=False)
    raw["src"] = clean_rep(raw["repeater_name"])
    raw["m"] = pd.to_numeric(raw.mjd_400, errors="coerce")
    for c in ("ra", "dec", "dm_fitb", "dm_fitb_err", "bonsai_dm", "event_id"):
        raw[c] = pd.to_numeric(raw[c], errors="coerce")
    d = raw[(raw.src != "") & raw.ra.notna() & raw.dec.notna()
            & raw.dm_fitb.notna() & raw.m.notna()].copy()
    d = d.drop_duplicates(["src", "event_id"])
    return d


def floor_s(ref_dm, ref_D, dm_target):
    return local_scale(dm_target, ref_dm, ref_D, REF_FRAC)[0]


class FloorGrid:

    def __init__(self, ref_dm, ref_D, n_grid=220, lo=50.0, hi=2400.0):
        self.g = np.exp(np.linspace(np.log(lo), np.log(hi), n_grid))
        self.v = np.array([floor_s(ref_dm, ref_D, float(x)) for x in self.g], float)
        ok = np.isfinite(self.v)
        self.g, self.v = self.g[ok], self.v[ok]

    def __call__(self, dm_target):
        if not np.isfinite(dm_target):
            return np.nan
        return float(np.interp(dm_target, self.g, self.v))


def haversine_arcsec(ra1, de1, ra2, de2):
    r = np.radians
    a = np.sin(r(de2 - de1) / 2) ** 2 + np.cos(r(de1)) * np.cos(r(de2)) * np.sin(r(ra2 - ra1)) ** 2
    return 2 * np.arcsin(np.sqrt(a)) * 3600 / (np.pi / 180)


def per_source_stats(d, floor):
    out = []
    for s, g in d.groupby("src"):
        ra, de = g.ra.to_numpy(float), g.dec.to_numpy(float)
        dm, err = g.dm_fitb.to_numpy(float), g.dm_fitb_err.to_numpy(float)
        n = len(g)
        if n < 2:
            continue
        pairs = []
        for i in range(n):
            for j in range(i + 1, n):
                th = haversine_arcsec(ra[i], de[i], ra[j], de[j])
                ddm = abs(dm[i] - dm[j])
                ef = float(np.hypot(err[i], err[j]))
                sfl = floor(0.5 * (dm[i] + dm[j]))
                el = float(np.sqrt(2.0) * sfl) if np.isfinite(sfl) and sfl > 0 else np.nan
                pairs.append(dict(theta=th, ddm=ddm,
                                  rF=ddm / ef if ef > 0 else np.nan,
                                  rL=ddm / el if np.isfinite(el) and el > 0 else np.nan,
                                  dm_mean=0.5 * (dm[i] + dm[j])))
        if not pairs:
            continue
        pa = min(pairs, key=lambda p: p["theta"])
        pbF = min(pairs, key=lambda p: p["rF"] if np.isfinite(p["rF"]) else 9e9)
        pbL = min(pairs, key=lambda p: p["rL"] if np.isfinite(p["rL"]) else 9e9)
        recovered = any(p["theta"] < 600 and p["ddm"] < 5.0 for p in pairs)
        out.append(dict(src=s, n=n, recovered=bool(recovered),
                        A_theta=pa["theta"], A_ddm=pa["ddm"], A_rF=pa["rF"], A_rL=pa["rL"],
                        A_dmmean=pa["dm_mean"],
                        B_rF=pbF["rF"], B_theta=pbF["theta"],
                        B_rL=pbL["rL"], B_theta_L=pbL["theta"]))
    return pd.DataFrame(out)


def summarize(t, label):
    def st(col):
        v = pd.to_numeric(t[col], errors="coerce").dropna().to_numpy(float)
        return dict(n=int(v.size), median=round(float(np.median(v)), 3),
                    lt1_frac=round(float(np.mean(v < 1)), 4),
                    p90=round(float(np.quantile(v, .9)), 3))
    rec = t[t.recovered]
    return {
        "pool": label, "n_sources": int(len(t)), "n_recovered": int(t.recovered.sum()),
        "recovery_rate": round(float(t.recovered.mean()), 4),
        "A_thetaClosest_formal": st("A_rF"), "A_thetaClosest_floor": st("A_rL"),
        "B_minRatio_formal": st("B_rF"), "B_minRatio_floor": st("B_rL"),
        "min_theta_arcsec": {"median": round(float(t.A_theta.median()), 1),
                             "lt30_frac": round(float((t.A_theta < 30).mean()), 4),
                             "lt600_frac": round(float((t.A_theta < 600).mean()), 4)},
        "A_formal_recoveredOnly": st_on(rec, "A_rF"), "A_floor_recoveredOnly": st_on(rec, "A_rL"),
    }


def st_on(t, col):
    v = pd.to_numeric(t[col], errors="coerce").dropna().to_numpy(float)
    if not v.size:
        return {"n": 0}
    return dict(n=int(v.size), median=round(float(np.median(v)), 3),
                lt1_frac=round(float(np.mean(v < 1)), 4))


def main() -> int:
    raw = pd.read_csv(ROOT / "data/raw/chime/cat2_25/chimefrbcat2.csv", low_memory=False)
    ref_dm = pd.to_numeric(raw.dm_fitb, errors="coerce").to_numpy(float)
    ref_b = pd.to_numeric(raw.bonsai_dm, errors="coerce").to_numpy(float)
    ref_D = ref_dm - ref_b

    d = burst_pool_cat2()
    floor = FloorGrid(ref_dm, ref_D)
    t = per_source_stats(d, floor)
    S = summarize(t, "P1 Cat2 localized, dedup by event_id")

    ef = float(np.hypot(CAND["e1"], CAND["e2"]))
    sfl = floor(0.5 * (CAND["dm1"] + CAND["dm2"]))
    el = float(np.sqrt(2.0) * sfl)
    ddmc = abs(CAND["dm1"] - CAND["dm2"])
    candF, candL = ddmc / ef, ddmc / el
    vF = pd.to_numeric(t.A_rF, errors="coerce").dropna().to_numpy(float)
    vL = pd.to_numeric(t.A_rL, errors="coerce").dropna().to_numpy(float)
    S["candidate"] = dict(theta=CAND["theta"], ddm=round(ddmc, 4),
                          sigma_formal=round(ef, 4), sigma_floor=round(el, 4),
                          r_formal=round(candF, 4), r_floor=round(candL, 4),
                          percentile_formal=round(float(np.mean(vF < candF)), 4),
                          percentile_floor=round(float(np.mean(vL < candL)), 4))
    S["note_percentile_direction"] = ("percentile = 分布中比该候选更松的比例; "
                                      "越小表示该候选比典型真重复源更紧")

    rep_names = set(clean_rep(raw["repeater_name"])) - {""}
    cat1 = pd.read_parquet(ROOT / "data/interim/cat1_bursts.parquet")
    c1rep = set(clean_rep(cat1["repeater_name"] if "repeater_name" in cat1 else cat1["src"])) - {""}
    only1 = sorted(c1rep - rep_names - {""})
    S["positive_control_denominators"] = dict(
        cat2_sources=int(len(rep_names)), cat1_only_sources=only1,
        total_85=int(len(rep_names) + len(only1)),
        note="85 = 83 (Catalog 2) + 2 (Catalog 1 only)，与 §3.3 分母构造一致")

    S["positive_control_miss_reasons"] = dict(
        counts=miss_reasons(d, floor),
        convention="事件层池（与候选搜索同池：28 对是事件级去重结果，非行级 44）",
        note=("38 个漏检中 37 个因源内无任何 θ<600″ 的对，恰 1 个因 θ<600″ 内 ΔDM 超 ±5；"
              "故'全部漏检均因定位散布、无一因 DM'须软化为'37/38 因定位散布'。"
              "行级口径召回 61/83、合并去重口径 55/83，差异仅来自是否计入子暴行。"))

    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "known_repeater_control.json").write_text(
        json.dumps(S, ensure_ascii=False, indent=2), encoding="utf-8")
    t.to_csv(OUT / "known_repeater_control_rows.csv", index=False)
    print(json.dumps(S, ensure_ascii=False, indent=2))
    return 0


def miss_reasons(d, floor):
    import collections
    cnt = collections.Counter()
    for s, g in d.groupby("src"):
        ra, de = g.ra.to_numpy(float), g.dec.to_numpy(float)
        dm = g.dm_fitb.to_numpy(float)
        n = len(g)
        if n < 2:
            cnt["目录内事件<2（无法配对）"] += 1
            continue
        th = np.array([haversine_arcsec(ra[i], de[i], ra[j], de[j])
                       for i in range(n) for j in range(i + 1, n)])
        dd = np.array([abs(dm[i] - dm[j])
                       for i in range(n) for j in range(i + 1, n)])
        if np.any((th < 600) & (dd < 5.0)):
            cnt["已召回"] += 1
        elif not np.any(th < 600):
            cnt["角距>600″（逐暴定位散布）"] += 1
        elif not np.any((th < 600) & (dd < 5.0)):
            cnt["DM 差>±5 容差"] += 1
        else:
            cnt["其他"] += 1
    return dict(cnt)


if __name__ == "__main__":
    raise SystemExit(main())
