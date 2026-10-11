"""Full ranking table: pair-level evidence columns plus union-find components.

Usage: python scripts/rank_c2_candidates.py
Output: data/processed/upgrade/c2_ranked_{crossday,sameday,hypotheses}_v1.{parquet,csv}
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from astropy import units as u
from astropy.coordinates import SkyCoord

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw" / "chime" / "cat2_25" / "chimefrbcat2.csv"
V1 = ROOT / "data" / "processed" / "upgrade" / "c2_candidates_v1.parquet"
V2 = ROOT / "data" / "processed" / "upgrade" / "c2_candidates_v2.parquet"
OUTD = ROOT / "data" / "processed" / "upgrade"

FLAGCOLS = ["excluded_flag", "sidelobe_flag", "intrachan_flag", "citizen_science_flag"]


def clean_rep(s):
    r = s.astype("string").str.strip()
    return r.fillna("").mask(r.isin(["-9999", "nan", "NaN", "None"]), "").to_numpy(dtype=object)


def _truthy(v) -> bool:
    return str(v).strip().lower() in ("1", "1.0", "true")


def build_event_frame(raw, names):
    rows = {}
    for nm in names:
        g = raw[raw.tns_name == nm]
        sub = pd.to_numeric(g.sub_num, errors="coerce")
        pr = g[sub.isna() | (sub == 0)]
        r = (pr if len(pr) else g).iloc[0]
        d = {"name": nm, "mjd": float(r.mjd_400),
             "ra": float(r.ra), "dec": float(r.dec),
             "dm": float(r.dm_fitb), "dm_err": float(r.dm_fitb_err),
             "bonsai_snr": float(r.bonsai_snr), "snr_fitb": float(pd.to_numeric(r.snr_fitb)),
             "flag_frac": float(r.flag_frac) if pd.notna(r.flag_frac) else np.nan,
             "n_rows_in_cat": int(len(g)),
             "src_rep": clean_rep(pd.Series([r.repeater_name]))[0]}
        fl = [c for c in FLAGCOLS if _truthy(r[c])]
        d["flags"] = "|".join(fl) if fl else "-"
        d["excluded"] = int("excluded_flag" in fl)
        rows[nm] = d
    return pd.DataFrame(rows.values()).set_index("name")


def nearest_repeater(raw, ev_names):
    d = raw[raw.ra.notna() & raw.dec.notna()].copy()
    d["src_rep"] = clean_rep(d["repeater_name"])
    srcs = {k: SkyCoord(v.ra.to_numpy(float) * u.deg, v.dec.to_numpy(float) * u.deg)
            for k, v in d[d.src_rep != ""].groupby("src_rep") if len(v) >= 2}
    out = {}
    for nm in ev_names:
        tp = SkyCoord(float(raw.loc[raw.tns_name == nm, "ra"].iloc[0]) * u.deg,
                      float(raw.loc[raw.tns_name == nm, "dec"].iloc[0]) * u.deg)
        best, bname = np.inf, ""
        for name, c in srcs.items():
            s = float(c.separation(tp).arcsec.min())
            if s < best:
                best, bname = s, name
        out[nm] = (best, bname)
    return out, len(srcs)


def components(pairs):
    parent = {}

    def find(x):
        parent.setdefault(x, x)
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for a, b in pairs:
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[ra] = rb
    grp = {}
    for x in list(parent):
        grp.setdefault(find(x), []).append(x)
    name2cid, cid2size = {}, {}
    for i, (_, ms) in enumerate(sorted(grp.items(), key=lambda kv: -len(kv[1])), 1):
        cid2size[i] = len(ms)
        for nm in ms:
            name2cid[nm] = i
    return name2cid, cid2size


def main() -> int:
    raw = pd.read_csv(RAW, low_memory=False)
    v1 = pd.read_parquet(V1)
    v2 = pd.read_parquet(V2)
    keep = ["name1", "name2", "bv_status1", "bv_status2", "verdict"]
    p = v1.merge(v2[keep], on=["name1", "name2"], how="left", validate="one_to_one")

    names = sorted(set(p.name1) | set(p.name2))
    ev = build_event_frame(raw, names)
    nr, n_src = nearest_repeater(raw, names)
    comp, sizes = components(list(zip(p.name1, p.name2)))
    p["comp_id"] = [comp[n] for n in p.name1]
    p["comp_size"] = [sizes[comp[n]] for n in p.name1]
    for side in (1, 2):
        nm = p[f"name{side}"]
        p[f"bonsai_snr{side}"] = nm.map(ev["bonsai_snr"]).astype(float)
        p[f"snr_fitb{side}"] = nm.map(ev["snr_fitb"]).astype(float)
        p[f"flag_frac{side}"] = nm.map(ev["flag_frac"]).astype(float)
        p[f"flags{side}"] = nm.map(ev["flags"])
        p[f"dm_err{side}"] = nm.map(ev["dm_err"]).astype(float)
        p[f"low_snr{side}"] = p[f"bonsai_snr{side}"] < 10
        p[f"rep_nearest_arcsec{side}"] = nm.map(lambda x: round(nr[x][0], 1))
        p[f"rep_nearest{side}"] = nm.map(lambda x: nr[x][1])
        p[f"rep_within_1deg{side}"] = p[f"rep_nearest_arcsec{side}"] < 3600
    p["any_flagged"] = p["flags1"].ne("-") | p["flags2"].ne("-")
    p["any_low_snr"] = p.low_snr1 | p.low_snr2
    p["any_near_repeater"] = p.rep_within_1deg1 | p.rep_within_1deg2
    p["is_cross_day"] = p.dt_days >= 1.0

    cols = ["comp_id", "comp_size", "cls", "name1", "name2", "dt_days", "theta_arcsec",
            "ddm", "ddm_err", "ddm_sigma",
            "bonsai_snr1", "bonsai_snr2", "snr_fitb1", "snr_fitb2",
            "flag_frac1", "flag_frac2", "flags1", "flags2",
            "rep_nearest_arcsec1", "rep_nearest1",
            "rep_nearest_arcsec2", "rep_nearest2",
            "bv_status1", "bv_status2", "verdict"]

    def ranked(sub):
        return sub.sort_values(["ddm_sigma", "theta_arcsec", "dt_days"]).round(3)[cols]

    cross = ranked(p[p.is_cross_day])
    same = ranked(p[~p.is_cross_day])

    hy = []
    for cid, g in p.groupby("comp_id"):
        gc = g[g.is_cross_day]
        pick = gc if len(gc) else g
        hy.append(dict(
            comp_id=int(cid), n_nodes=int(g.comp_size.iloc[0]),
            n_pairs=len(g), n_cross=int(g.is_cross_day.sum()),
            members=",".join(sorted(set(g.name1) | set(g.name2))),
            best_ddm_sigma=float(pick.ddm_sigma.min()),
            best_theta_arcsec=float(pick.theta_arcsec.min()),
            dt_range=[float(g.dt_days.min()), float(g.dt_days.max())],
            min_bonsai_snr=float(min(pick.bonsai_snr1.min(), pick.bonsai_snr2.min())),
            max_flag_frac=float(np.nanmax([pick.flag_frac1.max(), pick.flag_frac2.max()])),
            any_flagged=bool(pick.any_flagged.any()),
            any_near_repeater=bool(pick.any_near_repeater.any()),
            verdict=str(pick.iloc[np.argmin(pick.ddm_sigma.values)].verdict)))
    hyd = pd.DataFrame(hy).sort_values("best_ddm_sigma").reset_index(drop=True)
    hyd["contaminated"] = hyd.comp_id.isin(
        p.loc[(p.rep_nearest_arcsec1 < 600) | (p.rep_nearest_arcsec2 < 600), "comp_id"].unique())
    hyd["nearest_rep"] = [", ".join(sorted(set(
        p.loc[p.comp_id == c, ["rep_nearest1", "rep_nearest2"]].values.ravel())
        - {"", "nan"})) for c in hyd.comp_id]
    surv = hyd[~hyd.contaminated]

    cross.to_parquet(OUTD / "c2_ranked_crossday_v1.parquet", index=False)
    cross.to_csv(OUTD / "c2_ranked_crossday_v1.csv", index=False)
    same.to_parquet(OUTD / "c2_ranked_sameday_v1.parquet", index=False)
    hyd.to_parquet(OUTD / "c2_ranked_hypotheses_v1.parquet", index=False)
    hyd.to_csv(OUTD / "c2_ranked_hypotheses_v1.csv", index=False)

    pd.set_option("display.width", 260)
    print(f"候选 {len(p)} 对 / 事件 {len(names)} / 连通分量（独立假设）= {len(sizes)}；"
          f"重复源对照集 = {n_src} 个")
    print(f"两个 SNR 口径差异实测：bonsai 中位 {p[['bonsai_snr1','bonsai_snr2']].values.mean():.1f} "
          f"vs fitb 中位 {p[['snr_fitb1','snr_fitb2']].values.mean():.1f}")
    print("\n######## 跨日档（真升级候选）按 ΔDM/σ 升序 ########")
    print(cross.to_string(index=False))
    print("\n######## 同日档（同一风暴内部两两组合 = 非独立假设）########")
    print(same.to_string(index=False))
    print("\n######## 假设级汇总（独立检验单位）########")
    print(hyd.to_string(index=False))
    print("\n######## 存活的新链接假设（剔除已编目重复源污染）########")
    print(surv.to_string(index=False))
    print("\n== 层检查汇总 ==")
    print(f"  含 excluded/sidelobe/intrachan 成员的对 = {int(p.any_flagged.sum())}/{len(p)}")
    print(f"  含 bonsai_snr<10 成员的对 = {int(p.any_low_snr.sum())}/{len(p)}"
          f"（若改用 snr_fitb<10 则 = {int(((p.snr_fitb1 < 10) | (p.snr_fitb2 < 10)).sum())}）")
    print(f"  有成员距某已编目重复源 <1° 的对 = {int(p.any_near_repeater.sum())}/{len(p)}"
          f"（<600″ 的对 = {int(((p.rep_nearest_arcsec1 < 600) | (p.rep_nearest_arcsec2 < 600)).sum())}）")
    n_oo = int((clean_rep(raw["repeater_name"]) == "").sum())
    cone = np.pi * (600 / 3600.0) ** 2
    print(f"\n== 机会符合背景 ==")
    print(f"  one-off 行 {n_oo} × 600″ 锥面积 {cone:.5f} deg² / 23953 deg² "
          f"→ 期望 ≈ {n_oo * cone / 23953:.4f} 个事件落在**指定某一个**重复源附近")
    print(f"  受污染分量 = {int(hyd.contaminated.sum())}/{len(hyd)}"
          f"（涉及对 {int(p.comp_id.isin(hyd.loc[hyd.contaminated, 'comp_id']).sum())}/{len(p)}）")
    print(f"  → 存活的新链接假设 = {len(surv)}；存活的对 = "
          f"{int(p.comp_id.isin(surv.comp_id).sum())}")
    (OUTD / "c2_ranked_manifest.json").write_text(json.dumps(dict(
        n_pairs=len(p), n_events=len(names), n_components=len(sizes),
        n_repeater_sources=n_src, sort_key="ddm_sigma asc, theta asc, dt asc",
        split="cross_day / same_day 分表", low_snr_definition="bonsai_snr < 10",
        inputs=[str(V1.name), str(V2.name), str(RAW.name)]), ensure_ascii=False, indent=2),
        encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
