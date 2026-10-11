"""Recompute Delta-DM consistency with the empirical Cat2 DM error floor as denominator.

Usage: python scripts/dm_sigma_floor.py
Output: data/processed/upgrade/c2_ranked_dm_sigma_v1.{csv,parquet,md}, c2_ranked_*_v14.{parquet,csv}
"""
from __future__ import annotations

import datetime as dt
import json
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw" / "chime" / "cat2_25" / "chimefrbcat2.csv"
RANKED = ROOT / "data" / "processed" / "upgrade" / "c2_ranked_crossday_v1.parquet"
RANKED_SD = ROOT / "data" / "processed" / "upgrade" / "c2_ranked_sameday_v1.parquet"
OUTD = ROOT / "data" / "processed" / "upgrade"
FRACS = (0.15, 0.25, 0.40)
MIN_N = 60
LOGIC_VERSION = "1.4.0"


def local_scale(dm_target, ref_dm, ref_D, frac):
    if not np.isfinite(dm_target) or dm_target <= 0:
        return np.nan, np.nan, 0
    ok = np.isfinite(ref_dm) & np.isfinite(ref_D) & (ref_dm > 0)
    a, b = ref_dm[ok], ref_D[ok]
    f = float(frac)
    for _ in range(8):
        sel = np.abs(a - dm_target) <= f * dm_target
        n = int(sel.sum())
        if n >= MIN_N:
            sub = b[sel]
            ctr = np.median(sub)
            return float(1.4826 * np.median(np.abs(sub - ctr))), f, n
        f *= 1.5
    return np.nan, f, int((np.abs(a - dm_target) <= f * dm_target).sum())


def main() -> int:
    raw = pd.read_csv(RAW, low_memory=False)
    b = pd.to_numeric(raw.bonsai_dm, errors="coerce").to_numpy(float)
    f = pd.to_numeric(raw.dm_fitb, errors="coerce").to_numpy(float)
    ref_dm, ref_D = f, f - b
    m = np.isfinite(ref_dm) & np.isfinite(ref_D)
    print(f"全 Cat2 参考集：可比事件 n={int(m.sum())}；"
          f"|D| 全表中位={np.median(np.abs(ref_D[m])):.3f}  "
          f"形式误差中位={np.nanmedian(pd.to_numeric(raw.dm_fitb_err, errors='coerce')):.3f}")

    def member_dm(name, col):
        g = raw[raw.tns_name == name]
        if not len(g):
            return np.nan
        sub = pd.to_numeric(g.sub_num, errors="coerce")
        pr = g[sub.isna() | (sub == 0)]
        r = (pr if len(pr) else g).iloc[0]
        v = pd.to_numeric(r[col], errors="coerce")
        return float(v) if pd.notna(v) else np.nan

    cross = pd.read_parquet(RANKED)
    if RANKED_SD.exists():
        sameday = pd.read_parquet(RANKED_SD)
    else:
        sameday = pd.DataFrame()

    for lab, df in (("cross", cross), ("same", sameday)):
        if not len(df):
            continue
        for side in (1, 2):
            df[f"fitb_dm{side}"] = [member_dm(n, "dm_fitb") for n in df[f"name{side}"]]
            df[f"bonsai_dm{side}"] = [member_dm(n, "bonsai_dm") for n in df[f"name{side}"]]
            df[f"D{side}"] = df[f"fitb_dm{side}"] - df[f"bonsai_dm{side}"]
        df["dm_mean"] = df[["fitb_dm1", "fitb_dm2"]].mean(axis=1)
    if cross["fitb_dm1"].isna().any():
        raise SystemExit("有成员事件在 Cat2 原始表取不到 dm_fitb —— 中止，不猜数")

    REF_FRAC = 0.25
    for fr in FRACS:
        res = cross["dm_mean"].map(lambda x: local_scale(x, ref_dm, ref_D, fr))
        cross[f"s_F{int(fr * 100)}"] = [r[0] for r in res]
        cross[f"usedfrac_F{int(fr * 100)}"] = [r[1] for r in res]
        cross[f"nlocal_F{int(fr * 100)}"] = [r[2] for r in res]
        cross[f"sigB_F{int(fr * 100)}"] = np.sqrt(2.0) * cross[f"s_F{int(fr * 100)}"]
        cross[f"sigA_F{int(fr * 100)}"] = cross[f"s_F{int(fr * 100)}"]
    key = f"F{int(REF_FRAC * 100)}"
    cross["r_formal"] = cross["ddm"].abs() / cross["ddm_err"]
    cross["r_floorB"] = cross["ddm"].abs() / cross[f"sigB_{key}"]
    cross["r_floorA"] = cross["ddm"].abs() / cross[f"sigA_{key}"]
    for fr in FRACS:
        k = f"F{int(fr * 100)}"
        cross[f"rB_{k}"] = cross["ddm"].abs() / cross[f"sigB_{k}"]
        cross[f"rA_{k}"] = cross["ddm"].abs() / cross[f"sigA_{k}"]

    cont = set(cross.loc[cross.rep_nearest_arcsec1.lt(600) | cross.rep_nearest_arcsec2.lt(600),
                         "comp_id"])
    cross["contaminated"] = cross.comp_id.isin(cont)
    live = cross[~cross.contaminated].copy()

    cols = ["name1", "name2", "dt_days", "theta_arcsec", "ddm", "dm_mean",
            "s_F25", "nlocal_F25", "r_formal", "r_floorA", "r_floorB",
            "rB_F15", "rB_F25", "rB_F40", "rA_F15", "rA_F25", "rA_F40",
            "rep_nearest_arcsec1", "rep_nearest1", "flag_frac1", "flag_frac2",
            "bonsai_snr1", "bonsai_snr2"]
    out = live.sort_values("r_floorB")[cols].round(3).reset_index(drop=True)
    out.insert(0, "rank_floor", np.arange(1, len(out) + 1))
    out["rank_formal"] = out.r_formal.rank(method="min").astype(int)

    pd.set_option("display.width", 260)
    print("\n######## 存活新链接候选：ΔDM 一致度按经验地板重排（局部门窗 ±0.25·DM）########")
    print(out[["rank_floor", "name1", "name2", "dt_days", "theta_arcsec", "ddm", "dm_mean",
               "s_F25", "nlocal_F25", "r_formal", "r_floorB", "r_floorA", "rank_formal"]].to_string(index=False))

    print("\n######## 门窗敏感性（保守口径 r_floorB，±0.15/0.25/0.40·DM）########")
    print(out[["name1", "name2", "rB_F15", "rB_F25", "rB_F40",
               "rA_F15", "rA_F25", "rA_F40"]].to_string(index=False))

    ranks = {f"F{int(fr * 100)}":
             live.sort_values(f"rB_F{int(fr * 100)}")
             .apply(lambda r: f"{r.name1}x{r.name2}", axis=1).tolist()
             for fr in FRACS}
    stable = all(v == ranks["F25"] for v in ranks.values())
    verdict = "完全一致" if stable else "不一致 —— 需报告差异"
    print(f"\n排序在局部门窗 frac=0.15/0.25/0.40 下{verdict}")
    if not stable:
        for k, v in ranks.items():
            print(f"  {k}: {v}")

    top = out.iloc[0]
    second = out.iloc[1]
    print(f"\n== 领先幅度 ==")
    print(f"  #1 {top.name1}×{top.name2}: 形式 {top.r_formal:.2f} → 地板(保守) {top.r_floorB:.3f}"
          f" → 地板(宽松) {top.r_floorA:.3f}")
    print(f"  #2 {second.name1}×{second.name2}: 形式 {second.r_formal:.2f} → "
          f"地板(保守) {second.r_floorB:.3f}")
    print(f"  比值 #2/#1：形式口径 {second.r_formal / top.r_formal:.1f}× → "
          f"地板口径 {second.r_floorB / top.r_floorB:.1f}×  "
          f"（{'领先扩大' if second.r_floorB / top.r_floorB > second.r_formal / top.r_formal else '领先缩小'}）")

    drop = cross[cross.contaminated]
    if len(drop):
        print(f"\n剔除（F45 已判属已编目源）：{len(drop)} 对，"
              f"最近已编目源 = {sorted(set(drop.rep_nearest1.dropna()))}")

    OUTD.mkdir(parents=True, exist_ok=True)
    out.to_csv(OUTD / "c2_ranked_dm_sigma_v1.csv", index=False)
    out.to_parquet(OUTD / "c2_ranked_dm_sigma_v1.parquet", index=False)

    full = live.sort_values("r_floorB").reset_index(drop=True)
    full.insert(0, "rank_floor", np.arange(1, len(full) + 1))
    full["rank_formal"] = full.r_formal.rank(method="min").astype(int)
    full.to_parquet(OUTD / "c2_ranked_crossday_v14.parquet", index=False)
    full.to_csv(OUTD / "c2_ranked_crossday_v14.csv", index=False)
    if len(sameday):
        for fr in FRACS:
            k = f"F{int(fr * 100)}"
            sameday[f"s_{k}"] = sameday["dm_mean"].map(
                lambda x: local_scale(x, ref_dm, ref_D, fr)[0])
            sameday[f"rB_{k}"] = sameday["ddm"].abs() / (np.sqrt(2.0) * sameday[f"s_{k}"])
        sameday["r_formal"] = sameday["ddm"].abs() / sameday["ddm_err"]
        sameday.to_parquet(OUTD / "c2_ranked_sameday_v14.parquet", index=False)
    (OUTD / "c2_ranked_dm_sigma_manifest.json").write_text(json.dumps(dict(
        built_at=dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds"),
        logic_version=LOGIC_VERSION,
        change="ΔDM 一致度分母：形式误差 → Cat2 自测经验地板 s(DM)=1.4826·MAD(dm_fitb−bonsai_dm)，局部窗 ±0.25·DM",
        n_reference_events=int(m.sum()), local_fracs=list(FRACS), min_local_n=MIN_N,
        sigma_conventions=dict(conservative_sqrt2="sigma_comb=sqrt(2)*s",
                               loose="sigma_comb=s"),
        ranking_stable_across_windows=bool(stable),
        null_battery_affected=False,
        note="F49 置换显著性建在绝对 DM 窗上，本改动不触及；只影响以 σ 表述的 ΔDM 一致度",
        inputs=[RANKED.name, RAW.name]), ensure_ascii=False, indent=2), encoding="utf-8")

    md = ["# C2 候选 ΔDM 一致度（经验地板口径，logic v1.4.0）", "",
          f"分母 `s(DM)=1.4826·MAD(dm_fitb−bonsai_dm)`，局部门窗 ±0.25·DM，"
          f"参考集 n={int(m.sum())} 个 Cat2 事件。",
          "σ_comb 取 √2·s（保守）。同表给宽松口径。", "",
          "| # | 对 | Δt (d) | θ (″) | ΔDM | DM̄ | s(DM̄) | 旧 ΔDM/σ_形式 | **新（保守）** | 新（宽松） |",
          "|---|---|---|---|---|---|---|---|---|---|"]
    for _, r in out.iterrows():
        md.append(f"| {int(r['rank_floor'])} | {r['name1']}×{r['name2']} | {r['dt_days']:.2f} | "
                  f"{r['theta_arcsec']:.1f} | {r['ddm']:.3f} | {r['dm_mean']:.1f} | "
                  f"{r['s_F25']:.3f} | {r['r_formal']:.2f} | **{r['r_floorB']:.3f}** | "
                  f"{r['r_floorA']:.3f} |")
    md += ["", f"排序对局部门窗 frac（0.15/0.25/0.40）**{'稳定' if stable else '不稳定'}**。",
           "F49 置换显著性不受影响（其统计量建在绝对 DM 窗上）。"]
    (OUTD / "c2_ranked_dm_sigma_v1.md").write_text("\n".join(md) + "\n", encoding="utf-8")
    print(f"\n落盘 {OUTD.relative_to(ROOT)}/c2_ranked_dm_sigma_v1.{{csv,parquet,md}} + manifest")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
