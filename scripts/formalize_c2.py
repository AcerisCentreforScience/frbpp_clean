"""Freeze the C2 candidate set from the G4b grid.

Usage: python scripts/formalize_c2.py
Output: data/processed/upgrade/c2_candidates_v1.{parquet,csv}, c2_formalization.md
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from frbpp.upgrade.crossmatch import clean_rep                       # noqa: E402

RAW2 = ROOT / "data" / "interim" / "cat2_bursts.parquet"
REPORT = ROOT / "data" / "processed" / "upgrade" / "g4b_report.json"
OUT = ROOT / "data" / "processed" / "upgrade"

_SEARCHED = {"FRB20200511C", "FRB20211128E", "FRB20200919A", "FRB20210111A",
             "FRB20230814G", "FRB20230825D", "FRB20220607C", "FRB20220621D"}


def pub_status(n1, n2):
    if n1 in _SEARCHED or n2 in _SEARCHED:
        return "通用网络检索无结果（2026-10-08）；待 Blinkverse/TNS 权威核查（v2）"
    return "未外部核查（v1）"


def _f(x):
    try:
        v = float(x)
        return v if np.isfinite(v) else np.nan
    except Exception:  # noqa: BLE001
        return np.nan


def haversine_arcsec(ra1, dec1, ra2, dec2):
    if not all(np.isfinite([ra1, dec1, ra2, dec2])):
        return np.nan
    r1, d1, r2, d2 = map(np.radians, [ra1, dec1, ra2, dec2])
    x = np.clip(np.sin(d1) * np.sin(d2) + np.cos(d1) * np.cos(d2) * np.cos(r1 - r2), -1, 1)
    return float(np.degrees(np.arccos(x)) * 3600.0)


def main() -> int:
    rep = json.loads(REPORT.read_text(encoding="utf-8"))
    c = rep["candidates"]
    cat2 = pd.read_parquet(RAW2)
    cat2 = cat2.assign(rep_clean=clean_rep(cat2.repeater_name))

    prim = {}
    for tns, g in cat2.groupby("tns_name"):
        g0 = g.sort_values("sub_num")
        r0 = g0.iloc[0]
        prim[str(tns)] = {
            "ra": _f(r0["ra"]), "dec": _f(r0["dec"]),
            "dm": _f(r0["dm_fitb"]), "dm_err": _f(r0["dm_fitb_err"]),
            "snr": _f(r0["snr_fitb"]), "fluence": _f(r0["fluence"]),
            "width": _f(r0["width_fitb"]), "mjd": _f(r0["mjd_400"]),
            "rep": str(r0["rep_clean"]) if r0["rep_clean"] not in (None, "") else "",
            "c1": int(_f(r0["catalog1_flag"]) == 1) if np.isfinite(_f(r0["catalog1_flag"])) else -1,
            "nsub": int(len(g0)),
        }

    CRITERION_DM_WINDOW = 5.0
    rows, missing, edge = [], [], []
    for grp, arr in (("same_day", c["same_day"]), ("cross_day", c["cross_day"])):
        for p in arr:
            a, b = prim.get(p["name1"]), prim.get(p["name2"])
            if a is None or b is None:
                missing.append((p["name1"], p["name2"]))
                continue
            if float(p["dm_window"]) > CRITERION_DM_WINDOW:
                edge.append({"group": grp, "name1": p["name1"], "name2": p["name2"],
                             "dt_days": round(float(p["dt_days"]), 4),
                             "theta_arcsec": p["theta"], "dm_window": float(p["dm_window"]),
                             "ddm": round(a["dm"] - b["dm"], 4)})
                continue
            ddm = a["dm"] - b["dm"]
            e1, e2 = a["dm_err"], b["dm_err"]
            ddm_err = float(np.sqrt(e1**2 + e2**2)) if np.isfinite(e1) and np.isfinite(e2) else np.nan
            theta = haversine_arcsec(a["ra"], a["dec"], b["ra"], b["dec"])
            grouped = bool(a["rep"] and b["rep"] and a["rep"] == b["rep"])
            ungrouped = bool(not a["rep"] and not b["rep"])
            if grp == "cross_day":
                cls = ("跨日强候选" if (p["dt_days"] >= 10.0 and ungrouped)
                       else ("跨日候选（兼容并档）" if grouped else "跨日候选"))
            else:
                cls = "同日档（风暴/双探测）"
            rows.append({
                "group": grp, "cls": cls,
                "name1": p["name1"], "name2": p["name2"],
                "dt_days": round(float(p["dt_days"]), 4),
                "theta_cfg_arcsec": p["theta"], "theta_arcsec": (round(theta, 1)
                                                                 if np.isfinite(theta) else None),
                "dm1": round(a["dm"], 4), "dm2": round(b["dm"], 4),
                "ddm": round(float(ddm), 4),
                "ddm_err": (round(ddm_err, 4) if np.isfinite(ddm_err) else None),
                "ddm_sigma": (round(abs(ddm) / ddm_err, 2)
                              if np.isfinite(ddm_err) and ddm_err > 0 else None),
                "snr1": a["snr"], "snr2": b["snr"],
                "fluence1": a["fluence"], "fluence2": b["fluence"],
                "nsub1": a["nsub"], "nsub2": b["nsub"],
                "official_grouped": grouped,
                "official_ungrouped": ungrouped,
                "catalog1_flag_1": a["c1"], "catalog1_flag_2": b["c1"],
                "rep1": a["rep"], "rep2": b["rep"],
                "ra1": a["ra"], "dec1": a["dec"], "ra2": b["ra"], "dec2": b["dec"],
                "public_status": pub_status(p["name1"], p["name2"]),
            })

    df = pd.DataFrame(rows).sort_values(["group", "dt_days"], ascending=[True, True])
    df.to_parquet(OUT / "c2_candidates_v1.parquet", index=False)
    df.to_csv(OUT / "c2_candidates_v1.csv", index=False)
    df_edge = pd.DataFrame(edge).sort_values(["group", "dt_days"], ascending=[True, True])
    if len(df_edge):
        df_edge.to_csv(OUT / "c2_candidates_edge_v1.csv", index=False)

    n_grouped = int(df["official_grouped"].sum())
    n_ung = int(df["official_ungrouped"].sum())
    strong = df[df["cls"] == "跨日强候选"]
    print(f"正式化 {len(df)} 对（缺 {len(missing)}）| 官方已归组 {n_grouped} | 官方未归组 {n_ung}")
    print(f"跨日强候选 {len(strong)} 对：")
    for r in strong.itertuples():
        s = f"{r.ddm_sigma:.1f}σ" if r.ddm_sigma is not None else "err-NA"
        print(f"   {r.name1} × {r.name2}  Δt={r.dt_days:.1f}d ΔDM={r.ddm:+.2f} ({s}) "
              f"θ={r.theta_arcsec}\"")

    lines = ["# C2 候选正式化 v1（G4b 尾巴闭环）", "",
             f"- 输入：`g4b_report.json` 九档并集 {len(df) + len(df_edge)} 对"
             f"（同日 {int((df['group'] == 'same_day').sum())} + 跨日 "
             f"{int((df['group'] == 'cross_day').sum())} + 边缘 "
             f"{len(df_edge)}）；本表取 DM 窗 <= {CRITERION_DM_WINDOW:.0f} 的判据池 {len(df)} 对",
             f"- 仅由更松档（DM 窗 > {CRITERION_DM_WINDOW:.0f}）支持的边缘对另存 "
             f"`c2_candidates_edge_v1.csv`；null bar 仍按九档并集重算",
             f"- 官方分组对照（Cat2 `repeater_name` 双侧同非空）：已归组 **{n_grouped}**、"
             f"未归组 **{n_ung}**",
             f"- 分类：同日档 {int((df['group']=='same_day').sum())} / 跨日候选 "
             f"{int((df['group']=='cross_day').sum())}（其中强 {len(strong)}）", "",
             "## 跨日候选（按 Δt 降序）", "",
             "| # | 对 | Δt (d) | ΔDM | ΔDM/σ | θ (″) | SNR | 官方归组 | 类别 |",
             "|---|---|---|---|---|---|---|---|---|"]
    for i, r in enumerate(df[df["group"] == "cross_day"].itertuples(), 1):
        s = f"{r.ddm_sigma:.1f}σ" if r.ddm_sigma is not None else "—"
        snr = f"{r.snr1:.0f}/{r.snr2:.0f}"
        lines.append(f"| {i} | {r.name1} × {r.name2} | {r.dt_days:.1f} | {r.ddm:+.2f} | {s} "
                     f"| {r.theta_arcsec} | {snr} | {'是' if r.official_grouped else '否'} "
                     f"| {r.cls} |")
    lines += ["", "> 注：20230814G × 20230825D–I 六配对同属一个物理假设（单一源的风暴群成员）；"
              "θ=600″ 档为搜索配置边缘。官方分组对照 = Cat2 `repeater_name` 双侧同非空。", "",
              "## 同日档（风暴/双探测，15 对）", "",
              "| # | 对 | Δt (h) | ΔDM | θ (″) | 官方归组 |",
              "|---|---|---|---|---|---|"]
    for i, r in enumerate(df[df["group"] == "same_day"].itertuples(), 1):
        lines.append(f"| {i} | {r.name1} × {r.name2} | {r.dt_days*24:.1f} | {r.ddm:+.2f} "
                     f"| {r.theta_arcsec} | {'是' if r.official_grouped else '否'} |")
    lines += ["", "## 下一步（v2 协议）", "",
              "1. 跨日强候选：外部目录复查（Blinkverse/TNS/社区公告）——public_status 列待填",
              "2. 候选升级源：加入 C3 实时核验的监控清单（新事件即自动升级判定）",
              "3. 观测侧：v2 对接 FAST/MeerKAT 申请流程（超出本仓库范围，登记）"]
    (OUT / "c2_formalization.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"→ c2_candidates_v1.parquet/csv + c2_formalization.md")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())