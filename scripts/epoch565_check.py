"""Recheck Cat2 coverage around the pair epochs within a 2-degree cone and DM window.

Usage: python scripts/epoch565_check.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from astropy.time import Time

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
CAT2 = ROOT / "data" / "interim" / "cat2_bursts.parquet"

NAME1, NAME2 = "FRB20200511C", "FRB20211128E"
MJD1, MJD2 = 58980.748549148826, 59546.20229763866
DT = MJD2 - MJD1


def mjd(s):
    return float(Time(s, scale="utc").mjd)


cat = pd.read_parquet(CAT2)
r1 = cat[cat.tns_name == NAME1].iloc[0]
r2 = cat[cat.tns_name == NAME2].iloc[0]
ra = 0.5 * (float(r1.ra) + float(r2.ra))
dec = 0.5 * (float(r1.dec) + float(r2.dec))
dm0 = 0.5 * (float(r1.dm_fitb) + float(r2.dm_fitb))
print(f"目标: RA={ra:.5f} Dec={dec:.5f} DM={dm0:.2f} | Cat2 span "
      f"{mjd('2018-07-25'):.1f}..{mjd('2023-09-15'):.1f}")

print("\n== 565.45 d 族历元（t0=第一暴）==")
for k in range(0, 7):
    e = MJD1 + k * DT
    print(f"  k={k}: MJD {e:9.2f}  UTC {Time(e, format='mjd', scale='utc').iso[:10]}")

T = mjd("2023-06-16")
print(f"\n== Cat2 在 2023-06-16 ±30 d（MJD {T:.1f} ±30）==")
w = cat[(cat.mjd_400 >= T - 30) & (cat.mjd_400 <= T + 30)]
print(f"  该窗口 Cat2 总事件数 = {len(w)}")
if len(w):
    print(w[["tns_name", "mjd_400", "ra", "dec", "dm_fitb", "snr_fitb"]]
          .sort_values("mjd_400").to_string(index=False))

d2r = np.radians
c = np.cos(d2r(dec))
sep = np.degrees(np.arccos(np.clip(
    np.sin(d2r(cat.dec.to_numpy(float))) * np.sin(d2r(dec))
    + np.cos(d2r(cat.dec.to_numpy(float))) * np.cos(d2r(dec))
    * np.cos(d2r(cat.ra.to_numpy(float) - ra)), -1, 1)))
near = cat[(sep < 2.0) & np.isfinite(sep)]
print(f"\n== 全 Cat2：2° 锥内事件 = {len(near)} ==")
if len(near):
    out = near[["tns_name", "mjd_400", "ra", "dec", "dm_fitb"]].copy()
    out["sep_arcsec"] = (sep[near.index.to_numpy()] * 3600).round(1)
    out["dDM"] = (out.dm_fitb - dm0).round(2)
    print(out.sort_values("mjd_400").to_string(index=False))

fp = cat[(cat.dm_fitb >= dm0 - 5) & (cat.dm_fitb <= dm0 + 5)]
print(f"\n== DM 指纹 {dm0:.0f}±5 全天空 = {len(fp)} 条 ==")
if len(fp):
    o = fp[["tns_name", "mjd_400", "ra", "dec", "dm_fitb"]].copy()
    o["sep_arcsec"] = (sep[fp.index.to_numpy()] * 3600).round(1)
    print(o.sort_values("mjd_400").to_string(index=False))

echo = cat[(sep < 0.1) & np.isfinite(sep)]
print(f"\n== 6′ 内（≈该源定位尺度）事件 = {len(echo)} ==")
if len(echo):
    o = echo[["tns_name", "mjd_400", "ra", "dec", "dm_fitb", "snr_fitb"]].copy()
    o["sep_arcsec"] = (sep[echo.index.to_numpy()] * 3600).round(1)
    print(o.sort_values("mjd_400").to_string(index=False))


import json  # noqa: E402

RAWD = ROOT / "data" / "raw"


def _sep(ras, decs):
    c_ = np.cos(np.radians(dec))
    return np.degrees(np.arccos(np.clip(
        np.sin(np.radians(decs)) * np.sin(np.radians(dec))
        + np.cos(np.radians(decs)) * np.cos(np.radians(dec))
        * np.cos(np.radians(np.asarray(ras) - ra)), -1, 1))) * 3600.0


print("\n===== 跨目录交叉（该位置 / 该历元）=====")

ap = json.loads((RAWD / "apertif" / "apertif_frbs.json").read_text(encoding="utf-8"))
ar = ap["rows"]
ra_a = np.array([r["RAJ2000"] for r in ar], float)
de_a = np.array([r["DEJ2000"] for r in ar], float)
mj_a = np.array([r["MJD"] for r in ar], float)
sp_a = _sep(ra_a, de_a)
print(f"[Apertif] n={len(ar)}  MJD {mj_a.min():.1f}→{mj_a.max():.1f} "
      f"({Time(mj_a.min(), format='mjd').iso[:10]} → {Time(mj_a.max(), format='mjd').iso[:10]})")
print(f"  2° 锥内 = {int((sp_a < 7200).sum())}；最近 3 个:")
for i in np.argsort(sp_a)[:3]:
    print(f"    {ar[i]['TNS']:14s} sep={sp_a[i]:9.1f}\"  MJD={mj_a[i]:.1f}  "
          f"DM={ar[i]['DM']}  z={ar[i].get('z')}")
print(f"  覆盖 2023-06-16 的 Apertif 事件数 = {int(((mj_a > T - 30) & (mj_a < T + 30)).sum())}")

fa = json.loads((RAWD / "fast" / "fast_20220912A_bursts.json").read_text(encoding="utf-8"))
print(f"[FAST] 表={fa['table']} n={fa['n_rows']} → 仅 FRB20220912A 单源，"
      f"**对本位置无覆盖**（无坐标列）")

gcn = (RAWD / "gcn" / "history_20260905_20261008.jsonl")
print(f"[GCN] 留存窗 2026-09-05→10-08（仅 ~1 月）→ **结构上无法覆盖 2023-06-16**")

print("[Blinkverse] 仅源级表（无逐暴）→ 无法查具体事件（F36）")
print("[TNS] 未接入（bot 待审核）→ **这是唯一能覆盖非 CHIME 探测的官方渠道，目前空缺**")


RAWCAT2 = ROOT / "data" / "raw" / "chime" / "cat2_25" / "chimefrbcat2.csv"
raw = pd.read_csv(RAWCAT2, low_memory=False)
print(f"\n===== Cat2 原始表：{len(raw)} 行，旗标普查 =====")
for c in ["excluded_flag", "sidelobe_flag", "intrachan_flag",
          "citizen_science_flag", "catalog1_flag", "catalog1_param_flag"]:
    if c in raw.columns:
        print(f"  {c}: {dict(raw[c].astype(str).str.strip().value_counts().head(4))}")
snr = pd.to_numeric(raw.bonsai_snr, errors="coerce")
print(f"  bonsai_snr: min={np.nanmin(snr):.2f} p1={np.nanquantile(snr, 0.01):.2f} "
      f"median={np.nanmedian(snr):.2f} | SNR<10 的暴 = {int((snr < 10).sum())}")

sep_raw = _sep(raw.ra.to_numpy(float), raw.dec.to_numpy(float))
for R in (7200.0, 10800.0):
    n = raw[sep_raw < R]
    print(f"\n  位置 {R/3600:.0f}° 锥内（含被剔行）= {len(n)} 条")
    if len(n):
        cols = [c for c in ["tns_name", "mjd_400", "dm_fitb", "bonsai_snr", "snr_fitb",
                            "excluded_flag", "sidelobe_flag", "intrachan_flag",
                            "flag_frac", "catalog1_flag", "notes_fitb"] if c in raw.columns]
        o = n[cols].copy()
        o["sep_arcsec"] = sep_raw[n.index.to_numpy()].round(1)
        print(o.sort_values("mjd_400").to_string(index=False))


