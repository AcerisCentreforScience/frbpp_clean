"""Fetch TNS per-object details for third-party attribution checks.

Usage: python scripts/tns_object_details.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from astropy import units as u
from astropy.coordinates import SkyCoord

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from tns_search_565d import TNS  # noqa: E402

OUT = ROOT / "data" / "processed" / "upgrade"

PAIR_MID = SkyCoord(20.185253 * u.deg, 51.148742 * u.deg)
STORM = ["FRB20230814G", "FRB20230825D", "FRB20230825E", "FRB20230825F",
         "FRB20230825G", "FRB20230825H", "FRB20230825I"]
PAIR = ["FRB20200511C", "FRB20211128E"]
SURV = ["FRB20211127E", "FRB20220623D", "FRB20200919A", "FRB20210111A",
        "FRB20200814C", "FRB20210829A", "FRB20211116C", "FRB20211129A",
        "FRB20220607C", "FRB20220621D", "FRB20200220I", "FRB20230721C"]
DMFIN = ["FRB20181017D", "FRB20200514D", "FRB20200620B", "FRB20200924A",
         "FRB20220616A", "FRB20230724A", "FRB20181017D"]
REF = ["FRB20220912A", "FRB20180916B"]

FIELDS = ["objid", "name_prefix", "objname", "object_type", "radeg", "decdeg",
          "radeg_err", "decdeg_err", "dm", "dm_err", "discoverydate", "discoverer",
          "discovery_data_source", "repeater_of", "internal_names", "hostname",
          "redshift", "reporter"]


def f(v):
    try:
        return float(v)
    except (TypeError, ValueError):
        return np.nan


def main() -> int:
    cfg = json.loads((ROOT / "config" / "tns.local.json").read_text(encoding="utf-8"))
    tns = TNS(cfg.get("tns_id", 0), cfg["bot_name"], cfg["api_key"], sleep=1.3)

    targets = []
    for grp, lst in [("pair", PAIR), ("storm_F45", STORM), ("surviving", SURV),
                     ("dm_fingerprint", DMFIN), ("reference", REF)]:
        for nm in lst:
            if nm not in [t[1] for t in targets]:
                targets.append((grp, nm))

    rows = []
    for grp, nm in targets:
        code, o = tns.get_object(nm)
        if code != 200 or not isinstance(o, dict):
            print(f"  [{grp}] {nm}: code={code} —— 无记录或失败")
            rows.append(dict(group=grp, query=nm, code=code, found=False))
            continue
        rep = o.get("repeater_of")
        repname = ""
        if isinstance(rep, dict):
            repname = f"{rep.get('name_prefix', '')}{rep.get('objname', '')}"
        elif isinstance(rep, list) and rep and isinstance(rep[0], dict):
            repname = f"{rep[0].get('name_prefix', '')}{rep[0].get('objname', '')}"
        elif rep not in (None, "", 0, "0"):
            repname = str(rep)
        ins = o.get("internal_names")
        if isinstance(ins, list):
            ins = "; ".join(str(x.get("internal_name") if isinstance(x, dict) else x)
                            for x in ins[:4])
        rows.append(dict(group=grp, query=nm, code=code, found=True,
                         objid=o.get("objid"), objtype=o.get("object_type"),
                         ra=f(o.get("radeg")), dec=f(o.get("decdeg")),
                         ra_err=f(o.get("radeg_err")), dec_err=f(o.get("decdeg_err")),
                         dm=f(o.get("dm")), dm_err=f(o.get("dm_err")),
                         disc_date=str(o.get("discoverydate"))[:10],
                         discoverer=str(o.get("discoverer"))[:24],
                         data_src=str(o.get("discovery_data_source"))[:14],
                         repeater_of=repname, internal=str(ins)[:60],
                         reporter=str(o.get("reporter"))[:20]))
        print(f"  [{grp}] {nm}: objid={o.get('objid')} DM={o.get('dm')} "
              f"repeater_of='{repname}' date={str(o.get('discoverydate'))[:10]}")

    df = pd.DataFrame(rows)
    for c in ("ra", "dec"):
        if c not in df:
            df[c] = np.nan
    df["ra"] = pd.to_numeric(df["ra"], errors="coerce")
    df["dec"] = pd.to_numeric(df["dec"], errors="coerce")
    df["sep_pair_arcsec"] = [
        PAIR_MID.separation(SkyCoord(r * u.deg, d * u.deg)).arcsec
        if np.isfinite(r) and np.isfinite(d) else np.nan
        for r, d in zip(df["ra"], df["dec"])]

    pd.set_option("display.width", 250)
    print("\n== 关键检验 1：F45 风暴簇在 TNS 里的 repeater_of ==")
    st = df[df.group == "storm_F45"]
    if "repeater_of" in st:
        print(st[["query", "dm", "sep_pair_arcsec", "repeater_of", "data_src",
                  "disc_date"]].round(1).to_string(index=False))
        linked = st[st.repeater_of.astype(str).str.len() > 0]
        print(f"  → 有 repeater_of 的 = {len(linked)}/{len(st)}"
              f"；指向 20220912A 的 = "
              f"{int(st.repeater_of.astype(str).str.contains('20220912').sum())}")

    print("\n== 关键检验 2：565 d 对在 TNS 里是否被记为重复源（F43 的第二口径）==")
    p = df[df.group == "pair"]
    print(p[["query", "objtype", "dm", "disc_date", "repeater_of", "discoverer",
             "data_src"]].to_string(index=False))

    print("\n== 关键检验 3：DM 指纹源与本对的角度（新候选链接筛查）==")
    dm = df[df.group == "dm_fingerprint"].sort_values("sep_pair_arcsec")
    print(dm[["query", "dm", "disc_date", "sep_pair_arcsec", "repeater_of"]]
          .round(1).to_string(index=False))
    close = dm[dm.sep_pair_arcsec < 3600]
    print(f"  → 角度 < 1° 的 DM 相似源 = {len(close)}")

    print("\n== 存活候选 12 事件的 TNS 归属 ==")
    sv = df[df.group == "surviving"]
    print(sv[["query", "dm", "disc_date", "sep_pair_arcsec", "repeater_of"]].round(1)
          .to_string(index=False))

    OUT.mkdir(parents=True, exist_ok=True)
    df.to_csv(OUT / "tns_object_details.csv", index=False)
    df.to_parquet(OUT / "tns_object_details.parquet", index=False)
    print(f"\n落盘 {OUT.relative_to(ROOT)}\\tns_object_details.{{csv,parquet}}  共 {len(df)} 行")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
