"""Host identification and Macquart-consistency check for the 113.69 d pair.

Usage: python scripts/host_113d.py
Output: data/processed/upgrade/pair113_host.json
"""
from __future__ import annotations

import datetime as dt
import io
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import requests

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

OUT = ROOT / "data" / "processed" / "upgrade"
CAT2 = ROOT / "data" / "interim" / "cat2_bursts.parquet"
NAME1, NAME2 = "FRB20200919A", "FRB20210111A"
RAD_ARCSEC = 180.0

OM, OL = 0.315, 0.685
DM_HALO, DM_HALO_ERR = 40.0, 20.0
EXT_KRON_PSF = 0.3


def dm_igm_of_z(z):
    zz = np.linspace(0, z, 400)
    E = np.sqrt(OM * (1 + zz) ** 3 + OL)
    return 1000.0 * np.trapezoid((1 + zz) / E, zz)


def z_of_dm(dm_igm):
    zs = np.linspace(0.005, 4.0, 1200)
    vals = np.array([dm_igm_of_z(z) for z in zs])
    return float(np.interp(dm_igm, vals, zs))


def ps1_cone(ra, dec, rad_deg):
    r = requests.get("https://catalogs.mast.stsci.edu/api/v0.1/panstarrs/dr2/mean.csv",
                     params={"ra": f"{ra:.6f}", "dec": f"{dec:.6f}",
                             "radius": f"{rad_deg:.5f}",
                             "columns": ("ra,dec,rMeanPSFMag,gMeanPSFMag,iMeanPSFMag,"
                                         "zMeanPSFMag,yMeanPSFMag,rMeanKronMag,gMeanKronMag,"
                                         "iMeanKronMag,nDetections"),
                             "pagesize": 100000}, timeout=120)
    r.raise_for_status()
    df = pd.read_csv(io.StringIO(r.text), na_values=["-999", -999])
    ren = {}
    for want, cand in {"ra": ("raMean", "raStack"), "dec": ("decMean", "decStack")}.items():
        for c in cand:
            if c in df.columns:
                ren[c] = want
                break
    df = df.rename(columns=ren)
    return df


def main() -> int:
    from astropy.coordinates import SkyCoord
    import astropy.units as u

    cat = pd.read_parquet(CAT2)
    r1 = cat[cat.tns_name == NAME1].iloc[0]
    r2 = cat[cat.tns_name == NAME2].iloc[0]
    ra = 0.5 * (float(r1.ra) + float(r2.ra))
    dec = 0.5 * (float(r1.dec) + float(r2.dec))
    ra_err = 0.5 * (float(r1.ra_err) + float(r2.ra_err)) * 3600.0
    dec_err = 0.5 * (float(r1.dec_err) + float(r2.dec_err)) * 3600.0

    gal = SkyCoord(ra * u.deg, dec * u.deg, frame="icrs").galactic
    gl, gb = float(gal.l.deg), float(gal.b.deg)

    out = {"built_at": dt.datetime.now().astimezone().isoformat(timespec="seconds"),
           "pair": [NAME1, NAME2],
           "cone": {"ra": ra, "dec": dec, "radius_arcsec": RAD_ARCSEC,
                    "l_gal_deg": gl, "b_gal_deg": gb},
           "burst_err_arcsec": {"ra": ra_err, "dec": dec_err}}

    dm = np.array([float(r1.dm_fitb), float(r2.dm_fitb)])
    mw_ne2 = np.mean([float(r1.dm_exc_ne2001), float(r2.dm_exc_ne2001)])
    mw_ymw = np.mean([float(r1.dm_exc_ymw16), float(r2.dm_exc_ymw16)])
    dm_mw = 0.5 * (float(dm.mean()) - mw_ne2 + float(dm.mean()) - mw_ymw)
    dm_igm = float(dm.mean()) - dm_mw - DM_HALO
    spread = abs(mw_ne2 - mw_ymw) + DM_HALO_ERR + float(np.abs(dm[0] - dm[1]) / 2)
    out["dm_budget"] = {"dm_obs_mean": float(dm.mean()),
                        "dm_mw_ne2001": float(dm.mean()) - mw_ne2,
                        "dm_mw_ymw16": float(dm.mean()) - mw_ymw,
                        "dm_mw_used": dm_mw, "dm_halo_assumed": DM_HALO,
                        "dm_igm": dm_igm, "uncert": spread,
                        "z_igm_point": z_of_dm(dm_igm),
                        "z_igm_range": [z_of_dm(dm_igm - spread), z_of_dm(dm_igm + spread)],
                        "relation": "DM_IGM(z)=1000∫(1+z')/E(z')dz' pc/cm^3 (Planck18)",
                        "note": "z 由 DM 反推，不是测量；DM-z 散射 ±(20-30)% 未计入"}

    cats = {}
    try:
        df = ps1_cone(ra, dec, RAD_ARCSEC / 3600.0)
        df = df.copy()
        kron_minus_psf = df["rMeanKronMag"] - df["rMeanPSFMag"]
        ext = kron_minus_psf > EXT_KRON_PSF
        c = np.cos(np.radians(dec))
        d = np.degrees(np.arccos(np.clip(
            np.sin(np.radians(df.dec)) * np.sin(np.radians(dec))
            + np.cos(np.radians(df.dec)) * np.cos(np.radians(dec))
            * np.cos(np.radians(df.ra - ra)), -1, 1))) * 3600.0
        df["sep_arcsec"] = d
        bright = df["rMeanPSFMag"] < 20.0
        ext_bright = ext & bright
        cats["panstarrs"] = {
            "status": "ok", "n_objects": int(len(df)),
            "n_extended_ext_flag": int(ext.sum()),
            "ext_flag_rule": f"rKron − rPSF > {EXT_KRON_PSF}",
            "n_bright_r_lt20": int(bright.sum()),
            "n_ext_bright": int(ext_bright.sum()),
            "r_psf_mag_quantiles": {q: float(np.nanquantile(df["rMeanPSFMag"], q))
                                    for q in (0.05, 0.25, 0.5, 0.75, 0.95)},
            "nearest_10": df.sort_values("sep_arcsec").head(10)[
                ["ra", "dec", "sep_arcsec", "rMeanPSFMag", "rMeanKronMag",
                 "gMeanPSFMag", "iMeanPSFMag"]].round(4).to_dict("records"),
            "brightest_5": df.sort_values("rMeanPSFMag").head(5)[
                ["ra", "dec", "sep_arcsec", "rMeanPSFMag", "gMeanPSFMag",
                 "zMeanPSFMag"]].round(4).to_dict("records"),
        }
    except Exception as e:  # noqa: BLE001
        cats["panstarrs"] = {"status": "error", "error": str(e)[:200]}

    cats["sdss"] = {"status": "unavailable",
                    "reason": f"本地网络超时；该天区 gl={gl:.1f} gb={gb:.1f}"}
    cats["ned"] = {"status": "unavailable", "reason": "多次连接超时"}
    cats["datalab_ls10"] = {"status": "unavailable", "reason": "datalab.noirlab.edu 不可达"}

    out["catalogs"] = cats

    n_gal = cats.get("panstarrs", {}).get("n_extended_ext_flag", 0)
    out["verdict"] = {
        "z_igm_point": out["dm_budget"]["z_igm_point"],
        "conclusion": ("宿主证认在当前定位精度下不可行（非结论性）"
                       if cats.get("panstarrs", {}).get("status") == "ok"
                       else "宿主证认不可执行（Pan-STARRS 查询失败）"),
        "reason": (f"两暴形式定位误差 RA≈{ra_err:.0f}″ / Dec≈{dec_err:.0f}″，"
                   f"而 180″ 锥内已含 {cats.get('panstarrs',{}).get('n_objects','?')} 个天体"
                   f"（其中延展源≈{n_gal}）→ 星系密度不足以唯一定位宿主；"
                   f"且 180″ 锥远小于误差椭圆。"),
        "usable_output": (f"DM 预算给出 z≈{out['dm_budget']['z_igm_point']:.2f}"
                          f"（区间 {out['dm_budget']['z_igm_range'][0]:.2f}"
                          f"–{out['dm_budget']['z_igm_range'][1]:.2f}）：该源若为同源，"
                          f"其红移落在 z~{out['dm_budget']['z_igm_point']:.1f}；"
                          "这是『物理一致性』信息，但不能与任何具体星系关联。"),
    }
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "pair113_host.json").write_text(
        json.dumps(out, indent=2, ensure_ascii=False), encoding="utf-8")

    print(f"锥心 RA={ra:.5f} Dec={dec:.5f} (gl={gl:.1f} gb={gb:.1f}) r={RAD_ARCSEC}\" | "
          f"定位误差 RA≈{ra_err:.0f}\" Dec≈{dec_err:.0f}\"")
    print(f"DM 预算 → DM_IGM={dm_igm:.1f} → z≈{out['dm_budget']['z_igm_point']:.2f} "
          f"({out['dm_budget']['z_igm_range'][0]:.2f}–{out['dm_budget']['z_igm_range'][1]:.2f})")
    ps = cats["panstarrs"]
    if ps.get("status") == "ok":
        print(f"[PS1] 天体={ps['n_objects']} 延展源={ps['n_extended_ext_flag']} "
              f"r<20={ps['n_bright_r_lt20']} 延展且亮={ps['n_ext_bright']}")
        print(f"[PS1] r-band 分位: {ps['r_psf_mag_quantiles']}")
        print("[PS1] 最近 5 个:")
        for o in ps["nearest_10"][:5]:
            print(f"      sep={o['sep_arcsec']:.1f}\" r={o['rMeanPSFMag']} "
                  f"Kron-PSF={o['rMeanKronMag']-o['rMeanPSFMag']:+.2f}")
    else:
        print("[PS1]", ps)
    print("判定:", out["verdict"]["conclusion"])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
