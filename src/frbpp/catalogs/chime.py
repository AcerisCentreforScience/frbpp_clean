"""CHIME/FRB Catalogue 1, Catalogue 2 and Injections readers returning contract-validated frames."""
from __future__ import annotations

import io
from pathlib import Path

import numpy as np
import pandas as pd
from astropy import units as u
from astropy.coordinates import EarthLocation, SkyCoord
from astropy.time import Time

from frbpp.contracts import FLAG_EXCLUDED, BurstCatalog

CHIME_SITE = EarthLocation.from_geodetic(-119.624 * u.deg, 49.320 * u.deg, 545.0 * u.m)

VIZIER_CAT1_MAP = {
    "Name": "tns_name", "OName": "previous_name", "RpName": "repeater_name",
    "RAJ2000": "ra", "e_RAJ2000": "ra_err",
    "DEJ2000": "dec", "e_DEJ2000": "dec_err",
    "GLON": "gl", "GLAT": "gb",
    "ExpUp": "exp_up", "e_ExpUp": "exp_up_err",
    "ExpLo": "exp_low", "e_ExpLo": "exp_low_err",
    "SNR": "bonsai_snr", "DM": "bonsai_dm",
    "SNRfitb": "snr_fitb", "DMfitb": "dm_fitb", "e_DMfitb": "dm_fitb_err",
    "DMeNE2001": "dm_exc_ne2001", "DMeYMW16": "dm_exc_ymw16",
    "BcWidth": "bc_width", "Scat": "scat_time", "e_Scat": "scat_time_err",
    "Flux": "flux", "e_Flux": "flux_err",
    "Fluence": "fluence", "e_Fluence": "fluence_err",
    "Nsb": "sub_num",
    "MJD400": "mjd_400", "e_MJD400": "mjd_400_err",
    "MJDinf": "mjd_inf", "e_MJDinf": "mjd_inf_err",
    "Widthfitb": "width_fitb", "e_Widthfitb": "width_fitb_err",
    "SpInd": "sp_idx", "e_SpInd": "sp_idx_err",
    "spRun": "sp_run", "e_spRun": "sp_run_err",
    "B_Freq": "high_freq", "b_Freq": "low_freq", "Fpk": "peak_freq",
    "chi2": "chi_sq", "dof": "dof", "Frac": "flag_frac", "Flag": "excluded_flag",
}


def mjd_to_bjd_tdb(mjd, ra_deg, dec_deg, location=None) -> np.ndarray:
    t, _ = mjd_to_bjd_tdb_masked(mjd, ra_deg, dec_deg, location=location)
    return t


def mjd_to_bjd_tdb_masked(mjd, ra_deg, dec_deg, location=None):
    location = location if location is not None else CHIME_SITE
    t = Time(np.asarray(mjd, dtype=float), format="mjd", scale="utc",
             location=location)
    ra = np.asarray(ra_deg, dtype=float)
    dec = np.asarray(dec_deg, dtype=float)
    known = np.isfinite(ra) & np.isfinite(dec)
    lst_deg = t.sidereal_time("apparent", longitude=location.lon).to_value(u.deg)
    ra_use = np.where(known, ra, lst_deg)
    dec_use = np.where(known, dec, location.lat.to_value(u.deg))
    coord = SkyCoord(ra_use * u.deg, dec_use * u.deg)
    ltt = t.light_travel_time(coord, kind="barycentric", location=location)
    return np.asarray((t.tdb + ltt).jd, dtype=float), known


def parse_vizier_tsv(path) -> tuple[pd.DataFrame, dict[str, str]]:
    lines = Path(path).read_text(encoding="utf-8", errors="replace").splitlines()
    header, rows = None, []
    for ln in lines:
        if not ln.strip() or ln.startswith("#"):
            continue
        if set(ln.replace("\t", "").strip()) == {"-"}:
            continue
        if header is None:
            header = [h.strip() for h in ln.split("\t")]
            continue
        rows.append(ln)
    if header is None or not rows:
        raise ValueError(f"VizieR TSV 解析失败: {path}")
    df = pd.read_csv(io.StringIO("\n".join(["\t".join(header), *rows])), sep="\t",
                     na_values=["", " ", "-9999", "null", "NaN"])
    for c in df.columns:
        if pd.api.types.is_string_dtype(df[c]) or df[c].dtype == object:
            s = df[c].str.strip()
            df[c] = s.mask(s == "")
    df = df.rename(columns={c: VIZIER_CAT1_MAP.get(c, c) for c in df.columns})
    return df, VIZIER_CAT1_MAP


def _require(df: pd.DataFrame, cols, table: str) -> None:
    missing = [c for c in cols if c not in df.columns]
    if missing:
        raise ValueError(f"{table} 缺少必需列: {missing}; 实际列: {list(df.columns)}")


def _flags_from_excluded(df: pd.DataFrame) -> np.ndarray:
    exc = df["excluded_flag"].to_numpy(dtype=float)
    flags = np.zeros(len(df), dtype=int)
    flags[np.nan_to_num(exc, nan=0.0).astype(bool)] = FLAG_EXCLUDED
    return flags


def _build_catalog(df: pd.DataFrame, meta: dict) -> BurstCatalog:
    n = len(df)
    band_lo = df["low_freq"].to_numpy(float) if "low_freq" in df else np.full(n, np.nan)
    band_hi = df["high_freq"].to_numpy(float) if "high_freq" in df else np.full(n, np.nan)
    return BurstCatalog(
        t=df["t_bjd_tdb"].to_numpy(float),
        source=df["tns_name"].astype(str).to_numpy(),
        dm=df["dm_fitb"].to_numpy(float),
        dm_err=df["dm_fitb_err"].to_numpy(float),
        snr=df["bonsai_snr"].to_numpy(float),
        fluence=df["fluence"].to_numpy(float),
        band_lo=band_lo,
        band_hi=band_hi,
        ra=df["ra"].to_numpy(float),
        dec=df["dec"].to_numpy(float),
        flags=_flags_from_excluded(df),
        meta=meta,
    )


_CAT_REQUIRED = ("tns_name", "ra", "dec", "bonsai_snr", "dm_fitb", "dm_fitb_err",
                 "fluence", "mjd_400", "excluded_flag")


def load_cat1(tsv_path) -> tuple[pd.DataFrame, BurstCatalog]:
    df, _ = parse_vizier_tsv(tsv_path)
    n_raw = len(df)
    df = df[pd.to_numeric(df["recno"], errors="coerce").notna()].copy()
    if "sub_num" not in df.columns:
        raise ValueError("Cat1 缺少 sub_num 列（期望 VizieR Nsb）")
    _require(df, _CAT_REQUIRED, "Cat1")
    df["t_bjd_tdb"] = mjd_to_bjd_tdb(df["mjd_400"], df["ra"], df["dec"])
    df["is_burst_primary"] = pd.to_numeric(df["sub_num"], errors="coerce") == 0
    burst = df[df["is_burst_primary"]].copy()
    meta = {
        "telescope": "CHIME/FRB",
        "catalog": "Catalog 1 (VizieR J/ApJS/257/59/table2; doi:10.26093/cds/vizier.22570059)",
        "native_time": "mjd_400 = MJD(UTC, topocentric) arrival at 400 MHz",
        "time_system": "BJD_TDB (barycentric light-travel applied; DRAO site)",
        "level": "burst-level (sub_num==0 of the sub-burst table)",
        "n_rows_raw": int(n_raw), "n_rows_sub": int(len(df)),
        "source_file": str(tsv_path),
    }
    return df, _build_catalog(burst, meta)


def load_cat2(csv_path) -> tuple[pd.DataFrame, BurstCatalog]:
    df = pd.read_csv(csv_path)
    _require(df, _CAT_REQUIRED, "Cat2")
    df = df.copy()
    t, known = mjd_to_bjd_tdb_masked(df["mjd_400"], df["ra"], df["dec"])
    df["t_bjd_tdb"] = t
    df["ltt_position_known"] = known
    meta = {
        "telescope": "CHIME/FRB",
        "catalog": "Catalog 2 (CANFAR doi:10.11570 via CISTI.CANFAR/25.0066)",
        "native_time": "mjd_400 = MJD(UTC, topocentric) arrival at 400 MHz",
        "time_system": "BJD_TDB (barycentric light-travel applied; DRAO site)",
        "n_ltt_fallback": int((~known).sum()),
        "ltt_fallback": "无坐标行用本地天顶方向近似(~秒级); 仅日尺度统计可用",
        "source_file": str(csv_path),
    }
    return df, _build_catalog(df, meta)


R3_RA_DEG, R3_DEC_DEG = 29.503, 65.717


def load_r3_180916_toa(path) -> tuple[pd.DataFrame, dict]:
    vals = []
    for ln in Path(path).read_text(encoding="utf-8").splitlines():
        ln = ln.strip()
        if not ln or ln.startswith("#"):
            continue
        vals.append(float(ln))
    mjd = np.asarray(vals, dtype=float)
    df = pd.DataFrame({"mjd_native": mjd,
                       "t_bjd_tdb": mjd_to_bjd_tdb(mjd, R3_RA_DEG, R3_DEC_DEG)})
    meta = {
        "source": "CHIME/FRB Release 03 (CANFAR CISTI.CANFAR/20.0002; R3_CHIME_TOA_20Feb4.txt)",
        "native_time": "MJD (UTC, topocentric)",
        "time_system": "BJD_TDB (barycentric light-travel applied; coords = FRB 20180916B)",
        "radec_used": (R3_RA_DEG, R3_DEC_DEG),
    }
    return df, meta


def load_injections(pickle_path) -> pd.DataFrame:
    return pd.read_pickle(pickle_path)