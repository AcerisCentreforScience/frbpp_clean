"""Exposure model and observation-log assembly."""
from __future__ import annotations

import numpy as np
from astropy import units as u
from astropy.time import Time

from frbpp.catalogs.chime import CHIME_SITE, mjd_to_bjd_tdb
from frbpp.contracts import ObservationLog

SIDEREAL_DAY_D = 0.99726958
_LST_RATE_PER_UTC_HOUR = 1.00273790935
_LOWER_TRANSIT_DEC_MIN = 90.0 - 49.320


def load_r3_exposure_series(npz_path):
    z = np.load(npz_path, allow_pickle=True)
    offsets = np.asarray(z["expo_unit_day_1s_step"], dtype=float)
    ts = z["tStart"]
    obj = ts.item() if getattr(ts, "shape", None) == () else ts
    recovered = "time_mjd"
    try:
        t_start = float(obj.mjd)
    except Exception:
        d = obj._time.__dict__
        t_start = float(d["jd1"]) + float(d["jd2"]) - 2400000.5
        recovered = "pickle_instance_dict"
    mjd = t_start + offsets
    meta = {
        "source": "CHIME/FRB Release 03 (R3_CHIME_exposure.npz)",
        "t_start_mjd_utc": t_start,
        "t_start_recovered_from": recovered,
        "n_samples": int(offsets.size),
        "sample_stride_median_s": float(np.median(np.diff(offsets)) * 86400.0),
    }
    return mjd, meta


def series_to_windows(mjd_times, max_gap_s: float = 120.0):
    t = np.asarray(mjd_times, dtype=float)
    if t.size == 0:
        return []
    breaks = np.where(np.diff(t) * 86400.0 > max_gap_s)[0]
    edges = np.concatenate(([0], breaks + 1, [t.size]))
    return [(float(t[i]), float(t[j - 1])) for i, j in zip(edges[:-1], edges[1:])]


def transit_grid(ra_deg: float, dec_deg: float, t0_mjd: float, t1_mjd: float,
                 location=None, n_iter: int = 5) -> np.ndarray:
    location = location if location is not None else CHIME_SITE
    n = int((t1_mjd - t0_mjd) / SIDEREAL_DAY_D) + 2
    t = t0_mjd + np.arange(n) * SIDEREAL_DAY_D
    ra_h = float(ra_deg) / 15.0
    for _ in range(n_iter):
        lt = Time(t, format="mjd", scale="utc", location=location)
        lst_h = lt.sidereal_time("apparent", longitude=location.lon).hour
        diff_h = (lst_h - ra_h + 12.0) % 24.0 - 12.0
        t = t - diff_h / (24.0 * _LST_RATE_PER_UTC_HOUR)
    m = (t >= t0_mjd) & (t <= t1_mjd)
    return np.sort(t[m])


def obslog_from_series(mjd_times, ra_deg: float, dec_deg: float,
                       pad_s: float = 2.0, meta_extra: dict | None = None):
    windows = series_to_windows(mjd_times)
    if not windows:
        raise ValueError("曝光序列为空, 无法组装 ObservationLog")
    t0 = np.array([w[0] for w in windows])
    t1 = np.array([w[1] for w in windows]) + pad_s / 86400.0
    log = ObservationLog(
        t_start=mjd_to_bjd_tdb(t0, ra_deg, dec_deg),
        t_end=mjd_to_bjd_tdb(t1, ra_deg, dec_deg),
        sens=np.ones(len(windows)),
        tgt_ra=float(ra_deg), tgt_dec=float(dec_deg),
        beam_ra=np.full(len(windows), np.nan),
        beam_dec=np.full(len(windows), np.nan),
        meta={"builder": "obslog_from_series", "pad_s": pad_s,
              "n_windows": len(windows), **(meta_extra or {})},
    )
    return log


def obslog_from_teff(ra_deg: float, dec_deg: float, t0_mjd: float, t1_mjd: float,
                     t_eff_s: float, window_hours: float = 0.2, location=None,
                     extra_meta: dict | None = None):
    location = location if location is not None else CHIME_SITE
    grids = [transit_grid(ra_deg, dec_deg, t0_mjd, t1_mjd, location=location)]
    if np.isfinite(dec_deg) and dec_deg > _LOWER_TRANSIT_DEC_MIN:
        grids.append(transit_grid((float(ra_deg) + 180.0) % 360.0, dec_deg,
                                  t0_mjd, t1_mjd, location=location))
    transit = np.sort(np.concatenate(grids))
    half_d = window_hours / 24.0 / 2.0
    w_lo = transit - half_d
    w_hi = transit + half_d
    span_d = float(np.sum(w_hi - w_lo))
    if span_d <= 0 or t_eff_s <= 0:
        raise ValueError(f"曝光不足：T_eff={t_eff_s:.1f}s, 窗口总长={span_d:.2f}d")
    sens_value = (t_eff_s / 86400.0) / span_d
    log = ObservationLog(
        t_start=mjd_to_bjd_tdb(w_lo, ra_deg, dec_deg),
        t_end=mjd_to_bjd_tdb(w_hi, ra_deg, dec_deg),
        sens=np.full(len(transit), sens_value),
        tgt_ra=float(ra_deg), tgt_dec=float(dec_deg),
        beam_ra=np.full(len(transit), np.nan),
        beam_dec=np.full(len(transit), np.nan),
        meta={"builder": "obslog_from_teff", "t_eff_hours": t_eff_s / 3600.0,
              "window_hours": window_hours, "sens": sens_value,
              **(extra_meta or {})},
    )
    return log, log.meta


def obslog_from_map(ra_deg: float, dec_deg: float, t0_mjd: float, t1_mjd: float,
                    exposure_L=None, exposure_U=None, window_hours: float = 0.2,
                    location=None):
    from astropy_healpix import HEALPix

    hpix = HEALPix(nside=4096, order="ring")
    ipix = int(hpix.lonlat_to_healpix(float(ra_deg) * u.deg, float(dec_deg) * u.deg))
    t_eff_s = 0.0
    used = []
    for name, arr in (("L", exposure_L), ("U", exposure_U)):
        if arr is not None and ipix < len(arr):
            t_eff_s += float(arr[ipix]) * 4.0
            used.append(name)
    return obslog_from_teff(ra_deg, dec_deg, t0_mjd, t1_mjd, t_eff_s,
                            window_hours=window_hours, location=location,
                            extra_meta={"builder": "obslog_from_map", "maps_used": used,
                                        "ipix": ipix})


class ExposureModel:

    def __init__(self, obslog: ObservationLog):
        self.log = obslog

    def epsilon(self, t_bjd, dm=None, snr=None):
        t = np.asarray(t_bjd, dtype=float)
        out = np.zeros(t.shape, dtype=float)
        for a, b, s in zip(self.log.t_start, self.log.t_end, self.log.sens):
            out[(t >= a) & (t < b)] = s
        return out

    def integral(self, lam_fn, t0, t1, dt_days: float = 60.0 / 86400.0) -> float:
        total = 0.0
        for a, b, s in zip(self.log.t_start, self.log.t_end, self.log.sens):
            lo, hi = max(float(a), float(t0)), min(float(b), float(t1))
            if hi <= lo:
                continue
            n = max(int(np.ceil((hi - lo) / dt_days)) + 1, 2)
            tt = np.linspace(lo, hi, n)
            lam = np.asarray(lam_fn(tt), dtype=float)
            total += float(np.sum(0.5 * (lam[1:] + lam[:-1]) * np.diff(tt)) * s)
        return total