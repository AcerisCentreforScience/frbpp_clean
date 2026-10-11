"""Tests for the exposure model and observation-log assembly."""
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from frbpp.contracts import ObservationLog  # noqa: E402
from frbpp.exposure import model as em  # noqa: E402

RAW_R3 = ROOT / "data/raw/chime/03_180916/R3_CHIME_exposure.npz"
RAW_TOA = ROOT / "data/raw/chime/03_180916/R3_CHIME_TOA_20Feb4.txt"
RAW_CAT1 = ROOT / "data/raw/chime/05_cat1/cat1_table2_vizier.tsv"
RAW_LMAP = (ROOT / "data/raw/chime/05_cat1/"
            "exposure_int_20180828_20190702_transit_L_beam_FWHM-600_res_4s_0.86_arcmin.npz")


def test_exposure_model_integral_window_guard():
    log = ObservationLog(
        t_start=np.array([100.0, 102.0]), t_end=np.array([101.0, 103.0]),
        sens=np.array([0.5, 0.5]), tgt_ra=10.0, tgt_dec=20.0,
        beam_ra=np.array([np.nan, np.nan]), beam_dec=np.array([np.nan, np.nan]))
    m = em.ExposureModel(log)
    lam = lambda t: np.full(np.shape(t), 2.0)  # noqa: E731
    assert np.isclose(m.integral(lam, 99.0, 105.0), 2.0)
    assert np.isclose(m.integral(lam, 101.2, 101.8), 0.0)
    assert m.epsilon(np.array([100.5, 101.5, 102.5])).tolist() == [0.5, 0.0, 0.5]


def test_transit_grid_sidereal_spacing_and_lst():
    from astropy.time import Time

    t = em.transit_grid(29.503, 65.717, 58377.0, 58387.0)
    assert 9 <= len(t) <= 11
    assert np.all(np.abs(np.diff(t) - em.SIDEREAL_DAY_D) < 1e-3)
    lt = Time(t, format="mjd", scale="utc", location=em.CHIME_SITE)
    lst = lt.sidereal_time("apparent", longitude=em.CHIME_SITE.lon).hour
    diff = (lst - 29.503 / 15.0 + 12.0) % 24.0 - 12.0
    assert np.all(np.abs(diff) < 1e-3)


@pytest.mark.skipif(not RAW_R3.exists(), reason="R3 曝光数据缺失")
def test_r3_series_windows_and_coverage():
    from frbpp.catalogs import chime

    mjd, meta = em.load_r3_exposure_series(RAW_R3)
    assert mjd.size == 106477
    assert abs(meta["t_start_mjd_utc"] - 58358.480118) < 5e-4
    wins = em.series_to_windows(mjd)
    assert 280 <= len(wins) <= 300
    total_h = sum(b - a for a, b in wins) * 24.0
    assert 65.0 < total_h < 76.0
    est_sample_h = mjd.size * meta["sample_stride_median_s"] / 3600.0
    assert 58.0 < est_sample_h < 61.0
    df, _ = chime.load_r3_180916_toa(RAW_TOA)
    toa = df["mjd_native"].to_numpy()
    covered = sum(any(a <= x <= b for a, b in wins) for x in toa)
    assert int(covered) == 17


@pytest.mark.skipif(not RAW_R3.exists(), reason="R3 曝光数据缺失")
def test_obslog_from_series_contract():
    mjd, _ = em.load_r3_exposure_series(RAW_R3)
    log = em.obslog_from_series(mjd, 29.503, 65.717)
    n = len(log.t_start)
    assert n == log.sens.size == log.beam_ra.size == log.beam_dec.size
    assert np.all(log.sens == 1.0)
    assert np.all(log.t_end > log.t_start)
    assert 2458300 < log.t_start.min() and log.t_start.max() < 2458800
    assert np.all(np.isnan(log.beam_ra))


@pytest.mark.skipif(not (RAW_LMAP.exists() and RAW_CAT1.exists()),
                    reason="L 曝光图或 Cat1 数据缺失")
def test_l_map_vs_catalog_anchor():
    from astropy import units as u
    from astropy_healpix import HEALPix

    from frbpp.catalogs import chime

    L = np.asarray(np.load(RAW_LMAP, allow_pickle=True)["exposure"], dtype=float)
    df, _ = chime.load_cat1(RAW_CAT1)
    sub = df[df["is_burst_primary"]]
    hpix = HEALPix(nside=4096, order="ring")
    ip = hpix.lonlat_to_healpix(sub["ra"].to_numpy(dtype=float) * u.deg,
                                sub["dec"].to_numpy(dtype=float) * u.deg)
    assert np.all(ip < len(L))
    val_h = L[ip] * 4.0 / 3600.0
    dec = sub["dec"].to_numpy(dtype=float)
    low = sub["exp_low"].to_numpy(dtype=float)
    m = np.isfinite(low) & (dec > 41.0)
    ratio = val_h[m] / low[m]
    med = float(np.median(ratio))
    assert 0.7 < med < 1.8, f"L×exp_low 比值中位 {med:.2f} 偏离锚点 1.13"