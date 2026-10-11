"""Stage 1 acceptance gate against the raw catalogue; skipped when the data is absent."""
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw" / "chime"
_NEEDED = [
    "05_cat1/cat1_table2_vizier.tsv",
    "cat2_25/chimefrbcat2.csv",
    "03_180916/R3_CHIME_TOA_20Feb4.txt",
    "06_injections/chimefrb_catalog1_injections.p",
]
pytestmark = pytest.mark.skipif(
    not all((RAW / p).exists() for p in _NEEDED),
    reason="原始数据缺失：先运行 scripts/fetch_chime.py（见 docs/NOTES.md E7）",
)

sys.path.insert(0, str(ROOT / "src"))
from frbpp.catalogs import chime  # noqa: E402

TOA = RAW / "03_180916" / "R3_CHIME_TOA_20Feb4.txt"


def test_cat1_parse_count_and_time():
    df, cat = chime.load_cat1(RAW / "05_cat1" / "cat1_table2_vizier.tsv")
    assert len(df) == 600
    assert df["tns_name"].nunique() == 536
    assert cat.n_bursts == 536
    assert np.isfinite(cat.t).all()
    assert 2458000 < cat.t.min() and cat.t.max() < 2460000
    sub = df[df["is_burst_primary"]]
    off = cat.t - sub["mjd_400"].to_numpy(float)
    assert np.all(np.abs(off - 2400000.5) < 0.02)
    assert np.isnan(cat.dm_err).sum() == int(sub["dm_fitb_err"].isna().sum())


def test_cat2_parse_count():
    df, cat = chime.load_cat2(RAW / "cat2_25" / "chimefrbcat2.csv")
    assert len(df) == 5045
    assert cat.n_bursts == 5045
    assert np.isfinite(cat.t).all()
    assert cat.meta["n_ltt_fallback"] == 37


def test_g1_180916_identity():
    from astropy.time import Time

    df, _ = chime.load_r3_180916_toa(TOA)
    mjd = df["mjd_native"].to_numpy()
    assert len(mjd) == 38, "论文: 38 bursts"
    assert abs(mjd.min() - 58377.429721) < 1e-6
    assert abs(mjd.max() - 58883.0552353) < 1e-6
    assert abs((mjd.max() - mjd.min()) - 505.6255143) < 1e-3
    assert Time(mjd.min(), format="mjd").utc.iso[:10] == "2018-09-16"
    assert Time(mjd.max(), format="mjd").utc.iso[:10] == "2020-02-04"


def test_g1_injections_fields_and_monotonic():
    inj = chime.load_injections(RAW / "06_injections" / "chimefrb_catalog1_injections.p")
    need = {"beam_x", "beam_y", "bonsai_snr", "dm_det", "dm_inj", "fluence_jy_ms_inj",
            "l2_rfi_grade", "t_utc_det", "t_utc_inj", "pulse_width_ms_inj",
            "spectral_index_inj"}
    assert need.issubset(set(inj.columns))
    assert len(inj) > 80_000

    det = (inj["bonsai_snr"].notna() & (inj["bonsai_snr"] > 9)).to_numpy()
    flu = inj["fluence_jy_ms_inj"].to_numpy(float)
    ok = np.isfinite(flu)
    qs = np.quantile(flu[ok], np.linspace(0, 1, 7))
    fracs = []
    for lo, hi in zip(qs[:-1], qs[1:]):
        m = ok & (flu >= lo) & (flu <= hi)
        fracs.append(float(det[m].mean()))
    fracs = np.asarray(fracs)
    assert fracs[-1] > fracs[0]
    assert np.all(np.diff(fracs) > -0.02)
    from scipy.stats import spearmanr
    rho, _ = spearmanr(np.arange(len(fracs)), fracs)
    assert rho > 0.9, f"效率曲线非单调: fracs={fracs}"