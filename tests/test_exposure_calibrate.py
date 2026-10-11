"""Tests for injection-recovery calibration."""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from frbpp.exposure import calibrate  # noqa: E402

RAW_INJ = ROOT / "data/raw/chime/06_injections/chimefrb_catalog1_injections.p"


def test_wilson_interval_bounds():
    lo, hi = calibrate.wilson_interval(np.array([50.0]), np.array([100.0]))
    assert 0.0 <= lo[0] < 0.5 < hi[0] <= 1.0
    lo2, hi2 = calibrate.wilson_interval(np.array([0.0]), np.array([0.0]))
    assert (lo2[0], hi2[0]) == (0.0, 0.0)
    lo3, hi3 = calibrate.wilson_interval(np.array([1.0]), np.array([1.0]))
    assert lo3[0] <= 1.0 <= hi3[0] + 1e-12


def test_efficiency_curve_synthetic():
    flu = np.array([1.0] * 10 + [10.0] * 10 + [100.0] * 10)
    det = np.array([False] * 8 + [True] * 2 + [False] * 5 + [True] * 5 + [True] * 10)
    inj = pd.DataFrame({
        "fluence_jy_ms_inj": flu,
        "bonsai_snr": np.where(det, 20.0, np.nan),
        "l2_rfi_grade": np.where(det, 8.0, np.nan),
    })
    cur = calibrate.efficiency_curve(inj, "fluence_jy_ms_inj",
                                     bins=[0.5, 5.0, 50.0, 200.0])
    assert len(cur) == 3
    assert np.allclose(cur["eff"], [0.2, 0.5, 1.0])
    assert np.all(cur["eff_lo"] <= cur["eff"] + 1e-12)
    assert np.all(cur["eff"] <= cur["eff_hi"] + 1e-12)


@pytest.mark.skipif(not RAW_INJ.exists(), reason="注入原始数据缺失")
def test_injections_efficiency_monotone_fluence():
    from scipy.stats import spearmanr

    from frbpp.catalogs import chime

    inj = chime.load_injections(RAW_INJ)
    cur = calibrate.efficiency_curve(inj, "fluence_jy_ms_inj", bins=8)
    assert cur["n"].sum() > 50_000
    rho, _ = spearmanr(cur["bin_lo"], cur["eff"])
    assert rho > 0.9
    assert cur["eff"].iloc[-1] > cur["eff"].iloc[0]
    assert 0.02 < float(cur["eff"].min()) and float(cur["eff"].max()) <= 1.0