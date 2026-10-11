"""Unit tests for the data contracts: length checks, NaN handling, routing invariants, coordinate identity."""
import numpy as np
import pytest

from frbpp.contracts import (
    BurstCatalog,
    SourceRecord,
    angular_separation_deg,
    enters_periodicity_likelihood,
    match_coordinates,
    pipeline_for,
)


def _kwargs(n=3, **overrides):
    kw = dict(
        t=np.arange(n, dtype=float) + 2459000.0,
        source=np.array(["FRB20180916B"] * n),
        dm=np.full(n, 348.8),
        dm_err=np.full(n, np.nan),
        snr=np.full(n, 10.0),
        fluence=np.full(n, np.nan),
        band_lo=np.full(n, 400.0),
        band_hi=np.full(n, 800.0),
        ra=np.full(n, np.nan),
        dec=np.full(n, np.nan),
        flags=np.zeros(n, dtype=int),
        meta={"telescope": "CHIME", "catalog_version": "unit-test"},
    )
    kw.update(overrides)
    return kw


def test_burst_catalog_length_check():
    with pytest.raises(ValueError):
        BurstCatalog(**_kwargs(dm=np.zeros(2)))


def test_missing_values_stay_nan():
    cat = BurstCatalog(**_kwargs())
    assert np.all(np.isnan(cat.dm_err))
    assert np.all(np.isnan(cat.fluence))
    assert cat.n_bursts == 3
    assert cat.meta["telescope"] == "CHIME"


def test_source_record_category_validation():
    rec = SourceRecord(name="FRB20180916B", ra=29.5, dec=65.7,
                       category="repeater", n_bursts=40, n_episodes=12)
    assert rec.min_episodes == 5
    with pytest.raises(ValueError):
        SourceRecord(name="X", ra=0.0, dec=0.0, category="unknown",
                     n_bursts=1, n_episodes=1)


def test_layered_routing_invariants():
    assert not enters_periodicity_likelihood("one_off")
    assert not enters_periodicity_likelihood("suspected")
    assert enters_periodicity_likelihood("repeater")
    assert "periodicity" not in pipeline_for("one_off")
    assert "periodicity" not in pipeline_for("suspected")
    assert "periodicity" in pipeline_for("repeater")
    assert pipeline_for("one_off") == ("upgrade_side_channel", "population_baseline")


def test_coordinate_identity_10arcsec():
    ra, dec = 29.5, 65.7
    assert match_coordinates(ra, dec, ra, dec + 5e-5)
    assert not match_coordinates(ra, dec, ra, dec + 1.0 / 60.0)
    sep = float(angular_separation_deg(ra, dec, ra, dec + 10.0 / 3600.0))
    assert abs(sep - 10.0 / 3600.0) < 1e-9