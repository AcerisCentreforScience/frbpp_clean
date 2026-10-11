"""Tests for the upgrade side-channel: geometry, de-duplication, clustering, name census."""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from frbpp.upgrade.crossmatch import (            # noqa: E402
    EventTable, candidate_groups, cone_pairs, geometry_pairs, name_census,
    repeater_recall, same_event_mask)


def _table(ra, dec, mjd, dm, name, cat, rep=None, multi=None, ex=None):
    n = len(ra)
    return EventTable(np.array(ra, float), np.array(dec, float), np.array(mjd, float),
                      np.array(dm, float), np.array(name), np.array(cat),
                      np.array(rep if rep is not None else [""] * n),
                      np.array(multi if multi is not None else [False] * n),
                      np.array(ex if ex is not None else [False] * n))


def test_same_event_dedup():
    ev = _table([10, 10, 10, 10], [20, 20, 20, 20], [100, 100.01, 101, 100.0],
                [100, 100, 100, 100], ["A", "A", "A", "B"],
                ["cat2", "cat2", "cat2", "cat1"])
    i, j = geometry_pairs(ev, 60.0)
    m = same_event_mask(ev, i, j)
    kept = {(a, b) for a, b, k in zip(i, j, m) if not k}
    assert (0, 1) not in kept
    assert (0, 2) in kept and (1, 2) in kept
    assert (0, 3) not in kept


def test_cone_dm_clustering():
    ev = _table([10, 10.001, 50], [20, 20.0, 20], [100, 150, 200],
                [100, 100.5, 100], ["A", "B", "C"], ["cat2"] * 3)
    i, j = cone_pairs(ev, 60.0, 1.0)
    groups = candidate_groups(ev, i, j)
    assert len(groups) == 1 and set(groups[0].tolist()) == {0, 1}
    i2, j2 = cone_pairs(ev, 60.0, 0.2)
    assert i2.size == 0


def test_repeater_recall():
    ev = _table([10, 10.001, 50], [20, 20, 20], [100, 150, 200],
                [100, 100, 100], ["x1", "x2", "y"], ["cat1", "cat2", "cat2"],
                rep=["SRC", "SRC", "SRC"])
    i, j = cone_pairs(ev, 60.0, 5.0)
    known, rec = repeater_recall(ev, i, j)
    assert known == ["SRC"] and rec == ["SRC"]


def test_name_census_upgrade_detect():
    c1 = pd.DataFrame({"tns_name": ["N1", "N2"], "repeater_name": ["-9999", "-9999"],
                       "is_burst_primary": [True, True], "mjd_400": ["58000", "58010"]})
    c2 = pd.DataFrame({"tns_name": ["N1", "N2"], "mjd_400": ["58000.5", "58500"]})
    cen = name_census(c1, c2)
    names = [u["name"] for u in cen["upgrades"]]
    assert names == ["N2"] and cen["n_cat1_oneoff"] == 2