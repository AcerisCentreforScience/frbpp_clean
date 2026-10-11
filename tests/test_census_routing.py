"""Tests for stratified routing, sub-burst folding, hygiene fields."""
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from frbpp.census.routing import build_event_frame, route_sources   # noqa: E402


def _cat2(rows):
    return pd.DataFrame(rows, columns=["tns_name", "repeater_name", "sub_num",
                                       "mjd_400", "ra", "dec", "dm_fitb"])


def test_routing_categories():
    df = _cat2([
        ("FRB20180001A", "R1", 0, "58000", "10", "20", "100"),
        ("FRB20180002A", "R1", 0, "58002", "10", "20", "100"),
        ("FRB20180003A", "R1", 0, "58004", "10", "20", "100"),
        ("FRB20180004A", "R1", 0, "58007", "10", "20", "100"),
        ("FRB20180005A", "R1", 0, "58011", "10", "20", "100"),
        ("FRB20180006A", "R2", 0, "58100", "30", "40", "200"),
        ("FRB20180007A", "R2", 0, "58112", "30", "40", "200"),
        ("FRB20180008A", "R2", 0, "58130", "30", "40", "200"),
        ("FRB20180009A", None, 0, "58200", "50", "60", "300"),
        ("FRB20180010A", None, 0, "58300", None, None, "400"),
    ])
    ev = build_event_frame(df)
    tab, eps = route_sources(ev)
    cats = dict(zip(tab.src, tab.category))
    assert cats["R1"] == "A" and len(eps["R1"]) == 5
    assert cats["R2"] == "B" and len(eps["R2"]) == 3
    assert cats["FRB20180009A"] == "C"
    assert cats["FRB20180010A"] == "D"


def test_sub_burst_folding_and_hygiene():
    df = _cat2([
        ("FRB20180001A", "R1", 0, "58000.0000001", "10", "20", "100"),
        ("FRB20180001A", "R1", 1, "58000.0000002", "10", "20", "100"),
        ("FRB20180002A", "R1", 0, "58002.0", "10", "20", "100"),
    ])
    ev = build_event_frame(df)
    assert len(ev[ev.src == "R1"]) == 2
    tab, _ = route_sources(ev)
    r1 = tab[tab.src == "R1"].iloc[0]
    assert r1.n_bursts == 2 and r1.n_episodes == 2 and r1.category == "B"
    assert r1.rate1_per_day > 0 and r1.rate2_per_day >= 0