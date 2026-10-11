#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Assemble the reproduction package and assert every headline count.

Usage: python scripts/repro_package.py
Output: reproduction/{README.md,tables/,json/,logs/,MANIFEST.csv}
"""
from __future__ import annotations

import hashlib
import json
import shutil
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))
sys.stdout.reconfigure(encoding="utf-8", errors="replace")

from frbpp.upgrade.crossmatch import build_event_table  # noqa: E402
from dm_sigma_floor import RAW, local_scale  # noqa: E402

UP = ROOT / "data" / "processed" / "upgrade"
OUT = ROOT / "reproduction"
PAIR = ("FRB20200511C", "FRB20211128E")
CONTAM_SEED = "FRB20230814G"
WIN_FRAC = 0.25
PRIMARY_FLOOR_RATIO = 0.119

README = """# Reproduction package

Derived data for "An isolated candidate same-source FRB pair hidden in the apparent one-off
population". Nothing here needs private data: every input is a public catalogue, and every
artefact is produced by a script in `scripts/` with a fixed seed.

## Frozen conventions (do not drift)

| Item | Value used for all headline numbers |
| --- | --- |
| Event table | v1.3 (2026-10-09) event-layer de-duplication: 5608 catalogue rows minus 594 repeated Catalogue 1 rows = **5014 distinct events** |
| Apparent one-off subspace | 3742 events |
| Declaration domain | **3665** events with finite DM and coordinates |
| Candidate set | **28** frozen criterion pairs (30-pair union is used only as the stricter null bar) |
| Independent hypotheses | 8 union-find components, **7 surviving** after the audit removes one |
| DM error floor | s(DM) = 1.4826 MAD(D), D = dm_fitb - bonsai_dm, local window +-0.25 DM, n = 4961 |
| Two-burst combination | sigma_comb = sqrt(2) s (conservative); the loose s is reported alongside |
| Nulls | B = 1000 per model; p is a resolution limit p <= 1/(B+1), never a measured value |
| Catalogue release | one file, `data/raw/chime/cat2_25/chimefrbcat2.csv`; floor statistics must be recomputed on any other release |

## Contents and generators

| File | Produced by |
| --- | --- |
| `tables/events_frozen.csv` | `scripts/frbpp/upgrade/crossmatch.build_event_table` via `scripts/repro_package.py` |
| `tables/candidates_28pairs.csv` | `scripts/formalize_c2.py` (pairs) + union-find and `scripts/dm_sigma_floor.local_scale` here |
| `tables/dm_floor_bins.csv` | `scripts/dm_floor_robust.py` |
| `json/dm_floor_robust.json` | `scripts/dm_floor_robust.py` (bootstrap CIs, three windows, rank stability) |
| `json/same_beam_audit.json` | `scripts/same_beam_test.py` (conditional pools, rho sweep, hour angles) |
| `json/isolation_background.json` | `scripts/build_isolation_data.py` (nearest-neighbour ECDF, cross-day curve) |
| `json/null_battery_summary.json` | `scripts/c2_null_battery.py` |
| `json/event_table_audit.json`, `json/known_repeater_control.json` | `scripts/event_table_dedup_audit.py`, `scripts/known_repeater_control.py` |
| `json/mc_settings.json` | draws and seeds of every Monte Carlo used in the text |
| `logs/tns_search_ledger.jsonl`, `logs/tns_search_summary.md` | `scripts/tns_search_565d.py` (read-only TNS Search API) |

`candidates_28pairs.csv` carries, per pair: component id, same/cross-day flag, measured theta,
absolute ddm, sigma_floor, ddm/sigma_floor, ddm/sigma_formal, kept or removed, and the reason.
The removed component is the seven-event group audited as a catalogue grouping omission of
FRB20220912A; the 21 removed pair entries are retained here for audit.

## How to recompute the manuscript's numbers

```
python scripts/event_table_dedup_audit.py      # Table 1 counts: 5608 -> 5014 -> 3742 -> 3665
python scripts/formalize_c2.py                 # the 28 frozen candidate pairs
python scripts/c2_null_battery.py              # Table: four permutation models, B = 1000
python scripts/dm_sigma_floor.py               # ranked pairs under the empirical floor
python scripts/dm_floor_robust.py              # DM-floor appendix: bins, bootstrap, windows
python scripts/same_beam_test.py --draws 400 --seed 20261011   # same-beam audit appendix
python scripts/build_isolation_data.py         # Figure 8 background
python scripts/repro_package.py                # rebuild this package (asserts every count)
```

Figures: `python scripts/make_figs.py --only all --dpi 900 --outdir docs/figures/hires`
(joint isolation, sidereal quantization, scattering discriminant, null battery) and
`python scripts/make_layout_figs2.py` (nearest-neighbour ECDF and counting experiment,
flow, roadmap). Figure numbers in the manuscript are assigned by the class at compile time;
the filenames above are the stable identifiers. Design width is the
measured text block, 15.00 cm, with a 7 pt font floor; each generator re-runs a collision report
before saving.

## Seeds

Null battery and robustness pipeline: 20261009. Error-convolved same-beam audit and DM-floor
bootstrap: 20261011. Isolation background: recorded in `json/isolation_background.json`.
"""


def union_find_components(pairs: list[tuple[str, str]]) -> dict[str, str]:
    parent: dict[str, str] = {}

    def find(x: str) -> str:
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for a, b in pairs:
        parent.setdefault(a, a)
        parent.setdefault(b, b)
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[rb] = ra
    return {n: find(n) for n in parent}


def build_events() -> pd.DataFrame:
    c1 = pd.read_parquet(ROOT / "data/interim/cat1_bursts.parquet")
    c2 = pd.read_parquet(ROOT / "data/interim/cat2_bursts.parquet")
    ev = build_event_table(c1, c2, event_dedup=True)
    raw = pd.read_csv(RAW, low_memory=False)
    sub = raw[(raw.sub_num == 0) | raw.sub_num.isna()].copy()
    sub["tns_name"] = sub["tns_name"].astype("string")
    keep_cols = ["tns_name", "bonsai_dm", "dm_fitb", "dm_fitb_err", "scat_time",
                 "excluded_flag", "sidelobe_flag", "citizen_science_flag", "intrachan_flag",
                 "catalog1_flag", "catalog1_param_flag", "fluence_win_extended_flag",
                 "mjd_400", "t_bjd_tdb", "repeater_name"]
    keep_cols = [c for c in keep_cols if c in sub.columns]
    first = sub.dropna(subset=["tns_name"]).drop_duplicates("tns_name", keep="first")
    d = pd.DataFrame({
        "name": np.asarray(ev.name),
        "mjd": np.asarray(ev.mjd, float),
        "ra_deg": np.asarray(ev.ra, float),
        "dec_deg": np.asarray(ev.dec, float),
        "dm_used": np.asarray(ev.dm, float),
        "src_rep": [str(s) for s in ev.src_rep],
    })
    d = d.merge(first[keep_cols], left_on="name", right_on="tns_name", how="left")
    for c in ("ra_err", "dec_err"):
        if c in first.columns:
            d[c] = d["name"].map(first.set_index("tns_name")[c])
    d["chime_beam_id"] = "not in catalogue"
    d["zenith_angle_deg"] = "not in catalogue"
    d["in_declared_domain"] = ((d["src_rep"] == "") & np.isfinite(d["dm_used"])
                               & np.isfinite(d["ra_deg"])).map({True: 1, False: 0})
    d["is_repeater_linked"] = (d["src_rep"] != "").map({True: 1, False: 0})
    return d.drop(columns=["tns_name"], errors="ignore")


def build_candidates(events: pd.DataFrame) -> pd.DataFrame:
    c = pd.read_csv(UP / "c2_candidates_v1.csv")
    assert len(c) == 28, f"冻结候选应为 28 对，实得 {len(c)}"
    comps = union_find_components(list(zip(c.name1, c.name2)))
    c["comp_id"] = [comps[n] for n in c.name1]
    order = {v: i for i, v in enumerate(dict.fromkeys(c.comp_id))}
    c["component_id"] = c.comp_id.map(order)
    c = c.rename(columns={"theta_arcsec": "theta_arcsec_measured",
                          "theta_cfg_arcsec": "theta_arcsec_config",
                          "ddm": "ddm_abs"})
    out = c[["component_id", "group", "name1", "name2", "dt_days",
             "theta_arcsec_measured", "ddm_abs", "dm1", "dm2"]].copy()
    out["same_or_cross_day"] = np.where(c["group"].astype(str).str.startswith("cross"),
                                         "cross-day", "same-day")
    out["dm_mean"] = (out.dm1 + out.dm2) / 2.0
    raw = pd.read_csv(RAW, low_memory=False)
    b = pd.to_numeric(raw.bonsai_dm, errors="coerce").to_numpy(float)
    f = pd.to_numeric(raw.dm_fitb, errors="coerce").to_numpy(float)
    ref_D = f - b
    m = np.isfinite(f) & np.isfinite(ref_D)
    ref_dm, ref_D = f[m], ref_D[m]

    def floor_ratio(dm_mean: float, ddm: float) -> tuple[float, float, int]:
        s, used, nloc = local_scale(dm_mean, ref_dm, ref_D, WIN_FRAC)
        sig = float("nan") if not np.isfinite(s) else np.sqrt(2.0) * s
        r = float("nan") if not np.isfinite(sig) or sig <= 0 else ddm / sig
        return r, sig, int(nloc)

    fr = [floor_ratio(float(x), float(y)) for x, y in zip(out.dm_mean, out.ddm_abs)]
    out["sigma_floor_comb_pc_cm3"] = [round(x[1], 4) for x in fr]
    out["ddm_over_sigma_floor"] = [round(x[0], 4) for x in fr]
    out["n_local_floor_window_F25"] = [x[2] for x in fr]
    if "ddm_sigma" in c:
        out["ddm_over_sigma_formal"] = c["ddm_sigma"].round(4).values
    contam_id = int(c.loc[c.name1 == CONTAM_SEED, "component_id"].iloc[0])
    out["kept_for_claims"] = np.where(out.component_id == contam_id, "removed", "kept")
    out["reason"] = np.where(
        out.kept_for_claims == "removed",
        "catalogue grouping omission of FRB20220912A: component lies 1.6-93.6 arcsec from that "
        "named source's per-burst cloud (Section 3.4 audit); retained for audit, excluded from "
        "new-link claims",
        "surviving new-link hypothesis (one independent physical hypothesis per component)")
    return out.drop(columns=["dm1", "dm2"])


def summarize_null() -> dict:
    j = json.loads((UP / "c2_null_battery.json").read_text(encoding="utf-8"))
    rows = {}
    for k, v in (j.get("perm") or {}).items():
        if not isinstance(v, dict):
            continue
        rows[k] = {kk: v.get(kk) for kk in
                   ("name", "n_reps", "cand_mean", "cand_p95", "cand_max", "real_cand",
                    "n_ge_real", "single_p", "global_p_bonferroni") if kk in v}
    return {"config": j.get("config", {}), "models": rows,
            "reference_bar": "30-pair union over the five significant grid configurations "
                             "(stricter than the 28-pair set claimed)",
            "single_test_p": "p <= 1/(B+1) = 0.001 at B=1000 (resolution limit, not a measurement)",
            "global_p": "p <= 0.009 after Bonferroni x 9 grid configurations",
            "known_gap": "the battery's own JSON does not embed its seed; see mc_settings.json"}


def main() -> int:
    (OUT / "tables").mkdir(parents=True, exist_ok=True)
    (OUT / "json").mkdir(parents=True, exist_ok=True)
    (OUT / "logs").mkdir(parents=True, exist_ok=True)

    ev = build_events()
    assert len(ev) == 5014, f"冻结事件表应为 5014 行，实得 {len(ev)}"
    assert int(ev.in_declared_domain.sum()) == 3665, "声明域应为 3665"
    assert int(ev.is_repeater_linked.sum()) == 1272, "重复源关联事件应为 1272"
    ev.to_csv(OUT / "tables/events_frozen.csv", index=False)

    cand = build_candidates(ev)
    assert len(cand) == 28
    assert cand.component_id.nunique() == 8, f"分量应为 8 个，实得 {cand.component_id.nunique()}"
    assert (cand.kept_for_claims == "removed").sum() == 21, "被审计掉的应是 21 对"
    assert (cand.kept_for_claims == "kept").sum() == 7, "留下的应是 7 对"
    prim = cand[((cand.name1 == PAIR[0]) & (cand.name2 == PAIR[1]))
                | ((cand.name1 == PAIR[1]) & (cand.name2 == PAIR[0]))]
    assert len(prim) == 1, "主候选在 28 对表里应恰好出现一次"
    got = float(prim.ddm_over_sigma_floor.iloc[0])
    assert abs(got - PRIMARY_FLOOR_RATIO) < 5e-3, \
        f"主候选地板一致度 {got:.4f} 与 c2_dm_floor_robust_v1.json 的 {PRIMARY_FLOOR_RATIO} 不符"
    cand.to_csv(OUT / "tables/candidates_28pairs.csv", index=False)

    shutil.copyfile(UP / "c2_dm_floor_bins_v1.csv", OUT / "tables/dm_floor_bins.csv")
    for src, dst in ((UP / "c2_dm_floor_robust_v1.json", "dm_floor_robust.json"),
                     (UP / "same_beam_test.json", "same_beam_audit.json"),
                     (UP / "c2_isolation_figure_v1.json", "isolation_background.json"),
                     (UP / "event_table_audit.json", "event_table_audit.json"),
                     (UP / "known_repeater_control.json", "known_repeater_control.json")):
        shutil.copyfile(src, OUT / "json" / dst)
    (OUT / "json/null_battery_summary.json").write_text(
        json.dumps(summarize_null(), ensure_ascii=False, indent=2), encoding="utf-8")
    sb = json.loads((UP / "same_beam_test.json").read_text(encoding="utf-8"))
    df = json.loads((UP / "c2_dm_floor_robust_v1.json").read_text(encoding="utf-8"))
    (OUT / "json/mc_settings.json").write_text(json.dumps({
        "same_beam_error_convolution": {"draws": sb["error_convolved_counts"]["draws"],
                                        "seed": sb["error_convolved_counts"]["seed"]},
        "dm_floor_bootstrap": {"draws": df["window_frac_0.25"]["boot_draws"],
                               "seed": df["window_frac_0.25"]["seed"]},
        "null_battery": {"B_per_model": 1000, "seed_constant": 20261009,
                         "rerun": "python scripts/c2_null_battery.py --seed 20261009",
                         "seed_recorded_in": "the --seed default of argparse in "
                                             "scripts/c2_null_battery.py",
                         "known_gap": "the script writes its products without the seed, so the "
                                      "seed is recovered from the source, not from the JSON"},
        "generators": ["scripts/formalize_c2.py", "scripts/c2_null_battery.py",
                       "scripts/dm_sigma_floor.py", "scripts/dm_floor_robust.py",
                       "scripts/same_beam_test.py", "scripts/build_isolation_data.py",
                       "scripts/event_table_dedup_audit.py", "scripts/tns_search_565d.py"]},
        ensure_ascii=False, indent=2), encoding="utf-8")
    shutil.copyfile(UP / "tns_search_ledger.jsonl", OUT / "logs/tns_search_ledger.jsonl")
    shutil.copyfile(UP / "tns_search_summary.md", OUT / "logs/tns_search_summary.md")

    (OUT / "README.md").write_text(README, encoding="utf-8", newline="\n")

    rows = []
    for p in sorted(OUT.rglob("*")):
        if p.is_file() and p.name != "MANIFEST.csv":
            rows.append({"path": str(p.relative_to(OUT)), "bytes": p.stat().st_size,
                         "sha256": hashlib.sha256(p.read_bytes()).hexdigest()})
    pd.DataFrame(rows).to_csv(OUT / "MANIFEST.csv", index=False)
    print(f"复现包 {len(rows)} 个文件 -> {OUT.relative_to(ROOT)}")
    print(f"  events_frozen.csv   {len(ev)} 行 / 声明域 {int(ev.in_declared_domain.sum())}")
    print(f"  candidates_28pairs  {len(cand)} 对 / 分量 {cand.component_id.nunique()} "
          f"/ kept {(cand.kept_for_claims == 'kept').sum()} "
          f"/ removed {(cand.kept_for_claims == 'removed').sum()}")
    print(f"  主候选 ΔDM/σ_floor  = {got:.4f}（与产物一致）")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
