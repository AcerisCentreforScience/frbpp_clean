# Reproduction package

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
| `json/coupled_null.json` | `scripts/same_beam_coupled_null.py` (chance rate under a coupling-preserving DM permutation) |
| `json/dm_fingerprint.json` | `scripts/dm_fingerprint_check.py` (offline recomputation of the dispersion-measure fingerprint distances) |
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
