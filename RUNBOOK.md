# RUNBOOK

Ordered reproduction of every number in "An isolated candidate same-source FRB pair hidden in the
apparent one-off population". One command: `bash scripts/run_all.sh`. Each step is a plain script
that reads the artefacts of the steps before it and writes its own; `OUT/chain/<step>.log` keeps
the console output.

Runtime measured in a clean copy of this tree (Python 3.13, 8 cores, catalogue on local disk):
steps 1-9 and 12-17 are seconds, the null battery is ~15 minutes, the same-beam audit ~11
minutes, and the whole chain runs in ~28 minutes.

## Pool conventions the numbers depend on

| Quantity | Value used by the text |
| --- | --- |
| Catalogue rows | 5608 after the Cat1/Cat2 join |
| Distinct events (event-layer de-dup, logic v1.3) | 5014 |
| Apparent one-offs | 3742 |
| Declaration domain (finite DM and coordinates) | 3665 |
| Candidate pairs | 28, the frozen criterion set; the 30-pair union is the stricter null bar |
| Events in those pairs / independent hypotheses | 21 / 8 (7 surviving after the catalogue-grouping audit) |
| DM error floor | s(DM) = 1.4826 MAD(D), D = dm_fitb - bonsai_dm, local window +-0.25 DM |
| Two-burst sigma | sqrt(2) s (conservative); the loose s is reported alongside |
| Permutation draws | B = 1000 per model; p is a resolution limit p <= 1/(B+1), never a measurement |
| Seeds | 20261009 (null battery, robustness), 20261011 (error-convolved same-beam audit, floor bootstrap) |

## Steps

20 steps.

Each step reads the artefacts of the steps before it and writes its own. Exit codes and
seconds are taken from one clean-room invocation of `scripts/run_all.sh` after every derived file
was deleted, so this table cannot drift from the runner: it is regenerated from
`run_all.sh` and `OUT/chain` logs. Full descriptions of each step are in
[TUTORIAL.md § 4](TUTORIAL.md#4).

| # | Script (as invoked by `run_all.sh`) | Exit | Measured |
| --- | --- | --- | --- |
| 01 | `build_interim.py` | exit 0 | 3 s |
| 02 | `event_table_dedup_audit.py` | exit 0 | 0 s |
| 03 | `rerun_v13_event_dedup.py` | exit 0 | 2 s |
| 04 | `gate_g4b.py` | exit 0 | 4 s |
| 05 | `formalize_c2.py` | exit 0 | 3 s |
| 06 | `c2_external_v2.py` | exit 0 | 0 s |
| 07 | `rank_c2_candidates.py` | exit 0 | 3 s |
| 08 | `dm_sigma_floor.py` | exit 0 | 1 s |
| 09 | `dm_floor_robust.py` | exit 0 | 1 s |
| 10 | `c2_null_battery.py` | exit 0 | 872 s |
| 11 | `same_beam_test.py --draws 400 --seed 20261011` | exit 0 | 639 s |
| 12 | `same_beam_coupled_null.py --draws 20000 --seed 20261011` | exit 0 | 8 s |
| 13 | `dm_fingerprint_check.py` | exit 0 | 1 s |
| 14 | `build_isolation_data.py` | exit 0 | 4 s |
| 15 | `known_repeater_control.py` | exit 0 | 2 s |
| 16 | `isolation_conv_565d.py` | exit 0 | 3 s |
| 17 | `scatter_565d.py` | exit 0 | 0 s |
| 18 | `period_565d.py` | exit 0 | 1 s |
| 19 | `dm_z_inversion.py` | exit 0 | 1 s |
| 20 | `repro_package.py` | exit 0 | 1 s |

Step 20 is the guard: it aborts unless the candidate table is exactly the 28 frozen pairs, the
component decomposition is exactly 8, the audited-away group is exactly the 21 pairs of the
seven-event component, and the primary pair's floor consistency matches step 08.

Steps 02 and 04 are not optional: `formalize_c2.py` reads `g4b_report.json`, and
`rank_c2_candidates.py` reads the `c2_candidates_v2.parquet` written by step 06.

## Figures

```bash
python scripts/make_figs.py --only fig2,fig3,fig4,fig5,fig6 --dpi 900 --outdir docs/figures/hires
python scripts/make_fast_advantage.py --dpi 900 --outdir docs/figures/hires
python scripts/make_layout_figs.py  --outdir docs/figures/vector2     # dual.pdf, geometry.pdf
python scripts/make_layout_figs2.py                                   # scatter/flow/roadmap.pdf
```

They read the artefacts above (`c2_null_battery.json`, `event_table_audit.json`,
`pair565_isolation.json`, `pair565_scatter.json`, `c2_isolation_figure_v1.json`,
`g4b_rerun_20261009_v13.json` and the interim tables), so run the chain first. Each generator
draws at a 15.0 cm design width with a **7.30 pt type floor**, which is 7.06 pt once the paper
scales it into its 14.5 cm text block. The floor is enforced, not requested: `fit()` clamps at
the floor instead of shrinking further and the run exits non-zero if any string cannot hold it,
while the shared gate in `fig_gate.py` exits non-zero if any block overprints a neighbouring
label - text inside inset axes included - crosses the canvas or sits on top of data points,
and writes nothing in that case.
Rendering is deterministic: two independent runs produce byte-identical PNGs. The vector
PDFs embed a creation timestamp, so their reproducibility check is per-span text and size
equality against the copies the paper embeds (5/5 match), not bytes. Fonts are the
bundled DejaVu Sans; nothing here needs a system font.

The manuscript embeds `docs/figures/hires/*.png` and `docs/figures/vector2/*.pdf`. The submission
conversion path (`convert_figs.py`, which writes `figures_submission/*.pdf` and `*.tif`) also
feeds `build_retext.py`; the retext step is the one part of the figure pipeline that needs
Windows fonts (`msyh.ttc`, `times.ttf`) and is not on the reproduction path.

## Optional online checks (not needed to reproduce the text)

`reproduction/logs/` needs `tns_search_ledger.jsonl` and `tns_search_summary.md`, which are the
archived outputs of the credentialed TNS query; they ship under `data/processed/upgrade/` so that
step 18 completes offline. They cannot be regenerated without credentials.

`fetch_chime.py` and `fetch_external_catalogs.py` redownload the public inputs. `ads_search_565d.py`,
`rm_ads_565d.py`, `verify_refs_ads2.py`, `tns_search_565d.py`, `tns_object_details.py`,
`host_565d.py` and `gcn_live.py` contact ADS, TNS, MAST and GCN; they need `config/*.local.json`
(see the `.example.json` files). Their archived outputs ship with this release, so the audit trail
is readable without credentials.

`scripts/_tier.txt` marks which scripts are on the reproduction path and which are one-off probes
kept for the audit trail.


## Checking the text against the recomputed artefacts

```bash
python scripts/audit_manuscript_numbers.py --tex path/to/raa_manuscript_revised.tex
```

It reads each headline quantity from the artefacts above and looks for it near a topic anchor in
the manuscript, comparing at the number of decimals the text prints (so `23.86` in the artefact
and `23.9` in the text agree, while `5014` and `5045` do not). Quantities are reported as
`ok`, `MISMATCH` (the text states a different value), or `NOT-STATED`.

## Criterion guard

`formalize_c2.py` keeps the grid configurations up to DM +-5 as the candidate pool and archives
pairs supported only by the DM +-20 corner in `c2_candidates_edge_v1.csv`. With that filter the
chain reproduces the pool the text claims: 28 pairs, 21 events, 8 components, 7 cross-day pairs,
8 hypotheses. Without it the nine-configuration union (30 pairs, 25 events, 10 components) is
written into `c2_candidates_v1.parquet` and `repro_package.py` stops with
`AssertionError: expected 28 pairs`. The stricter null bar is unaffected either way:
`c2_null_battery.py` rebuilds the union from the grid itself.

## Resolved: the DM-z relation behind z = 1.52

`dm_z_inversion.py` used to invert the naive fully ionized integral with its own prefactor
`3 c H0 Omega_b / (8 pi G m_p) = 1118.7 pc cm-3`, which maps the adopted `DM_IGM = 1671.15` to
z = 1.35. The manuscript writes the relation it actually used (Section 4.4: "DM_IGM(z) = 1000 *
integral of (1+z')/E(z')", Planck 2018, macquart2020). With that prefactor the inversion returns
**z = 1.517 -> 1.52**, and the propagated bands **[1.21, 1.84]** (+-20%) and **[1.05, 2.01]**
(+-30%) reproduce digit for digit, together with the statement that the catalogue's +-0.5 dmpc
formal error moves z by <0.001. The script now uses K = 1000 and reports the derived prefactor
beside it, so the 12% difference between the two normalizations is visible rather than silent.

Both consequences were then corrected in the manuscript (2026-10-11):

1. The DM budget table now states `formal DM error <0.001` instead of `formal +-0.02`, matching both
   the inversion and the body of the same section.
2. Section 4.4 now states the calibration domain: the homogeneous relation sits 45 / 38 / 32% above
   the measured `DM_IGM` of the localised bursts at z = 0.189 / 0.321 / 0.660, which biases the
   inversion toward smaller redshifts, so z ~ 1.5 is the lower edge of a cosmological dispersion
   measure rather than a measured redshift, and H2 is excluded by the local-medium and environment
   tests rather than by this number. The per-anchor deviations are stored in
   `dm_z_inversion.json`, and `audit_manuscript_numbers.py` checks all three as text claims.

## Manuscript corrections this release accompanies

This repository ships code and derived data only, not the manuscript; the corrections below were made
to the accompanying paper and are recorded here so the two stay linked.

| Item | Audit claim | Checked against | Outcome |
| --- | --- | --- | --- |
| Petroff et al. 2022 | `A&ARv 30, 25` is wrong | ADS `2022A&ARv..30....2P` (volume 30, article 2) | corrected to `30, 2` |
| Zhang 2020 cited for "29 hits" of RM literature | a mechanisms review cannot back a literature count | ADS `2020Natur.587...45Z` = Zhang B., Nature 587, 45, "The physical mechanisms of fast radio bursts"; the count comes from `pair565_rm_ads.json` (`full:"FRB20201124A" AND full:"rotation measure"` = 29) | the record is real and correctly described; the citation was in the wrong place, so it is dropped and the counts stand as our own ADS audit |
| "Pan-STARRS DR2 (Chambers et al. 2016)" | DR2 should be Pan-STARRS1, or add Flewelling 2020 | the cone search runs on MAST `panstarrs/dr2/mean.csv`, i.e. the DR2 stacked catalogue, so the release name is right; ADS `2020ApJS..251....7F` = Flewelling et al. 2020, ApJS 251, 7, "The Pan-STARRS1 Database and Data Products" | kept "DR2", added Flewelling et al. 2020 beside Chambers et al. 2016 |
| Chen et al. 2026 `RAA 26, 115016` | volume/article not independently confirmed | ADS `2026RAA....26.115016C` | confirmed published; no change |
| Catalogue 2 "at burst level (5045 rows)" | 4539 bursts in the official abstract | our own CSV: 5045 rows = 4539 bursts + 506 sub-burst rows, 37 unlocalized | reworded to rows/bursts/sub-bursts |
| formal DM error on z | table said `+-0.02` | `dm_z_inversion.py` gives 0.0005; the section's own prose says `<0.001` | corrected to `<0.001` |

<a id="verification"></a>

## Verification actually performed

Run in a copy of this tree at a different path with a fresh virtual environment built from
`requirements-verified.txt` (Python 3.13, numpy 2.5.3, pandas 3.0.6, astropy 8.0.1):

| Check | Result |
| --- | --- |
| `pytest -q` | 50 passed, 1 skipped (the skipped test needs the 190 MB HEALPix exposure map) |
| full chain, 20 steps | 20/20 exit 0 (`CHAIN_RC=0`) in one invocation; `repro_package.py` asserts 28 pairs / 8 components / 7 kept / 21 removed and the floor ratio 0.119 |
| pair-level chance rate under a coupling-preserving null | cross-day joint box lambda = 0.0011 +/- 0.0331, P(>=1) = 0.0011 over 20000 DM permutations; randomising right ascension instead gives lambda = 0, so the coupling-preserving null is the conservative one |
| rebuilt `reproduction/` vs the shipped one | tables byte-identical (`events_frozen.csv` 5014 rows, `candidates_28pairs.csv` 28, `dm_floor_bins.csv` 12) and both TNS log files byte-identical; 794 JSON leaves compared, the only differences are 4 `built_at` timestamps |
| `audit_manuscript_numbers.py` | 74 of 74 quantities agree with the recomputed artefacts, including the anchor deviations and the catalogue row/burst split |
| figures regenerated from the recomputed artefacts | two independent renders of the 6 raster panels at 900 dpi are byte-identical (sha256 match), and those are the files the 38-page paper embeds; every raster string holds the 7.30 pt floor (7.06 pt as printed) with 0 overprints, and the 5 vector panels measure 7.3 pt per PDF glyph. Holding the floor re-laid the rasters out on purpose - fig2 heading and note wording, fig3 schematic headings plus the pair label moved into the empty band, fig4 legend side and the P(...) wrap, fig8 headroom above the longest bar - and no drawn number changed |
| timings | null battery 872 s, same-beam audit 639 s, coupling-preserving null 8 s, fingerprint recheck 1 s, everything else seconds; whole chain ~27 min |
