# frbpp — reproduction release

Point-process statistics for one-off FRB pairs in CHIME/FRB Catalogue 2, packaged so that every
number in the manuscript can be recomputed on a clean machine from public catalogue files only.

## Install

```bash
python -m venv .venv
.venv/bin/python -m pip install -r requirements-verified.txt   # versions everything was run with
.venv/bin/python -m pip install -e .                            # or: pip install -e ".[test]"
.venv/bin/python -m pytest -q
```

`requirements-verified.txt` pins the exact build used to produce the manuscript (Python 3.13,
numpy 2.5.3, pandas 3.0.6, astropy 8.0.1, matplotlib 3.11.2). `pyproject.toml` alone gives
unpinned minimum versions; use the lock if you want bit-identical output.

On Windows set `PYTHONIOENCODING=utf-8` before running the chain (some steps print `m2` units).

## Run

```bash
bash scripts/run_all.sh              # full chain, logs per step under OUT/chain/
```

`TUTORIAL.md` is the detailed bilingual (EN/中文) hands-on guide: every step, every quantity, troubleshooting and limits.
`RUNBOOK.md` lists the steps one by one, what each reads and writes, its seed, and its measured
runtime. The derived data land in `data/interim/` and `data/processed/upgrade/`;
`scripts/repro_package.py` then collects the artefacts the text cites into `reproduction/` and
aborts if any headline count disagrees.

## Layout

| Path | Contents |
| --- | --- |
| `src/frbpp/` | library: catalogue readers, point-process layer, exposure, upgrade side-channel, validation |
| `scripts/` | one script per step of the analysis plus the figure generators; `run_all.sh` is the chain |
| `data/raw/` | the public inputs actually used, with `chime/SOURCES.json` recording URL, bytes and sha256 |
| `data/interim/`, `data/processed/` | created by the chain |
| `tests/` | unit and smoke tests; one test skips when the HEALPix exposure map is absent |
| `config/` | credentials for the optional online checks; only `*.example.json` are shipped |

## What is and is not included

Every input of the headline chain is a public catalogue file and is shipped here (28 MB):
CHIME Catalogue 2 (`chimefrbcat2.csv`), Catalogue 1 Table 2, the R3 TOA/exposure files, the
Catalogue 1 injections pickle, and the Blinkverse source table used for the external check.

Not shipped: the CHIME Catalogue 1 HEALPix exposure maps (190 MB) and the full injections HDF5
(1.7 GB). They feed the stage-2 exposure and efficiency work, not the manuscript's numbers.
Recover them with `scripts/fetch_chime.py`.

The online checks (`ads_search_565d.py`, `tns_search_565d.py`, `host_565d.py`, `gcn_live.py`)
reach ADS, TNS, MAST and GCN and need credentials from `config/*.example.json`. Their outputs are
archived under `data/processed/upgrade/` and `reproduction/logs/`, so the chain does not have to
re-run them.
