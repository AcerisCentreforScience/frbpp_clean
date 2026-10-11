# frbpp — reproduction package for an isolated candidate same-source FRB pair

Reference implementation and reproduction package for the analysis behind
*"An isolated candidate same-source FRB pair hidden in the apparent one-off population"*.
It holds the analysis library, one script per step of the calculation, the unit tests, the public
catalogue inputs the study depends on, and the derived tables the study cites. One command
reproduces every published number from public data, with no private inputs and no network access.

## Documentation

| Document | Contents |
| --- | --- |
| **[TUTORIAL.md](TUTORIAL.md)** | full bilingual (English / 中文) guide: environment, the 18 calculation steps one by one, where each published number comes from, figure generation, troubleshooting, scope and limits |
| [RUNBOOK.md](RUNBOOK.md) | compact step table with inputs, outputs, seeds and measured runtimes; the frozen analysis conventions; the verification record |
| [TUTORIAL.md § 1](TUTORIAL.md#1) | three-minute start on Linux/macOS and Windows |
| [TUTORIAL.md § 6](TUTORIAL.md#6) | trace one specific published number back to its artefact field |
| [TUTORIAL.md § 11](TUTORIAL.md#11) | what this release does *not* establish |
| [scripts/_tier.txt](scripts/_tier.txt) | which scripts are on the reproduction path and which are archived one-off probes |
| [INPUTS.sha256](INPUTS.sha256) | hashes of every shipped input file |
| [LICENSE](LICENSE) · [CITATION.cff](CITATION.cff) | licence, and the scope note for third-party data · software citation metadata |
| [data/raw/chime/SOURCES.json](data/raw/chime/SOURCES.json) | download URL, size, hash and retrieval date of each public catalogue file |
| [reproduction/](reproduction/) | derived tables and JSON artefacts cited by the study, with `MANIFEST.csv` |

## Quick start

```bash
git clone https://github.com/AcerisCentreforScience/frbpp_clean.git && cd frbpp_clean
python3 -m venv .venv
.venv/bin/python -m pip install -r requirements-verified.txt
.venv/bin/python -m pip install -e .
.venv/bin/python -m pytest -q          # 50 passed, 1 skipped
bash scripts/run_all.sh                # 20 steps, ~27 min, per-step logs under OUT/chain/
```

`requirements-verified.txt` pins the exact build the results were produced with (Python 3.13,
numpy 2.5.3, pandas 3.0.6, astropy 8.0.1, matplotlib 3.11.2); `pyproject.toml` alone gives
unpinned minimum versions. On a Windows console set `PYTHONIOENCODING=utf-8` before running
individual scripts — some print `m²` units (the runner already exports it).

## Layout

| Path | Contents |
| --- | --- |
| `src/frbpp/` | library: catalogue readers, point-process engine, exposure model, upgrade side-channel, validation |
| `scripts/` | one script per analysis step, the figure generators, the manuscript auditors, and `run_all.sh` |
| `tests/` | unit and smoke tests; one test skips when the HEALPix exposure map is absent |
| `data/raw/` | the public inputs actually used (28 MB), with recorded provenance |
| `data/interim/`, `data/processed/` | created by the chain |
| `reproduction/` | shipped derived tables, JSON artefacts and their manifest |
| `config/` | `*.example.json` templates for the optional online checks; no credentials are shipped |

## What is and is not included

Shipped: CHIME/FRB Catalogue 2 (`chimefrbcat2.csv`), Catalogue 1 Table 2, the R3 TOA and exposure
files, the Catalogue 1 injections pickle, and the Blinkverse source table used for the external
cross-check — every input the headline chain reads.

Not shipped: the Catalogue 1 HEALPix exposure maps (190 MB) and the full injections HDF5 (1.7 GB).
They serve the exposure and efficiency work, not the published numbers; recover them with
`scripts/fetch_chime.py`.

The literature and object-catalogue checks (`ads_search_565d.py`, `tns_search_565d.py`,
`host_565d.py`, `gcn_live.py`) contact ADS, TNS, MAST and GCN and need credentials copied in from
the `config/*.example.json` templates. Their archived outputs ship with this release, so the
chain never has to reach the network.

## Verification

Rebuilt on a clean machine from a fresh environment: the 20 steps complete with exit code 0, the
regenerated `reproduction/` matches the shipped tables byte for byte and the shipped JSON
artefacts in every value except four `built_at` timestamps, the six raster figures regenerate
byte-for-byte across independent renders while holding a 7.30 pt type floor (7.06 pt as
printed), and `scripts/audit_manuscript_numbers.py`
confirms 74 of 74 published quantities against the recomputed artefacts. Details in
[RUNBOOK.md](RUNBOOK.md#verification) and [TUTORIAL.md § 12](TUTORIAL.md#12).

## Licence and citation

The code, tests and documentation are released under the [BSD 3-Clause licence](LICENSE), the same
licence used by the surrounding scientific Python stack; the third-party catalogue files under
`data/raw/` are not relicensed by it and remain subject to their own terms of use, as recorded in
[LICENSE](LICENSE) and in [data/raw/chime/SOURCES.json](data/raw/chime/SOURCES.json).

Software citation metadata is in [CITATION.cff](CITATION.cff); the release is tagged
[`v1.1.0`](https://github.com/AcerisCentreforScience/frbpp_clean/releases/tag/v1.1.0). The
companion paper reference is added here when it appears.
