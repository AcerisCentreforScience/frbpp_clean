# frbpp reproduction package — detailed guide
# frbpp 复现发布包 — 详细说明文档

This document is the long-form, hands-on companion to
[README.md](README.md) and the compact [RUNBOOK.md](RUNBOOK.md). It states, for every published
number, which script produces it, which artefact field holds it, and how to check it; it also
states plainly what this release does not establish.
本文件是 [README.md](README.md) 与简版 [RUNBOOK.md](RUNBOOK.md) 的详版配套文档：对每一个发表过的数字，
写明由哪个脚本产生、存在哪个产物字段、如何核对；同时也明确写出这份发布**不能**证明什么。

The repository contains code, the public inputs, and the derived artefacts only. The manuscript
itself is not part of this repository; the auditors in section 8 run against a local copy of it.
仓库只含代码、公开输入与派生产物；论文正文不在本仓库内，第 8 节的核对脚本是对本地论文副本运行的。

## Contents / 目录

0. [What this release reproduces / 能复现什么](#0)
1. [Quick start / 三分钟跑起来](#1)
2. [Layout / 目录结构](#2)
3. [Inputs / 输入数据](#3)
4. [The 18 steps, one by one / 十八步逐步详解](#4)
5. [Frozen conventions / 冻结口径](#5)
6. [Re-deriving one number / 想核对某一个数](#6)
7. [Regenerating the figures / 重出图](#7)
8. [Manuscript-side audits / 论文侧核对脚本](#8)
9. [Optional online checks / 在线核查（可选）](#9)
10. [Troubleshooting / 故障排查](#10)
11. [Honest limits / 边界](#11)
12. [How this release was verified / 验证方法](#12)
13. [Obtaining and archiving / 获取与存档](#13)
14. [Glossary / 术语对照](#14)
15. [FAQ / 常见问题](#15)
16. [Design notes / 设计与可复现性要点](#16)

---

<a id="0"></a>

## 0. 这份东西能复现什么 / What this release reproduces

**EN.** Everything the manuscript *An isolated candidate same-source FRB pair hidden in the
apparent one-off population* states as a number, except (a) results that require live credentials
(ADS / TNS / MAST / GCN) and (b) one denominator (58/85) that mixes Catalogue 1-only sources.
Headline quantities: the event pool 5608 → 5014 → 3742 → 3665; the candidate set 28 pairs / 21
events / 8 components / 7 surviving hypotheses; θ = 23.9″; |ΔDM| = 0.18 pc cm⁻³ = 0.12σ on the
empirical floor; B = 1000 permutation nulls with p ≤ 0.001 and Bonferroni×9 global p ≤ 0.009; the
same-beam and correlated-error audit; the scattering discriminant; the DM budget and its redshift;
the 565 d ↔ 567.002 sidereal-day quantization; and all 11 figures.

**中文.** 稿件里以"数字"出现的东西全部可复现，只有两类例外：需要凭证的在线核查（ADS / TNS / MAST /
GCN），以及混入了"仅 Cat1 源"的那个分母（58/85）。核心量包括：事件池 5608 → 5014 → 3742 → 3665；
候选集 28 对 / 21 事件 / 8 分量 / 7 个存活假设；θ = 23.9″；|ΔDM| = 0.18 pc cm⁻³，按经验地板 = 0.12σ；
B = 1000 的四类置换零假设（p ≤ 0.001，Bonferroni×9 全局 p ≤ 0.009）；同束与相关误差审计；散射判别量；
DM 预算及其红移；565 d ↔ 567.002 恒星日的量化；以及全部 11 张图。

**Verified / 实测状态.** 18/18 steps exit 0 in one invocation on a different drive with a fresh
virtual environment; the rebuilt `reproduction/` matches the shipped one byte-for-byte except four
`built_at` fields; the six raster figures regenerate pixel-identical to the compiled PDF (maxdiff 0);
`audit_manuscript_numbers.py` reports **67/67** agreeing quantities.
在换盘 + 全新虚拟环境下 18/18 步 exit 0；重建的 `reproduction/` 与随包件**只差 4 个 `built_at`**；
6 张位图重出后与论文 PDF 内嵌图逐像素相同；`audit_manuscript_numbers.py` **67/67 一致**。

---

<a id="1"></a>

## 1. 三分钟跑起来 / Quick start (three minutes)

### 1.1 前置 / Prerequisites

| | |
| --- | --- |
| Python | 3.11+ (verified on 3.13.15) / 已在 3.13.15 上验证 |
| Disk | ~1 GB free for `.venv`; the shipped inputs are 28 MB / 磁盘约 1 GB，随包输入 28 MB |
| CPU time | ~26 min single-run for the whole chain (16 min of it is one script) / 全链约 26 分钟（其中一步占 16 分钟） |
| Network | **not needed** for the chain; only for `fetch_*.py` and the online checks / 主链**不需要联网** |

### 1.2 Linux / macOS

```bash
tar xzf frbpp_clean.tar.gz && cd frbpp_clean
python3 -m venv .venv
.venv/bin/python -m pip install --upgrade pip
.venv/bin/python -m pip install -r requirements-verified.txt   # exact versions, bit-reproducible
.venv/bin/python -m pip install -e .
.venv/bin/python -m pytest -q                                   # 50 passed, 1 skipped
bash scripts/run_all.sh                                         # 18 steps, logs in OUT/chain/
```

### 1.3 Windows (PowerShell)

```powershell
cd frbpp_clean
python -m venv .venv
.\.venv\Scripts\Activate.ps1
pip install -r requirements-verified.txt
pip install -e .
python -m pytest -q
$env:PYTHONIOENCODING="utf-8"   # 必须：有些步骤会打印 m² 等字符 / required, some steps print m²
& 'C:\Program Files\Git\bin\bash.exe' scripts/run_all.sh
```

> **Why PYTHONIOENCODING / 为什么必须设它.** On a Chinese Windows console the default codec is GBK;
> `rank_c2_candidates.py` prints `m²` and dies with `UnicodeEncodeError` unless UTF-8 is forced.
> `scripts/run_all.sh` already exports it, so the failure only appears if you run scripts by hand.
> 中文 Windows 控制台默认 GBK，`rank_c2_candidates.py` 打印 `m²` 会崩；`run_all.sh` 已内置该变量，
> 手工单跑脚本时才需要自己设。

### 1.4 期望输出 / Expected output

```
--- [10_null_battery] exit=0  906s
--- [11_same_beam] exit=0  633s
--- [18_repro_package] exit=0  1s
CHAIN_RC=0
```
Step 18 prints the guard summary / 第 18 步会打印守卫结果：
```
candidates_28pairs  28 pairs / 8 components / kept 7 / removed 21
primary pair DeltaDM/sigma_floor = 0.1191 (matches the artefact)
```

---

<a id="2"></a>

## 2. 目录结构 / Layout

```
frbpp_clean/
├── README.md                  仓库首页：范围、上手、文档索引 / repository front page
├── RUNBOOK.md                 步骤表 + 口径 + 已修的稿件问题 / step table + conventions + corrections
├── TUTORIAL.md                本文件：详细说明 / this document
├── pyproject.toml             包定义（src 布局）/ package definition (src layout)
├── requirements-verified.txt  实测环境的精确版本锁 / exact versions everything was run with
├── INPUTS.sha256              10 个输入文件的哈希 / hashes of the 10 shipped inputs
├── .gitignore                 忽略派生数据与本地凭证 / ignores derived data and credentials
├── config/                    只有 *.example.json / examples only, no secrets
├── data/
│   ├── raw/                   公开目录输入（28 MB）+ chime/SOURCES.json 出处 / public inputs + provenance
│   ├── interim/               第 1 步生成 / produced by step 1
│   └── processed/upgrade/     第 2–17 步生成；随包附 2 个 TNS 归档件 / steps 2-17 + two archived TNS logs
├── src/frbpp/                 库：目录接入、点过程、暴露、升级侧信道、验证 / the library
├── scripts/                   72 个脚本 + run_all.sh + _tier.txt 分级表 / scripts + runner + tier list
└── tests/                     13 个测试文件 / 13 test modules
```

`scripts/_tier.txt` classifies every script as `pipeline` / `figures` / `online` / `stage-gate` /
`one-off`. Only `pipeline` + `figures` are on the reproduction path.
`scripts/_tier.txt` 把每个脚本分级为 `pipeline`（主链）/ `figures`（出图）/ `online`（需凭证或网络）/
`stage-gate`（引擎自有验收门）/ `one-off`（一次性探针）。只有前两类在复现路径上。

---

<a id="3"></a>

## 3. 输入数据 / Inputs

### 3.1 随包发布 / Shipped (28 MB)

| 文件 File | 大小 size | 用途 used by | 出处 provenance |
| --- | --- | --- | --- |
| `data/raw/chime/cat2_25/chimefrbcat2.csv` | 4,057,396 B (3.9 MiB) | 所有主表 / every headline table | CHIME/FRB Catalogue 2, CANFAR vault (see `SOURCES.json`) |
| `data/raw/chime/05_cat1/cat1_table2_vizier.tsv` | 289,394 B | Cat1 行、去重对照 / Cat1 rows, de-dup control | VizieR |
| `data/raw/chime/03_180916/R3_CHIME_TOA_20Feb4.txt` | 535 B | 第 1 步健全性检查 / step-1 sanity | CANFAR (AstroDataCitationDOI) |
| `data/raw/chime/03_180916/R3_CHIME_exposure.npz` | 853,072 B | 同上 / as above | CANFAR |
| `data/raw/chime/06_injections/chimefrb_catalog1_injections.p` | 18,614,740 B | 第 1 步读取形状 + 哈希 / read for shape + hash in step 1 | CANFAR |
| `data/raw/blinkverse/blinkverse_sources.json` | 2,559,111 B | 第 6 步外部核查 / step 6 external probe | Blinkverse mirror API; sha256 `3302a12f…` matches the working log |
| `data/raw/blinkverse/blinkverse_hosts.json` | 284,071 B | 同上 / as above | Blinkverse |
| `data/raw/fast/fast_20220912A_bursts.json` | 473,057 B | 一次性跨仪器 DM 偏移 / one-off DM-offset script | published FAST table |
| `data/raw/apertif/apertif_frbs.json` | 19,641 B | 一次性新批次扫描 / one-off batch scan | VizieR J/A+A/693/A279 |
| `data/raw/gcn/history_20260905_20261008.jsonl` | 87,797 B | 一次性 GCN 核查 / one-off GCN check | GCN Kafka retention |

Verify them / 逐文件核对：

```bash
python - <<'PY'
import hashlib, pathlib
for ln in pathlib.Path("INPUTS.sha256").read_text().split():
    pass
PY
python -c "import hashlib,pathlib;[print(hashlib.sha256(p.read_bytes()).hexdigest()[:16],p) for p in sorted(pathlib.Path('data/raw').rglob('*') ) if p.is_file()]"
```
Or diff against `INPUTS.sha256` line by line / 或直接与 `INPUTS.sha256` 对比。

### 3.2 不随包发布 / Not shipped (and why)

| 未包含 Not included | 大小 size | 影响 chain impact | 取回 how to get |
| --- | --- | --- | --- |
| Cat1 HEALPix exposure maps | 190 MB | 只影响 `test_exposure_model.py` 的一条测试（会 skip）/ one skipped test only | `python scripts/fetch_chime.py` |
| Full injections HDF5 | 1.7 GB | 只影响阶段 2 效率标定，不进正文 / stage-2 efficiency, not the text | `python scripts/fetch_chime.py` |

### 3.3 出处与许可 / Provenance and licensing

`data/raw/chime/SOURCES.json` records, per file, the download URL, byte count, sha256 and retrieval
timestamp. `INPUTS.sha256` does the same for the rest of the shipped inputs. All are public
catalogue releases; nothing here depends on private or embargoed data.
`SOURCES.json` 记录每个下载文件的 URL、字节数、sha256 与取回时间；其余输入由 `INPUTS.sha256` 覆盖。
全部为公开发表目录，不依赖任何私有或禁运数据。

---

<a id="4"></a>

## 4. 十八步链路逐步详解 / The 18 steps, one by one

Each step is an ordinary script: it reads the artefacts of earlier steps and writes its own.
`OUT/chain/<NN_name>.log` keeps its console output; `run_all.sh` records exit codes.
每一步都是普通脚本：读前序产物、写自己产物；`OUT/chain/<NN_name>.log` 是它的控制台输出。

### 步骤 01 — `build_interim.py` 规范化原始数据 / normalise raw files

| | |
| --- | --- |
| 读 reads | 4 个 raw 文件 / the four raw CHIME files |
| 写 writes | `data/interim/cat1_bursts.parquet`, `cat2_bursts.parquet`, `r3_180916_toa.parquet`, `_build_manifest.json` |
| 期望 expected | `cat1: 600 rows`, `cat2: 5045 rows`, `r3_toa: 38 rows`, `injections 84697 × 25` |
| 意义 why it matters | 契约校验在这里发生（列、dtype、NaN 不伪造）/ contract validation, NaN never fabricated |

Cat2 的 **5045 行**是发布 CSV 的行数：其中 **4539 是不同暴**，**506 是子暴行**，**37 行无坐标**。
这三个数在 §6 里逐条核对；稿件正文现在写的就是这个拆分。
The 5045 rows of Cat2 = 4539 distinct bursts + 506 sub-burst rows, 37 of them unlocalized.

### 步骤 02 — `event_table_dedup_audit.py` 事件池口径 / the event pool

| | |
| --- | --- |
| 读 reads | interim |
| 写 writes | `upgrade/event_table_audit.json` |
| 期望 expected | `row_level.all = 5608`, `dedup.all = 5014`, `dedup.one_off = 3742`, `dedup.one_off_finite_dm = 3665`, `cat1_rows_matched_in_cat2 = 594` |
| 口径 convention | 同名且 `\|Δt\| < 0.05 d` 判为同一暴 / same name and \|Δt\| < 0.05 d means the same burst |

这条给出正文 §2 的四个数。它同时说明**为什么行级计数会虚高**：Cat1 有 594 行在 Cat2 里重遇。
These are Section 2's counts, and the 594 repeated Cat1 rows explain why row-level counts inflate.

### 步骤 03 — `rerun_v13_event_dedup.py` v1.3 事件层去重 / event-layer de-duplication

Writes `g4b_rerun_20261009_v13.json`; expected `events 5014`, `one-off 3742`, `finite DM 3665`,
`sig cfg 9/9`, `union 30 pairs`. 这一步是"同一份判据在 v1.2 与 v1.3 两种去重层下的对照"。
It reports the same criterion under both de-duplication layers.

### 步骤 04 — `gate_g4b.py` 判据网格（**README 曾漏登记**）/ criterion grid (was missing from the old README)

Writes `g4b_report.json`, which step 05 reads. **跳过这一步，第 5 步会 `FileNotFoundError`。**
The 3×3 grid is θ ∈ {60″, 180″, 600″} × DM window ∈ {1, 5, 20} pc cm⁻³. Its union is 30 pairs.
Skipping this step makes step 5 crash; the union over the nine configurations is 30 pairs.

### 步骤 05 — `formalize_c2.py` 冻结候选池 / freeze the candidate pool

| | |
| --- | --- |
| 写 writes | `c2_candidates_v1.{parquet,csv}`（判据池 **28 对**）、`c2_candidates_edge_v1.csv`（**2 对**边缘）、`c2_formalization.md` |
| 判据 criterion | 保留 DM 窗 ≤ 5 的九档配置；只由 DM±20 支持的对进边缘文件 / keep configurations up to DM ±5 |
| 为什么 why | 论文主张 28 对，30 对并集只作"更严 null bar"；不加这道过滤，第 18 步断言会中止 |

被过滤掉的两对是 / the two archived edge pairs are
`FRB20221208D × FRB20230913A` (ΔDM = +9.25) 和 and `FRB20190519J × FRB20200809C` (ΔDM = −6.59).
它们的唯一支持配置是 θ = 600″ × DM±20。**`|ΔDM| < 5` 与"DM 窗 ≤ 5"给出完全相同的 28 对成员**（差集为空）。
Both are supported only by the loosest corner. `|ΔDM| < 5` and `DM window ≤ 5` select exactly the same 28 members.

### 步骤 06 — `c2_external_v2.py` 外部库核查（**README 曾漏登记**）/ external-catalogue probe

Writes `c2_candidates_v2.{parquet,csv}` + `c2_external_v2.md`; step 07 reads it. 内部有锚点自检：
180916B / 121102A / 20201124A 必须都是 `REPEAT_FIRST`，否则判为失败。
Contains a sanity anchor: three known repeaters must read `REPEAT_FIRST`.

### 步骤 07 — `rank_c2_candidates.py` 排序与连通分量 / ranking and components

Writes `c2_ranked_crossday_v1.*`, `c2_ranked_sameday_v1.*`, `c2_ranked_hypotheses_v1.*`.
Expected 期望：`hypotheses = 8 rows`, cross-day table `7 rows` after step 08.
这一步用并查集把成对合并成"独立物理假设"。/ Union-find merges pairs into independent hypotheses.

### 步骤 08 — `dm_sigma_floor.py` 经验地板重算一致度 / consistency under the empirical floor

Writes `c2_ranked_*_v14.*` 和 and `c2_ranked_dm_sigma_v1.{csv,parquet,md}`, `LOGIC_VERSION = 1.4.0`.
主候选 `r_floorB = 0.119186`（正文写 0.12σ），形式误差口径 `0.354`（正文写 0.35σ）。
σ 的定义见 §5；分母用 Cat2 自测的经验地板，**不用**目录给的形式误差。
The denominator is the catalogue's own empirical floor, never the formal per-burst error.

### 步骤 09 — `dm_floor_robust.py` 地板的稳健性 / floor robustness

| | |
| --- | --- |
| 写 writes | `c2_dm_floor_bins_v1.csv`（12 箱）、`c2_dm_floor_robust_v1.json` |
| 期望 expected | 可比样本 `n_comparable = 4961`；窗 ±0.25 时 `s = 1.0501`、`σ_comb = 1.4851`、比值 `0.1191`，bootstrap CI `[0.0928, 0.1421]` |
| 敏感性 sensitivity | 三窗 {0.15, 0.25, 0.40} × 两约定（保守 √2·s / 宽松 s）；只有**第一名**在所有格子中不变 |
| 形式误差检验 | `96.13%` 的样本 \|D\| 大于其自身形式误差；各箱中位数之比 **7.8–46.8**，全表 **27.2** |

Bootstrap 用 2000 次重抽、种子 `20261011`。/ 2000 bootstrap draws, seed 20261011.

### 步骤 10 — `c2_null_battery.py` 零假设电池（最耗时 16 分钟）/ the null battery

Writes `c2_null_battery.json`. 四类置换（RA 置换含 Dec 分层、DM 洗牌、RA+DM、完全随机），每类 B = 1000；
观测对照条数取九档并集 = **30**（比 28 更严）；`n_ge_real = 0` ⇒ `p ≤ 1/(B+1) = 0.001`，
Bonferroni×9 ⇒ 全局 `p ≤ 0.009`。**p 是分辨率下限，不是测量值**，这一点在正文与产物里都写死了。
`real_cand = 30` is the stricter bar; `p ≤ 0.001` is a resolution limit, never a measured value.

### 步骤 11 — `same_beam_test.py --draws 400 --seed 20261011` 同束与相关误差审计（11 分钟）

Writes `same_beam_test.json`. 三条结论 / three outputs:
1. **|HA| 分布随池与时间尺度变化**：5008 行（`mjd_400`，UTC）中位 `0.0488 h`、`<0.1 h` 占 `78.9%`；声明域
   3665 为 `0.05 h` / `76.0%`；已命名源内部两两中位 `0.0416 h`。改用 `t_bjd_tdb` 则近乎均匀（中位 `11.68 h`）
   ⇒ 赤经由 `mjd_400` 反推，UTC 下 |HA|≈0 属构造。引用任何 |HA| 数**必须同时给池与时间尺度**。
   Quote the pool *and* the timescale; under TDB the statistic is uniform, which is the control.
2. **ρ 扫描**：两暴相对 1σ 从 ρ=0 的 `993.7″` 到 ρ=0.99 的 `110.8″`，θ/σ_rel 只有 `0.024 → 0.215`，
   任何 ρ 下都远大于 θ ⇒ **23.9″ 不是被分辨的角距**，角距只能当选据，不能当证据。
   θ is unresolved at every ρ: it is a selector, not evidence.
3. **同束不是独立条件**：同束代理把跨日池里 θ<24″∧|ΔDM|<0.18 的条件率抬 `1.63×`，
   但**不同束池里 θ<24″ 的对数为 0** ⇒ 小 θ 已蕴含同凌日窗。误差卷积 400 次：全池中位 `3`，
   **跨日子集 `P(≥1) = 0`** ⇒ 小 θ 背景基本是同夜碎片。
   Same-beam is not an independent condition; the load-bearing restriction is cross-day.

### 步骤 12–16 — 背景、对照与三个专用量 / background, control, three specialised quantities

| 步 step | 脚本 script | 产物 output | 正文用在哪 where it is used |
| --- | --- | --- | --- |
| 12 | `build_isolation_data.py` | `c2_isolation_figure_v1.json` | 最近邻 ECDF、计数实验（图 8 背景）/ nearest-neighbour background |
| 13 | `known_repeater_control.py` | `known_repeater_control.json` + `_rows.csv` | §5.4：池 83、召回 45/83、候选在第 `20.5` 百分位 / known-repeater control |
| 14 | `isolation_conv_565d.py` | `pair565_isolation.json` | 图 3 联合孤立性 / joint isolation |
| 15 | `scatter_565d.py` | `pair565_scatter.json` | 图 6：控制细胞 n=83、33 可测（`40%`）、最小 `5.2 ms`、中位 `31 ms`；τ(1835) ≈ `24 ms`；R = `271` / scattering discriminant |
| 16 | `period_565d.py` | `pair565_period.json` | Δt = `565.4537 d` = `567.002` 恒星日（距整数 `164 s`）/ sidereal quantization |

**"565 d"不是周期。** CHIME 漂扫使任意两次探测的 Δt 自动接近整数恒星日，因此全文禁用 "periodicity"，
只说"N 个源中天周期"。别名守卫未把 ρ(565.45)=0.067 判为 flag，这否定的只是"由两暴量周期"。
565 d is not a period; the ban on the word is deliberate and the alias guard is reported honestly.

### 步骤 17 — `dm_z_inversion.py` DM 预算与红移 / DM budget and redshift

Writes `dm_z_inversion.json`. 关系按正文所用形式：`DM_IGM(z) = 1000 ∫₀ᶻ (1+z′)/E(z′) dz′`，Planck 2018；
`DM_obs = 1834.95 = MW 123.8 + halo 40 + IGM 1671.15` ⇒ **z = 1.517 → 正文写 1.52**；
±20% → `[1.21, 1.84]`，±30% → `[1.05, 2.01]`；形式 ±0.5 dmpc 只动 z **0.0005**（故表格写 `<0.001`）；
系统项网格给 z ∈ [1.493, 1.531]（正文写"shifts z by <0.04"）。
The script also reports the prefactor derived from the same cosmology, `1118.7`, and the z it would
give (`1.350`), so the 12% normalization gap is visible rather than silent.
Finally the low-z anchors: the measured DM_IGM of the localised bursts at z = 0.189 / 0.321 / 0.660
sits **45 / 38 / 32%** below this relation — stored per-anchor in the JSON and quoted in §4.4 of the
paper. That direction biases the inversion toward *smaller* z, so 1.5 is a lower edge, not a measurement.

### 步骤 18 — `repro_package.py` 装配复现包（**硬断言**）/ assemble the package, hard asserts

Writes `reproduction/{README.md, tables/, json/, logs/, MANIFEST.csv}`. 断言 / it aborts unless:
候选恰 28 对 / exactly 28 pairs; 分量恰 8 / 8 components; 被剔除恰是那 21 对（7 事件组）/ the 21
removed entries are the audited seven-event group; 主候选地板一致度与第 8 步产物相同 / the primary
ratio matches step 8. `MANIFEST.csv` 记每个文件的字节数与 sha256 / per-file bytes and sha256.

需要 `data/processed/upgrade/tns_search_ledger.jsonl` 与 `tns_search_summary.md`——它们是在线 TNS
查询的归档输出，**随包提供**，因此第 18 步离线也能完成；没有它们，这一步会在复制日志处失败。
The two archived TNS logs ship precisely so that step 18 runs offline.

---

<a id="5"></a>

## 5. 冻结口径（不要漂移）/ Frozen conventions (do not drift)

| 项 Item | 值 Value |
| --- | --- |
| 事件表 Event table | v1.3 事件层去重：5608 行 − 594 重复 Cat1 行 = **5014** 个独立事件 |
| 表观单次 Apparent one-offs | **3742** |
| 声明域 Declaration domain | **3665**（有限 DM + 坐标 / finite DM and coordinates） |
| 候选 Candidate set | **28** 对判据池；**30** 对九档并集仅作更严 null bar |
| 事件/假设 Events / hypotheses | **21** / **8**（污染审计后 **7** 存活） |
| DM 地板 Floor | `s(DM) = 1.4826·MAD(D)`，`D = dm_fitb − bonsai_dm`，局部窗 ±0.25·DM |
| 两暴合成 Two-burst σ | `σ_comb = √2·s`（保守）；宽松 `s` 并列报告 |
| 置换 Permutations | 每模型 B = 1000；`p ≤ 1/(B+1)` 是**分辨率下限**，不是测量值 |
| 种子 Seeds | 20261009（null battery、robustness）；20261011（同束卷积、地板 bootstrap） |
| 目录文件 Catalogue file | 只认 `chimefrbcat2.csv` 这一份；换发布版本必须重算地板统计 |

---

<a id="6"></a>

## 6. 想核对某一个数怎么下手 / How to re-derive one number

```bash
python scripts/audit_manuscript_numbers.py --tex path/to/raa_manuscript_revised.tex \
       --repro reproduction --processed data/processed/upgrade
```
它从产物取每个量、按**正文实际显示的位数**、在**语境锚点**附近比对，输出三态：
It reads each quantity from the artefacts and looks for it near a topic anchor, comparing at the
number of decimals the text prints. Three outcomes:

| 状态 | 含义 |
| --- | --- |
| `ok` | 产物值在该语境里出现（含四舍五入等价，如 23.86 ↔ 23.9） |
| `MISMATCH` | 语境在，但正文写的是另一个数（真不一致） |
| `NOT-STATED` | 语境在，但正文没有主张这个量 |

**为什么不裸查字符串 / why it is anchored.** 第一版我用全文搜串，结果 `1.35` 被表格列宽
`\hsize=1.35\hsize` 满足，成了假绿灯——这类假通过比红叉更危险。
An unanchored string search was satisfied by a table column spec (`\hsize=1.35\hsize`). That is the
failure mode to design against: a false green is worse than a red.

单个量的对应关系 / quantity → source:

| 正文数字 number | 产物 field | 命令 command |
| --- | --- | --- |
| 5608 / 5014 / 3742 / 3665 | `json/event_table_audit.json` | `python scripts/event_table_dedup_audit.py` |
| 28 对 / 21 事件 / 8 分量 | `c2_null_battery.json :: independence` | steps 05, 10 |
| 0.12σ、0.35σ | `c2_ranked_dm_sigma_v1.csv`, `c2_dm_floor_robust_v1.json :: window_frac_0.25` | steps 08, 09 |
| p ≤ 0.001 / 全局 ≤ 0.009 | `json/null_battery_summary.json :: models.*` | step 10 |
| 1.63×、993.7″、110.8″、0.024、0.215、400、中位 3 | `json/same_beam_audit.json` | step 11 |
| 45/83、20.5 百分位 | `json/known_repeater_control.json` | step 13 |
| 40%、5.2、31、24 ms、R = 271 | `pair565_scatter.json` | step 15 |
| 565.45 d、567.002 T_sid、164 s | `pair565_period.json`（恒星日换算 `T_sid = 0.99726958 d`） | step 16 |
| 1.52、[1.21,1.84]、[1.05,2.01]、<0.001 | `dm_z_inversion.json` | step 17 |
| 5045 / 4539 / 506 | `data/raw/.../chimefrbcat2.csv`（`sub_num` 为空或 0 即暴级） | step 01 |

---

<a id="7"></a>

## 7. 重出图 / Regenerating the figures

```bash
python scripts/make_figs.py --only fig2,fig3,fig4,fig5,fig6 --dpi 900 --outdir docs/figures/hires
python scripts/make_fast_advantage.py --dpi 900 --outdir docs/figures
python scripts/make_layout_figs.py  --outdir docs/figures/vector2      # dual.pdf, geometry.pdf
python scripts/make_layout_figs2.py                                    # scatter/flow/roadmap.pdf
```
- 必须先跑完链路：图读 `c2_null_battery.json`、`event_table_audit.json`、`pair565_isolation.json`、
  `pair565_scatter.json`、`c2_isolation_figure_v1.json`、`g4b_rerun_20261009_v13.json` 与 interim 表。
  Run the chain first; every panel is data-driven.
- 出图脚本一律按 **15.00 cm** 文本宽设计（`make_figs.py`、`make_layout_figs.py` 都是 `15.0/2.54`），
  但 `convert_figs.py` 与 `audit_figure_specs.py` 按 **16.2 cm** 折算有效尺寸 ⇒ 同一张图存在两套 dpi/pt 口径，
  两者读数不可直接互比；统一到版心实测值之前，这是已知不一致（见 §11.6）。
  每个生成器保存前重跑叠字/碰撞检查（`fig_gate.py`）；`--min-fs 0` 可退回旧字号行为。
  Generators design at 15.00 cm while conversion and the spec audit scale by 16.2 cm: reconcile before
  comparing their readings.
- 字体只用 matplotlib 自带 DejaVu Sans，**不需要系统字体** ⇒ Linux 上也能出同样的图。
  Bundled DejaVu only, so this works on Linux.
- 唯一需要 Windows 字体（`msyh.ttc` / `times.ttf`）的是 `build_retext.py` + `retext_figures.py`
  的换字支线，**不在正文复现路径上**。
  Only the retext branch needs Windows fonts, and it is off the reproduction path.
- 实测：在"换盘 + 全新环境"下重出的 6 张位图与论文 PDF 内嵌图**逐像素 maxdiff = 0**；
  5 张矢量图所画文字与磁盘版本零差异。
  Measured: 6 rasters pixel-identical to the PDF's (maxdiff 0), 5 vectors text-identical.

---

<a id="8"></a>

## 8. 稿件侧的审计脚本 / Manuscript-side audits

本仓库不含论文源文件（只有代码与派生数据），所以下面三条要对本地的论文副本运行：
This repository ships code and derived data, not the manuscript, so point these at a local copy of it.

```bash
python scripts/check_tex.py   docs/submission/raa_manuscript_revised.tex   # 结构一致性 / structure
python scripts/audit_refs.py  docs/submission/raa_manuscript_revised.tex   # 文献↔引用双向、字母序
python scripts/audit_submission.py raa_manuscript_revised                  # 九项期刊合规
python scripts/audit_manuscript_numbers.py --tex ...                       # 数字↔产物（§6）
```
`audit_submission.py` 需要 `docs/submission/<stem>.tex` 与 `docs/submission/_build/<stem>.pdf`，
即必须先编译（**XeLaTeX**，`pdflatex` 因 fontspec 会致命失败）：
It needs a compiled PDF; use XeLaTeX (pdflatex fails on fontspec):
```bash
PATH="/path/to/MiKTeX/miktex/bin/x64:$PATH"
xelatex -interaction=nonstopmode -output-directory=_build raa_manuscript_revised.tex   # ×3
```
当前状态 / current: 37 页、三遍 0 错误、`check_tex` ALL CHECKS PASSED、期刊合规审计**硬违规 0**、
`audit_refs` 33 条 ↔ 33 键无悬空无未引、字母序正确、`audit_manuscript_numbers` **67/67**。

---

<a id="9"></a>

## 9. 在线核查（可选，需凭证）/ Optional online checks needing credentials

```bash
cp config/ads.local.example.json  config/ads.local.json   # 填入真实值 / fill in real values
cp config/tns.local.example.json  config/tns.local.json
cp config/gcn_kafka.local.example.json config/gcn_kafka.local.json
```
`.gitignore` 已忽略 `config/*.local.json`，凭证不会进包 / never packaged。

三个已知的坑（都实测过）/ three traps, all measured:
1. **TNS 的 User-Agent 必须是注册返回的完整 `tns_marker`**，只用 bot 名会被 WAF 打回 403。
2. **TNS 检索省略 `radius` 时默认只搜 3″**，而 CHIME 定位误差是 709″ 量级 ⇒ 必然全漏；
   `include_frb=1` 也必填。`repeater_of` 恒为 `'null'`，**不能**用于归属判定。
3. ADS 需要个人 token；MAST 的 Pan-STARRS 锥检索走 `panstarrs/dr2/mean.csv`（**DR2** 堆叠表）。

这些核查的**结果已归档**（`data/processed/upgrade/pair565_rm_ads.json`、
`reproduction/logs/tns_search_*` 等），所以复现正文不需要重跑网络。
Their outputs ship as archives, so reproducing the text needs no network.

---

<a id="10"></a>

## 10. 故障排查 / Troubleshooting

| 症状 symptom | 原因 cause | 处置 fix |
| --- | --- | --- |
| 第 5 步 `FileNotFoundError: g4b_report.json` | 漏跑第 4 步 | 用 `run_all.sh`，或先跑 `gate_g4b.py` |
| 第 7 步 `FileNotFoundError: c2_candidates_v2.parquet` | 漏跑第 6 步 | 先跑 `c2_external_v2.py` |
| `UnicodeEncodeError: 'gbk' codec ... '\xb2'` | Windows 控制台编码 | `export PYTHONIOENCODING=utf-8`（`run_all.sh` 已设） |
| 第 18 步在复制 TNS 日志处失败 | 缺 2 个归档件 | 确认 `data/processed/upgrade/tns_search_*` 在包内 |
| `AssertionError: 冻结候选应为 28 对，实得 30` | 判据过滤被改动 | 恢复第 5 步的 `dm_window ≤ 5` 过滤；不要改正文 |
| astropy 提示 polar motion / IERS 过期 | IERS 数据只到装包那天 | 不影响本文数字（`astropy-iers-data` 已在锁里）；可离线设 `Conf.iers_auto_url` |
| 1 条测试 SKIP | 需要 190 MB 曝光图 | 属预期 / by design |
| 出图字体警告 | 用了系统字体名 | 本包只用 DejaVu；出现警告说明改过 `fig_gate`/rcParams |

---

<a id="11"></a>

## 11. 这份发布的边界 / Honest limits

1. **p 值不是测量值。** `p ≤ 0.001` 是 1000 次置换的分辨率下限；正文、产物、教程三处都这么写。
2. **23.9″ 未被分辨。** θ 只是选择器；正文措辞必须是"两个误差盒与重合相容"。
3. **|HA| 是池依赖量**，引用时必须同时给池与时间尺度。实测：5008 个有坐标行（`mjd_400`，UTC）中位 `0.0488 h`、
   `<0.1 h` 占 `78.9%`；声明域 3665（同为 UTC）`0.05 h` / `76.0%`；已命名源内部两两 `|ΔHA|` 中位 `0.0416 h`
   （124793 对、83 源）。**时间尺度对照**：同一批行改用 `t_bjd_tdb` 时 |HA| 近乎均匀（中位 `11.68 h`）
   ⇒ 目录赤经是按 `mjd_400`(UTC) 反推的，故 UTC 下 |HA|≈0 属构造结果，不是独立证据。
4. **排名只主张第一名。** 换地板口径后第 2 名可跌到第 6；第 4–7 名随局部门窗互换。
5. **DM–z 关系在低红移处系统性偏高 32–45%**（见步骤 17），所以 z 只作下沿；正文已把这层写进去。
6. **字号地板**：实测 5 张**矢量**版式图最小字号 **4.3–4.9 pt**，低于 7 pt 内部软线（RAA 位图分档
   **硬违规 0**）；位图字号无法从像素反读，故不在该项内。属稿件侧待办，不影响数字复现。
7. **期刊要件仍需作者提供**：仓库链接、AI 工具使用声明与日期、机构邮编、关键词体式。
8. **图与表的可复现性只到"重出与 PDF 逐像素相同"这一层**；期刊系统的色彩/分档要求由
   `audit_figure_specs.py` 判，它需要 `docs/figures` 目录在位。

---

<a id="12"></a>

## 12. 我怎么验证这份发布 / How this release was verified (and how to repeat it)

```bash
# 1) keep the shipped package aside, then copy the tree elsewhere / 先留随包件，再换盘换路径
cp -r frbpp_clean/reproduction /tmp/shipped_reproduction && cp -r frbpp_clean /tmp/reprotest/frbpp
cd /tmp/reprotest/frbpp
# 2) fresh environment from the lock / 全新环境 + 锁文件
python3 -m venv .venv && .venv/bin/pip install -r requirements-verified.txt && .venv/bin/pip install -e .
# 3) tests, then the whole chain / 先测试后整链
.venv/bin/python -m pytest -q && bash scripts/run_all.sh
# 4) rebuild vs the SHIPPED package, never vs the manifest it just rewrote
.venv/bin/python - <<'PY'
import json, pathlib
new, old = pathlib.Path("reproduction"), pathlib.Path("/tmp/shipped_reproduction")

def flat(d, p=""):
    out = {}
    if isinstance(d, dict):
        for k, v in d.items():
            out.update(flat(v, f"{p}/{k}"))
    elif isinstance(d, list):
        for i, v in enumerate(d):
            out.update(flat(v, f"{p}[{i}]"))
    else:
        out[p] = d
    return out

for f in sorted((new / "tables").glob("*.csv")):          # byte equality / 逐字节
    print("OK  " if f.read_bytes() == (old / f.relative_to(new)).read_bytes() else "DIFF", f)
for f in sorted((new / "json").glob("*.json")):           # leaf values / 叶子值
    a = flat(json.loads((old / f.relative_to(new)).read_text(encoding="utf-8")))
    b = flat(json.loads(f.read_text(encoding="utf-8")))
    diff = [k for k in set(a) | set(b) if a.get(k) != b.get(k)]
    print(f"{f.name:28s} leaves={len(a):4d} diff={len(diff):2d} {diff}")
PY
```
本机 2026-10-11 实测 / observed here: `pytest` 50 passed / 1 skipped / 0 failed；18 步单趟 `CHAIN_RC=0`；
三张表与两份 TNS 台账逐字节相同；794 个 JSON 叶子中只有 4 个 `built_at` 不同；6 张位图重出后与论文 PDF
内嵌图 maxdiff = 0；`audit_manuscript_numbers.py` 67/67。

**别拿新清单自比 / do not self-validate.** `repro_package.py` 会重写 `reproduction/MANIFEST.csv`，所以用
新生成的 MANIFEST 校验新文件必然全绿——那种"绿"没有信息量。比对基准必须是**随包那一份** `reproduction/`；
又因 `built_at` 每次运行都变，JSON 的正确判据是"差异叶子只能是 `built_at`"，而非"逐字节相同"。

还有一层平台差异要知道：`MANIFEST.csv` 里的 sha256 是在**产生它的那台机器**上按写出的字节算的，而
CSV/JSON 的行尾随平台而变（Windows 写 CRLF、Linux 写 LF）。所以逐字节校验只在同平台成立；跨平台的正确
判据是**值级比较**（上面的脚本按叶子值比对，并把 `built_at` 单列）。`.gitattributes` 因此对 `reproduction/**`
与 `data/raw/**` 关闭文本转换，保证克隆下来的输入与随包产物跟清单逐字节一致。

---

<a id="13"></a>

## 13. 获取与存档 / Obtaining and archiving

**Clone / 克隆**

```bash
git clone https://github.com/AcerisCentreforScience/frbpp_clean.git
cd frbpp_clean
```

**Deposit snapshot / 期刊或存档用的快照.** A deposit should freeze the source at one commit together
with the hashes of its inputs, and must exclude anything a build produces:

| 纳入 include | 排除 exclude |
| --- | --- |
| `src/`, `scripts/`, `tests/`, `config/` 下的 example 模板 | `OUT/`、`data/interim/` 里的 parquet 与 manifest |
| `data/raw/`（公开输入）与 `INPUTS.sha256` | `data/processed/upgrade/` 里由链路重算出的产物 |
| `reproduction/`（随包派生产物与 `MANIFEST.csv`）、三份文档、`pyproject.toml`、`requirements-verified.txt` | `.venv/`、`__pycache__/`、`.pytest_cache/`、任何本地凭证文件 |

`data/processed/upgrade/` 下的两个 `tns_search_*` 是在线核查的归档输出，**属于随包内容**：第 18 步要
复制它们才能完成打包，缺了就只能带着凭证重跑 TNS。`.gitignore` 已为它们设了例外。

**Self-check before publishing / 对外发布前自查**

```bash
git ls-files | wc -l                       # 期望 141
git ls-files | grep -E "local[.]json$|__pycache__|^OUT/|interim/"   # 期望只有 data/interim/.gitkeep
git ls-files | grep -c "tns_search_"               # 期望 3：两个随包归档件 + 产生它们的脚本
grep -rn "G:/FRB" --include="*.py" --include="*.sh" . ; echo rc=$?  # 期望 rc=1（无残留绝对路径）
sha256sum -c INPUTS.sha256                 # 期望 10 个输入文件全部 OK（Windows 用 Get-FileHash 比对）
```

最后一条最关键：哈希不符说明手里的目录数据不是产出这些结果的那一版发布，应当停止复现并重新取数
（见 [§3](#3)），而不是继续把差异解释掉。

---

<a id="14"></a>

## 14. 术语对照 / Glossary

| 符号 symbol | English | 中文 |
| --- | --- | --- |
| declaration domain | events with finite DM and coordinates, 3665 | 有限 DM 且有坐标的 3665 个事件 |
| apparent one-off | event not grouped by `repeater_name` | 未被 `repeater_name` 归组的表观单次事件 |
| event-layer de-duplication (v1.3) | collapse Cat1 rows that re-measure a Cat2 event | 把重测同一暴的 Cat1 行折回事件层 |
| criterion set / union | 28 pairs at DM ≤ 5 / 30 pairs over the nine grid cells | 判据池 28 对 / 九档并集 30 对（更严 null bar） |
| DM error floor | `s(DM) = 1.4826·MAD(dm_fitb − bonsai_dm)` in a local DM window | 局部 DM 窗内两口径差的稳健尺度 |
| σ_comb | √2·s (conservative) or s (loose) | 两暴合成误差（保守/宽松并列） |
| r_floorB / r_floorA | \|ΔDM\|/σ_comb with either convention | 两约定下的地板一致度 |
| θ | great-circle angle between the two localisations | 两暴大圆角距（**未被分辨**） |
| \|HA\| | absolute hour-angle difference | 中天角之差（漂移扫描的量化痕迹） |
| sidereal-day quantization | Δt ≈ integer × T_sid | 恒星日量化（漂扫的结构性后果） |
| null battery | four permutation families, B = 1000 | 四类置换零假设电池 |
| resolution limit | p ≤ 1/(B+1) | 分辨率下限（不是测量值） |
| union-find component | one independent physical hypothesis | 一个独立物理假设 |
| positive control | known repeaters replayed under the criterion | 已知重复源在同判据下的召回 |
| retext branch | erase-and-redraw figure path, Windows fonts | 换字支线（需 Windows 字体，不在复现路径） |

---

<a id="15"></a>

## 15. 常见问题 / FAQ

**Q: 为什么我看到 28 对，又看到 30 对？**
A: 28 是判据池（DM 窗 ≤ 5），30 是九档并集，只作更严的 null bar。第 5 步把前者写进
`c2_candidates_v1.csv`，把只在 DM±20 出现的两对写进 `c2_candidates_edge_v1.csv`。
**EN.** 28 is the criterion set, 30 the stricter union used only as the null bar; step 5 keeps them
in separate files, and step 18 asserts the 28.

**Q: 为什么不直接报一个 p 值？** A: B = 1000 下 0/1000 只能给 `p ≤ 1/(B+1)`；把它写成测量值是错的。
**Q: 为什么禁用 "periodicity"？** A: 565 d 恰为 567.002 个恒星日，是漂扫的算术后果，不是源的周期。
**Q: 目录里没找到第三个事件 = 该源没爆发？** A: 不等于。只能写"在我们可读记录内为空"。
**Q: 换了 Cat2 版本还能用这套数吗？** A: 不能。地板统计是这份 CSV 自测出来的，必须整链重算。
**Q: θ = 23.9″ 算证据吗？** A: 不算，两暴误差盒在任何 ρ 下都远大于 θ，它只是选择器。

---

<a id="16"></a>

## 16. 设计与可复现性要点 / Design notes

These are the mechanisms that make the numbers reproducible rather than incidental, plus the two
points where a reader has to be careful about interpretation.

| 要点 | 说明 |
| --- | --- |
| 判据守卫 | 第 5 步只保留 DM 窗 ≤ 5 的档位配置（即论文主张的 28 对判据池），仅由最松档 DM±20 支持的配对另存 `c2_candidates_edge_v1.csv`；更严的 null bar 由第 10 步自己从九档网格重建，两者不混用 |
| 打包断言 | 第 18 步对 28 对、8 分量、保留 7、剔除 21、主候选地板比 0.1191 逐条硬断言，任一条不符即中止，不产出"差不多"的包 |
| 误差地板 | ΔDM 一致度的分母是目录自测的经验地板 1.4826·MAD(D)，形式误差只作并列对照，永不作分母 |
| p 值语义 | B = 1000 时 0/1000 只能给 `p <= 1/(B+1)`；这是分辨率下限，不是测量值 |
| 池与时间尺度 | 任何 HA 差值必须同时给事件池与时间尺度：5008 行（UTC）与声明域 3665 给出的比例不同，改用 TDB 则分布近乎均匀 |
| 角距不是证据 | 23.9 角秒在任何相关系数 ρ 下都远小于合成 1σ 误差，只作选择器使用，正文措辞是两个误差盒与重合相容 |
| 关系标定域 | 第 17 步的 DM–z 关系在 z < 0.7 的已定位样本上高估 32–45%，故 z 只作"色散量是宇宙学的"下沿指示；脚本把朴素前因子与所用前因子并排打印，12% 的差异可见而非隐藏 |
| 可移植性 | 所有脚本以脚本自身位置定位仓库根目录，无绝对路径、无当前工作目录依赖；`run_all.sh` 强制 UTF-8 输出 |
| 依赖锁定 | `requirements-verified.txt` 锁到产出结果的那一套版本，`pyproject.toml` 只给下限 |

Two properties of the underlying catalogue constrain reading rather than computation, and are worth
repeating: hour angles cluster near transit because CHIME drift-scans, so an interval close to an
integer number of sidereal days is a property of the instrument, not of the source; and per-burst
localisation boxes overlap far more widely than the 23.9 arcsecond separation, so a small angle
selects candidates rather than demonstrates an association.
