"""Normalise raw catalogue files into data/interim Parquet plus a provenance manifest.

Usage: python scripts/build_interim.py
Output: data/interim/{cat1_bursts,cat2_bursts,r3_180916_toa}.parquet, data/interim/_build_manifest.json
"""
from __future__ import annotations

import datetime as dt
import hashlib
import json
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from frbpp.catalogs import chime  # noqa: E402

RAW = ROOT / "data" / "raw" / "chime"
INTERIM = ROOT / "data" / "interim"


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def git_sha() -> str:
    try:
        out = subprocess.run(["git", "-C", str(ROOT), "rev-parse", "HEAD"],
                             capture_output=True, text=True, timeout=10)
        return out.stdout.strip() if out.returncode == 0 else "no-git"
    except Exception:  # noqa: BLE001
        return "no-git"


def main() -> int:
    INTERIM.mkdir(parents=True, exist_ok=True)
    srcs = {
        "cat1": RAW / "05_cat1" / "cat1_table2_vizier.tsv",
        "cat2": RAW / "cat2_25" / "chimefrbcat2.csv",
        "r3_toa": RAW / "03_180916" / "R3_CHIME_TOA_20Feb4.txt",
        "injections": RAW / "06_injections" / "chimefrb_catalog1_injections.p",
    }
    for k, v in srcs.items():
        if not v.exists():
            raise SystemExit(f"缺少原始文件 {k}: {v}（先运行 scripts/fetch_chime.py）")

    df1, cat1 = chime.load_cat1(srcs["cat1"])
    p1 = INTERIM / "cat1_bursts.parquet"
    df1.to_parquet(p1, index=False)
    print(f"cat1: {len(df1)} 行 → {p1.name} | BJD_TDB {cat1.t.min():.3f}–{cat1.t.max():.3f}"
          f" | dm_err NaN {int(df1['dm_fitb_err'].isna().sum())}")

    df2, cat2 = chime.load_cat2(srcs["cat2"])
    p2 = INTERIM / "cat2_bursts.parquet"
    df2.to_parquet(p2, index=False)
    print(f"cat2: {len(df2)} 行 → {p2.name} | BJD_TDB {cat2.t.min():.3f}–{cat2.t.max():.3f}")

    dft, _ = chime.load_r3_180916_toa(srcs["r3_toa"])
    p3 = INTERIM / "r3_180916_toa.parquet"
    dft.to_parquet(p3, index=False)
    span = dft["mjd_native"].max() - dft["mjd_native"].min()
    print(f"r3_toa: {len(dft)} 行 → {p3.name} | 跨度 {span:.4f} d")

    inj = chime.load_injections(srcs["injections"])
    print(f"injections(检查): {inj.shape[0]} 行 × {inj.shape[1]} 列（阶段 2 用, 暂不落 parquet）")

    manifest = {
        "built_at": dt.datetime.now().astimezone().isoformat(timespec="seconds"),
        "git": git_sha(),
        "inputs": {k: {"path": str(v.relative_to(ROOT)), "sha256": sha256(v)}
                   for k, v in srcs.items()},
        "outputs": {p.name: {"rows": n, "sha256": sha256(p)}
                    for p, n in ((p1, len(df1)), (p2, len(df2)), (p3, len(dft)))},
        "notes": "见 data/raw/chime/SOURCES.json 与 docs/NOTES.md（数据源实测记录）",
    }
    (INTERIM / "_build_manifest.json").write_text(
        json.dumps(manifest, indent=2, ensure_ascii=False), encoding="utf-8")
    print("manifest:", INTERIM / "_build_manifest.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())