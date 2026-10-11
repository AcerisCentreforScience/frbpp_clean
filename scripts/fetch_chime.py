"""Download the public CHIME/FRB releases into data/raw.

Usage: python scripts/fetch_chime.py
"""
from __future__ import annotations

import datetime as dt
import hashlib
import json
import re
import sys
import time
from pathlib import Path

import requests

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw" / "chime"
BASE = "https://www.canfar.net/storage/vault"
HEADERS = {"User-Agent": "frbpp/0.1 (CHIME open data fetch; stage-1)"}

RELEASES = {
    "05_cat1": "AstroDataCitationDOI/CISTI.CANFAR/21.0007/data",
    "03_180916": "AstroDataCitationDOI/CISTI.CANFAR/20.0002/data",
    "06_injections": "AstroDataCitationDOI/CISTI.CANFAR/22.0005/data",
    "repeaters_23": "AstroDataCitationDOI/CISTI.CANFAR/23.0004/data",
    "cat2_25": "AstroDataCitationDOI/CISTI.CANFAR/25.0066/data",
}

DIRECT_URLS = {
    ("05_cat1", "cat1_table2_vizier.tsv"):
        "https://vizier.cds.unistra.fr/viz-bin/asu-tsv"
        "?-source=J/ApJS/257/59/table2&-out.max=unlimited&-out.all",
}

_ROW_RE = re.compile(r"rows\.push\((\[.*?\])\)", re.S)


def _get(url: str, *, stream: bool = False, tries: int = 4, timeout: int = 120):
    last = None
    for k in range(tries):
        try:
            r = requests.get(url, headers=HEADERS, stream=stream, timeout=timeout)
            if r.status_code == 200:
                return r
            r.close()
            last = f"HTTP {r.status_code}"
        except Exception as exc:  # noqa: BLE001
            last = f"{type(exc).__name__}: {exc}"
        time.sleep(1.5 * (k + 1))
    raise RuntimeError(f"GET 失败 {url}: {last}")


def list_dir(path: str) -> list[dict]:
    html = _get(f"{BASE}/list/{path}").text
    out = []
    for m in _ROW_RE.finditer(html):
        try:
            arr = json.loads(m.group(1))
        except json.JSONDecodeError:
            continue
        if len(arr) < 10:
            continue
        name, size, icon, full_path = arr[1], arr[2], arr[8], arr[9]
        out.append({
            "name": name,
            "size": size,
            "type": "dir" if "folder" in str(icon) else "file",
            "path": full_path,
        })
    return out


def probe(targets: list[str] | None = None) -> int:
    items = ([(t, RELEASES.get(t, t)) for t in targets] if targets
             else list(RELEASES.items()))
    for key, path in items:
        print(f"===== {key} :: {path}")
        try:
            rows = list_dir(path)
        except Exception as exc:  # noqa: BLE001
            print(f"  !! {type(exc).__name__}: {exc}")
            continue
        for row in rows:
            print(f"  [{row['type']:4}] {row['size']:>12}  {row['name']}")
    return 0


def _download_with_resume(url: str, dest: Path, max_attempts: int = 60):
    dest.parent.mkdir(parents=True, exist_ok=True)
    pos = dest.stat().st_size if dest.exists() else 0
    h = hashlib.sha256()
    if pos:
        with open(dest, "rb") as fh:
            for chunk in iter(lambda: fh.read(1 << 20), b""):
                h.update(chunk)
        print(f"    resume from {pos} bytes")
    total = None
    for attempt in range(1, max_attempts + 1):
        headers = dict(HEADERS)
        if pos:
            headers["Range"] = f"bytes={pos}-"
        try:
            r = requests.get(url, headers=headers, stream=True, timeout=120)
            if r.status_code not in (200, 206):
                r.close()
                raise RuntimeError(f"HTTP {r.status_code}")
            if pos and r.status_code == 200:
                pos, h = 0, hashlib.sha256()
                print("    (server ignored Range; restarting)")
            cr = r.headers.get("Content-Range")
            if cr and "/" in cr:
                total = int(cr.rsplit("/", 1)[1])
            elif r.status_code == 200 and r.headers.get("Content-Length"):
                total = int(r.headers["Content-Length"])
            with open(dest, "ab" if pos else "wb") as fh:
                for chunk in r.iter_content(chunk_size=1 << 20):
                    if chunk:
                        fh.write(chunk)
                        h.update(chunk)
                        pos += len(chunk)
            r.close()
            if total is None or pos >= total:
                return pos, h.hexdigest()
            print(f"    incomplete {pos}/{total}; retrying")
        except Exception as exc:  # noqa: BLE001
            print(f"    interrupted at {pos} bytes: {type(exc).__name__}; "
                  f"retry {attempt}/{max_attempts}")
        time.sleep(min(5 + attempt, 30))
    raise RuntimeError(f"下载失败（{max_attempts} 次尝试）: {url}")


def fetch(release: str, files: list[str]) -> int:
    if release not in RELEASES:
        raise SystemExit(f"未知 release: {release}（可选 {sorted(RELEASES)}）")
    dest_dir = RAW / release
    dest_dir.mkdir(parents=True, exist_ok=True)
    manifest_path = RAW / "SOURCES.json"
    manifest = {}
    if manifest_path.exists():
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    base_path = RELEASES[release]
    for fname in files:
        rel_path = f"{base_path}/{fname}"
        direct = DIRECT_URLS.get((release, fname))
        url = direct or f"{BASE}/file/{rel_path}"
        dest = dest_dir / Path(fname).name
        print(f"--> {url}")
        size, digest = _download_with_resume(url, dest)
        entry = {
            "url": url,
            "source": "vizier" if direct else "canfar",
            "release": release,
            "file": fname,
            "bytes": size,
            "sha256": digest,
            "retrieved_at": dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds"),
        }
        manifest[f"{release}::{fname}" if direct else rel_path] = entry
        print(f"    saved {dest} ({size} bytes) sha256={digest[:16]}...")
    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    manifest_path.write_text(json.dumps(manifest, indent=2, ensure_ascii=False),
                             encoding="utf-8")
    return 0


def main(argv: list[str]) -> int:
    cmd = argv[1] if len(argv) > 1 else "probe"
    if cmd == "probe":
        return probe(argv[2:] or None)
    if cmd == "fetch":
        if len(argv) < 4:
            raise SystemExit("用法: fetch <release> <file> [<file> ...]")
        return fetch(argv[2], argv[3:])
    raise SystemExit(f"未知命令: {cmd}（支持: probe / fetch）")


if __name__ == "__main__":
    sys.exit(main(sys.argv))