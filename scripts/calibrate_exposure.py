"""Build the detection-efficiency curve v1 from injections.

Usage: python scripts/calibrate_exposure.py
Output: data/processed/exposure/injection_efficiency_v1.json
"""
from __future__ import annotations

import datetime as dt
import hashlib
import json
import subprocess
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from frbpp.catalogs import chime          # noqa: E402
from frbpp.exposure import calibrate      # noqa: E402

RAW_INJ = ROOT / "data/raw/chime/06_injections/chimefrb_catalog1_injections.p"
OUT_DIR = ROOT / "data/processed/exposure"


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


def _print_curve(name: str, cur) -> None:
    print(f"— {name}")
    for _, r in cur.iterrows():
        print(f"    [{r['bin_lo']:>10.4g}, {r['bin_hi']:>10.4g}]  n={int(r['n']):>6}  "
              f"eff={r['eff']:.3f}  [{r['eff_lo']:.3f}, {r['eff_hi']:.3f}]")


def main() -> int:
    if not RAW_INJ.exists():
        raise SystemExit(f"缺少注入数据: {RAW_INJ}（先运行 scripts/fetch_chime.py）")
    inj = chime.load_injections(RAW_INJ)
    det = calibrate.detection_mask(inj)
    print(f"注入事件: {len(inj)}; 检出: {int(det.sum())}（{det.mean():.3%}）")
    print(f"判据: SNR>{calibrate.SNR_THRESHOLD} 且 (RFI等级>{calibrate.RFI_THRESHOLD} "
          f"或 SNR>{calibrate.HIGH_SNR_OVERRIDE}) —— Catalog 1 教程口径")

    cur_fluence = calibrate.efficiency_curve(inj, "fluence_jy_ms_inj", bins=8)
    cur_dm = calibrate.efficiency_curve(inj, "dm_inj", bins=8)
    _print_curve("fluence_jy_ms_inj", cur_fluence)
    _print_curve("dm_inj", cur_dm)

    from scipy.stats import spearmanr
    rho_f, _ = spearmanr(cur_fluence["bin_lo"], cur_fluence["eff"])
    rho_d, _ = spearmanr(cur_dm["bin_lo"], cur_dm["eff"])
    print(f"单调性（Spearman ρ）: fluence={rho_f:.3f}, dm={rho_d:.3f}")
    assert rho_f > 0.9, "fluence 效率曲线非单调, 检查判据或数据"
    assert cur_fluence["eff"].iloc[-1] > cur_fluence["eff"].iloc[0]

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    out_path = OUT_DIR / "injection_efficiency_v1.json"
    payload = {
        "built_at": dt.datetime.now().astimezone().isoformat(timespec="seconds"),
        "git": git_sha(),
        "input": {"path": str(RAW_INJ.relative_to(ROOT)), "sha256": sha256(RAW_INJ)},
        "criterion": {
            "snr_threshold": calibrate.SNR_THRESHOLD,
            "rfi_threshold": calibrate.RFI_THRESHOLD,
            "high_snr_override": calibrate.HIGH_SNR_OVERRIDE,
            "note": "Catalog 1 教程口径; 条件效率 P(检出|已送上注入)（NOTES F8）",
        },
        "summary": {"n_injected": int(len(inj)), "n_detected": int(det.sum()),
                    "overall_eff": float(det.mean())},
        "curves": {
            "fluence_jy_ms_inj": cur_fluence.to_dict(orient="records"),
            "dm_inj": cur_dm.to_dict(orient="records"),
        },
    }
    out_path.write_text(json.dumps(payload, indent=2, ensure_ascii=False),
                        encoding="utf-8")
    print(f"已写入: {out_path.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())