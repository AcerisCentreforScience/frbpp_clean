"""Detection efficiency epsilon estimated from CHIME injections by injection-recovery."""
from __future__ import annotations

import numpy as np
import pandas as pd

SNR_THRESHOLD = 9.0
RFI_THRESHOLD = 5.0
HIGH_SNR_OVERRIDE = 100.0


def detection_mask(inj: pd.DataFrame) -> np.ndarray:
    snr = inj["bonsai_snr"].to_numpy(dtype=float)
    grade = inj["l2_rfi_grade"].to_numpy(dtype=float)
    snr_ok = np.isfinite(snr) & (snr > SNR_THRESHOLD)
    rfi_ok = (np.isfinite(grade) & (grade > RFI_THRESHOLD)) | (snr > HIGH_SNR_OVERRIDE)
    return snr_ok & rfi_ok


def wilson_interval(k, n, z: float = 1.959963984540054):
    k = np.asarray(k, dtype=float)
    n = np.asarray(n, dtype=float)
    lo = np.zeros_like(n)
    hi = np.zeros_like(n)
    m = n > 0
    if np.any(m):
        kk, nn = k[m], n[m]
        p = kk / nn
        denom = 1.0 + z * z / nn
        centre = (p + z * z / (2.0 * nn)) / denom
        half = z * np.sqrt(p * (1.0 - p) / nn + z * z / (4.0 * nn * nn)) / denom
        lo[m] = np.clip(centre - half, 0.0, 1.0)
        hi[m] = np.clip(centre + half, 0.0, 1.0)
    return lo, hi


def efficiency_curve(inj: pd.DataFrame, value_col: str, bins,
                     detected=None) -> pd.DataFrame:
    detected = detection_mask(inj) if detected is None else np.asarray(detected, bool)
    v = inj[value_col].to_numpy(dtype=float)
    ok = np.isfinite(v)
    df = pd.DataFrame({"v": v[ok], "det": detected[ok]})
    if df.empty:
        raise ValueError(f"{value_col} 全为空, 无法分箱")
    if np.isscalar(bins):
        edges = np.unique(np.quantile(df["v"], np.linspace(0.0, 1.0, int(bins) + 1)))
        if edges.size < 3:
            raise ValueError("分位箱退化（值分布过于集中）")
    else:
        edges = np.asarray(bins, dtype=float)
    cut = pd.cut(df["v"], bins=edges, include_lowest=True)
    grp = df.groupby(cut, observed=True)["det"]
    sizes = grp.size()
    n = sizes.to_numpy(dtype=float)
    k = grp.sum().to_numpy(dtype=float)
    lo, hi = wilson_interval(k, n)
    eff = np.divide(k, n, out=np.zeros_like(k), where=n > 0)
    intervals = sizes.index
    return pd.DataFrame({
        "bin_lo": [iv.left for iv in intervals],
        "bin_hi": [iv.right for iv in intervals],
        "n": n.astype(int),
        "n_det": k.astype(int),
        "eff": eff,
        "eff_lo": lo,
        "eff_hi": hi,
    })