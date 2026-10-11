r"""Invert the DM-z relation quoted in the text for the pair's DM budget.

The manuscript states DM_IGM(z) = 1000 * integral_0^z (1+z')/E(z') dz' with Planck 2018
cosmology (macquart2020; scatter measured by baptista2024), so that prefactor is what the
inversion uses. The prefactor derived from the same cosmology's baryon density,
3 c H0 Omega_b / (8 pi G m_p), is reported next to it because the two differ by 12% and the
difference moves the inferred redshift by 0.17.

Usage: python scripts/dm_z_inversion.py
Output: data/processed/upgrade/dm_z_inversion.json
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
from scipy.integrate import quad

ROOT = Path(__file__).resolve().parents[1]

H0, OM, OB = 67.4, 0.315, 0.0493
G, MP, C, PC = 6.6743e-8, 1.6726219e-24, 2.99792458e10, 3.085677581e18

K_TEXT = 1000.0
DM_OBS = 1834.95
MW = {"NE2001": 119.11, "YMW16": 128.56, "mean": 123.8}
HALO = {"30": 30.0, "40": 40.0, "60": 60.0}
DM_FORMAL_ERR = 0.5
SCATTER = [0.20, 0.30]
ANCHORS = [("FRB181112", 0.189, 108.0), ("FRB180924", 0.321, 211.0),
           ("FRB190520B", 0.660, 496.0)]


def derived_k(h0: float, om: float, ob: float) -> float:
    h0_s = h0 * 1e5 / (PC * 1e6)
    n_e0 = 3.0 * h0_s ** 2 * ob / (8.0 * np.pi * G * MP)
    return n_e0 * (C / h0_s) / PC


def dm_igm(z: float, k: float = K_TEXT, om: float = OM, ol: float | None = None) -> float:
    ol = (1.0 - om) if ol is None else ol
    v, _ = quad(lambda x: (1.0 + x) / np.sqrt(om * (1.0 + x) ** 3 + ol), 0.0, z)
    return k * v


def z_of_dm(dm: float, **kw) -> float:
    lo, hi = 0.0, 12.0
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if dm_igm(mid, **kw) < dm:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


K_DERIVED = derived_k(H0, OM, OB)
print(f"正文前因子 K = {K_TEXT:.1f}   同宇宙学推导 K = {K_DERIVED:.1f} "
      f"(差 {100 * (K_DERIVED / K_TEXT - 1):+.1f}%)")
print(f"关系检验: DM_IGM(z=1) = {dm_igm(1.0):.0f} (K=1000) / "
      f"{dm_igm(1.0, k=K_DERIVED):.0f} (K=推导值) pc cm^-3")

print("\n== 低红移锚点（同一关系外推到 z<0.7 的一致性）==")
anchors = []
for nm, z, obs in ANCHORS:
    pred = dm_igm(z)
    dev = 100.0 * (obs - pred) / pred
    anchors.append({"name": nm, "z": z, "published_dm_igm": obs,
                    "relation_dm_igm": round(pred, 1), "dev_pct": round(dev, 1)})
    print(f"  {nm:12s} z={z:.3f}  关系预测={pred:7.1f}  发表值={obs:6.1f}  偏差={dev:+6.1f}%")
mean_abs_dev = float(np.mean([abs(x["dev_pct"]) for x in anchors]))
print(f"  平均绝对偏差 = {mean_abs_dev:.1f}%  → 该同族关系在 z<0.7 系统性高估 DM_IGM，"
      f"本反演只用于 z~1-2 的高红移预算")

print("\n== 反解 z（系统项组合，K=1000）==")
combos = []
base = None
for hk, hv in HALO.items():
    for mk, mv in MW.items():
        dig = DM_OBS - mv - hv
        zz = z_of_dm(dig)
        combos.append({"halo": hv, "galactic": mv, "dm_igm": round(dig, 2), "z": round(zz, 3)})
        if hk == "40" and mk == "mean":
            base = (dig, zz)
        print(f"  晕={hk:>2s} 星系={mk:8s} -> DM_IGM={dig:7.1f}  z={zz:.3f}")
dig0, z0 = base
print(f"\n采用（晕 40、星系 123.8）: DM_IGM={dig0:.1f}  **z={z0:.3f}**")
print(f"  同一预算在 K={K_DERIVED:.1f} 下给 z={z_of_dm(dig0, k=K_DERIVED):.3f}")

bands = {}
print("\n== 形式误差 vs 固有散布（加在 DM_IGM 上）==")
lo = z_of_dm(dig0 - DM_FORMAL_ERR)
hi = z_of_dm(dig0 + DM_FORMAL_ERR)
bands["formal_dm_error_pm"] = {"dm": DM_FORMAL_ERR, "z": [round(lo, 3), round(hi, 3)]}
print(f"  形式 ±{DM_FORMAL_ERR} dmpc        -> z in [{lo:.3f}, {hi:.3f}]")
for frac in SCATTER:
    lo = z_of_dm(dig0 * (1.0 - frac))
    hi = z_of_dm(dig0 * (1.0 + frac))
    bands[f"intrinsic_scatter_pm{int(frac * 100)}"] = {
        "dm_igm": [round(dig0 * (1 - frac), 1), round(dig0 * (1 + frac), 1)],
        "z": [round(lo, 2), round(hi, 2)]}
    print(f"  ±{frac:.0%} DM_IGM 散布      -> z in [{lo:.2f}, {hi:.2f}]")

print("\n== 宇宙学敏感性（DM_IGM 固定）==")
cosmo = []
for h0p, omp, obp, lab in [(70.0, 0.30, 0.048, "SH0R18 上界"),
                           (67.4, 0.315, 0.0493, "Planck 2018"),
                           (65.0, 0.32, 0.050, "低 H0")]:
    zz = z_of_dm(dig0, k=derived_k(h0p, omp, obp))
    cosmo.append({"label": lab, "h0": h0p, "omega_m": omp, "omega_b": obp,
                  "k_derived": round(derived_k(h0p, omp, obp), 1), "z": round(zz, 3)})
    print(f"  {lab:14s} H0={h0p} Om={omp} Ob={obp} -> z={zz:.3f}（用推导 K）")
for lab, k in [("正文 K=1000", K_TEXT), ("Planck 推导 K", K_DERIVED)]:
    print(f"  {lab:14s} -> z={z_of_dm(dig0, k=k):.3f}")

out = {
    "dm_obs": DM_OBS,
    "mw": MW,
    "halo_tested": list(HALO.values()),
    "relation": {"form": "DM_IGM(z) = K * integral_0^z (1+z')/E(z') dz'",
                 "k_used": K_TEXT, "k_derived_from_cosmology": round(K_DERIVED, 1),
                 "cosmology": {"h0": H0, "omega_m": OM, "omega_b": OB},
                 "cited_as": "macquart2020; scatter baptista2024"},
    "adopted": {"halo": 40.0, "galactic": MW["mean"], "dm_igm": round(dig0, 2),
                "z": round(z0, 3)},
    "z_with_derived_k": round(z_of_dm(dig0, k=K_DERIVED), 3),
    "bands": bands,
    "systematics_grid": combos,
    "cosmology_sensitivity": cosmo,
    "low_z_anchors": anchors,
    "low_z_anchor_mean_abs_dev_pct": round(mean_abs_dev, 1),
    "dm_igm_z1": {"k_used": round(dm_igm(1.0), 1), "k_derived": round(dm_igm(1.0, k=K_DERIVED), 1)},
}
(ROOT / "data/processed/upgrade/dm_z_inversion.json").write_text(
    json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")
print("\n写入 data/processed/upgrade/dm_z_inversion.json")
