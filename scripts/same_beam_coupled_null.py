r"""Coupling-preserving null for the pair-level joint box: permute dispersion measures only.

CHIME localisations are tied to transit time - catalogue right ascension is derived from the
arrival time - so a null that jitters positions independently, or randomises right ascension,
dilutes that coupling and can understate the chance rate of small separations. Here coordinates,
epochs and hour angles are left untouched and only the dispersion measures are permuted among
declaration-domain events, so the drift-scan coupling is preserved by construction.

A deliberately broken companion null (right ascension redrawn uniformly) is reported alongside, so
the price of the coupling assumption is measured rather than argued.

Pool definitions, cuts, the separation and the hour-angle routine are imported from
same_beam_test.py so the two audits cannot drift apart. Note that sep_arcsec takes RADIANS and
hour_angle_h takes DEGREES, as upstream does.

Usage: python scripts/same_beam_coupled_null.py --draws 20000 --seed 20261011
Output: data/processed/upgrade/same_beam_coupled_null.json
"""
from __future__ import annotations

import argparse
import datetime as dt
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

from frbpp.upgrade.crossmatch import build_event_table  # noqa: E402
from same_beam_test import (CAT1, CAT2, DDM_CUT, HA_SAME_BEAM_H, THETA_CUT,  # noqa: E402
                            hour_angle_h, sep_arcsec)

OUT = ROOT / "data/processed/upgrade/same_beam_coupled_null.json"
DECLARATION_DOMAIN = 3665


def box_counts(a_idx, b_idx, dm_vec, rng, draws: int) -> np.ndarray:
    """Number of pairs inside the dispersion-measure box, once per permutation of dm_vec."""
    out = np.empty(draws, dtype=np.int64)
    for d in range(draws):
        s = dm_vec[rng.permutation(dm_vec.size)]
        out[d] = int((np.abs(s[a_idx] - s[b_idx]) < DDM_CUT).sum())
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--draws", type=int, default=20000)
    ap.add_argument("--seed", type=int, default=20261011)
    a = ap.parse_args()
    rng = np.random.default_rng(a.seed)

    ev = build_event_table(pd.read_parquet(CAT1), pd.read_parquet(CAT2), event_dedup=True)
    mjd = np.asarray(ev.mjd, float)
    ra = np.radians(np.asarray(ev.ra, float))
    dec = np.radians(np.asarray(ev.dec, float))
    dm = np.asarray(ev.dm, float)
    usable = np.array([r == "" for r in ev.src_rep]) & np.isfinite(dm) & np.isfinite(ra)
    idx = np.where(usable)[0]
    assert idx.size == DECLARATION_DOMAIN, f"声明域应为 {DECLARATION_DOMAIN}，实得 {idx.size}"
    ra_d, dec_d, dm_d, t_d = ra[idx], dec[idx], dm[idx], mjd[idx]
    n = ra_d.size
    print(f"声明域事件 n = {n}（表观一次性、有限 DM 与坐标；坐标/历元/时角全部保持原样）")

    ii, jj = np.triu_indices(n, 1)
    th = sep_arcsec(ra_d[ii], dec_d[ii], ra_d[jj], dec_d[jj])
    small = th < THETA_CUT
    ia, ib = ii[small], jj[small]
    ddt = np.abs(t_d[ia] - t_d[ib])
    ha = hour_angle_h(t_d, np.degrees(ra_d))
    dha = np.abs(ha[ia] - ha[ib])
    dha = np.abs((dha + 12.0) % 24.0 - 12.0)
    print(f"θ < {THETA_CUT}″ 的对 = {len(ia)}（全池 {len(th):,} 对）；其中跨日 {int((ddt > 1.0).sum())}")

    pools = {
        "all pairs with theta < 24": np.ones(len(ia), bool),
        "cross-day (|dt| > 1 d) with theta < 24": ddt > 1.0,
        f"cross-day AND same-beam proxy (|dHA| < {HA_SAME_BEAM_H} h)":
            (ddt > 1.0) & (dha < HA_SAME_BEAM_H),
    }

    out = {"built_at": dt.datetime.now().astimezone().isoformat(timespec="seconds"),
           "null": "permute dispersion measures only; coordinates, epochs and hour angles untouched",
           "cuts": {"theta_arcsec": THETA_CUT, "ddm_pc_cm3": DDM_CUT,
                    "ha_same_beam_h": HA_SAME_BEAM_H},
           "pool_events": int(n), "pairs_theta_lt_cut": int(len(ia)),
           "draws": int(a.draws), "seed": int(a.seed), "pools": {}}

    for label, m in pools.items():
        x, y = ia[m], ib[m]
        if x.size == 0:
            out["pools"][label] = {"pairs_in_pool": 0}
            print(f"{label:52s} 无对可测")
            continue
        observed = int((np.abs(dm_d[x] - dm_d[y]) < DDM_CUT).sum())
        cnt = box_counts(x, y, dm_d, rng, a.draws)

        ra_rand = rng.uniform(0.0, 2.0 * np.pi, size=n)
        th_r = sep_arcsec(ra_rand[ii], dec_d[ii], ra_rand[jj], dec_d[jj])
        cr = (th_r < THETA_CUT) & (np.abs(t_d[ii] - t_d[jj]) > 1.0)
        cnt_broken = box_counts(ii[cr], jj[cr], dm_d, rng, a.draws) if int(cr.sum()) else np.zeros(1)

        out["pools"][label] = {
            "pairs_in_pool": int(x.size),
            "observed_in_box": observed,
            "lambda_mean": round(float(cnt.mean()), 4),
            "lambda_sd": round(float(cnt.std()), 4),
            "lambda_ci95": [round(float(np.percentile(cnt, 2.5)), 4),
                            round(float(np.percentile(cnt, 97.5)), 4)],
            "p_at_least_one": round(float((cnt >= 1).mean()), 5),
            "p_ge_observed": round(float((cnt >= observed).mean()), 5),
            "broken_coupling_pairs_in_pool": int(cr.sum()),
            "broken_coupling_lambda_mean": round(float(cnt_broken.mean()), 4),
            "broken_coupling_lambda_ci95": [round(float(np.percentile(cnt_broken, 2.5)), 4),
                                            round(float(np.percentile(cnt_broken, 97.5)), 4)],
        }
        print(f"{label:52s} 对内 {x.size:7,} 观测入盒 {observed} | 保耦合 λ={cnt.mean():.4f}"
              f"±{cnt.std():.4f} P(≥1)={(cnt >= 1).mean():.4f} P(≥观测)="
              f"{float((cnt >= observed).mean()):.4f} | 打散RA λ={cnt_broken.mean():.4f}")

    out["reading"] = ("lambda is the chance count of joint-box pairs under a null that keeps the "
                      "drift-scan RA-time coupling intact; p_ge_observed is the tail probability at or "
                      "above the observed count. The broken-coupling column measures how much that "
                      "assumption is worth and is a diagnostic, not a competing result.")
    OUT.write_text(json.dumps(out, indent=2, ensure_ascii=False), encoding="utf-8")
    print("写入", OUT.relative_to(ROOT))
    return 0


if __name__ == "__main__":
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    raise SystemExit(main())
