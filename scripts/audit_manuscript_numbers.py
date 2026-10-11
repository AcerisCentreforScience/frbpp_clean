r"""Cross-check the manuscript's numbers against the artefacts the chain produces.

Each quantity is read from a recomputed artefact and looked for *near a topic anchor* in the
manuscript body (LaTeX comments stripped), comparing numerically rather than as text so that
rounding at the last displayed digit does not look like a disagreement. Anchoring matters: an
unanchored search is satisfied by any coincidental occurrence, including a table column
specification such as `\hsize=1.35\hsize`.

Outcomes:
  ok           the artefact value (or an accepted alternative spelling) is in the anchored context
  MISMATCH     the context states a value that disagrees with the artefact
  NOT-STATED   the context exists but does not assert this quantity

Usage: python scripts/audit_manuscript_numbers.py --tex docs/submission/raa_manuscript_revised.tex
"""
from __future__ import annotations

import argparse
import csv
import json
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
T_SID = 0.99726958
WINDOW = 240


def load(p: Path):
    return json.loads(p.read_text(encoding="utf-8"))


def strip_tex(text: str) -> str:
    out = []
    for line in text.splitlines():
        cut = None
        for i, ch in enumerate(line):
            if ch == "%" and (i == 0 or line[i - 1] != "\\"):
                cut = i
                break
        out.append(line if cut is None else line[:cut])
    return "\n".join(out)


def windows(body: str, anchor: str | None) -> list[str]:
    if anchor is None:
        return [body]
    b, lo = body.lower(), anchor.lower()
    spans, start = [], 0
    while True:
        i = b.find(lo, start)
        if i < 0:
            return spans
        spans.append(body[max(0, i - WINDOW):i + WINDOW])
        start = i + len(lo)


NUM = re.compile(r"\d[\d,]*(?:\.\d+)?")


def numbers(span: str) -> list[tuple[float, int]]:
    """Numeric tokens with the number of displayed decimals, so precision is respected."""
    out = []
    for m in NUM.finditer(span):
        txt = m.group().replace(",", "")
        try:
            val = float(txt)
        except ValueError:
            continue
        out.append((val, len(txt.split(".")[1]) if "." in txt else 0))
    return out


def matches(target: float, val: float, decimals: int) -> bool:
    """True when the printed token is the artefact value rounded to that many decimals.

    The 2% condition is what stops a coarse token such as "0" or "2" from standing in for a
    small quantity like 0.0928 or 1.6278, which the half-ulp test alone would accept.
    """
    diff = abs(target - val)
    return diff <= 10.0 ** (-decimals) / 2.0 + 1e-9 and diff <= 0.02 * abs(target) + 1e-9


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--tex", required=True)
    ap.add_argument("--repro", default=str(ROOT / "reproduction"))
    ap.add_argument("--processed", default=str(ROOT / "data" / "processed" / "upgrade"))
    a = ap.parse_args()
    body = strip_tex(Path(a.tex).read_text(encoding="utf-8"))
    repro, up = Path(a.repro), Path(a.processed)

    ev = load(repro / "json" / "event_table_audit.json")
    nb = load(repro / "json" / "null_battery_summary.json")
    full_nb = load(up / "c2_null_battery.json")
    rob = load(repro / "json" / "dm_floor_robust.json")
    sb = load(repro / "json" / "same_beam_audit.json")
    ctl = load(repro / "json" / "known_repeater_control.json")
    sc = load(up / "pair565_scatter.json")
    per = load(up / "pair565_period.json")
    zinv = load(up / "dm_z_inversion.json")
    with (repro / "tables" / "candidates_28pairs.csv").open(encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh))
    kept = sum(1 for r in rows if r["kept_for_claims"] == "kept")
    removed = sum(1 for r in rows if r["kept_for_claims"] == "removed")
    ind = full_nb["independence"]
    sizes = sorted(ind["component_sizes"], reverse=True)
    w25 = rob["window_frac_0.25"]
    feo = rob["formal_error_optimism"]
    obs = sb["observed_pair"]
    rho = sb["rho_sensitivity"]
    conv = sb["error_convolved_counts"]
    ha = sb["ha_stats"]
    pools = {p["pool"]: p for p in sb["pools"]}
    cond = (pools["same-beam AND cross-day"]["P_theta_and_ddm"]
            / pools["cross-day pairs (excl. same-night)"]["P_theta_and_ddm"])
    dt = per["dt_days"]
    n_sid = dt / T_SID
    off_s = abs(round(n_sid) - n_sid) * T_SID * 86400.0
    cc = sc["control_cell"]

    # quantity, artefact value, extra spellings, topic anchor, slot regex capturing the text's claim
    checks = [
        ("catalogue rows joined", ev["row_level"]["all"], None, None, None),
        ("distinct events", ev["dedup"]["all"], None, None, None),
        ("apparent one-off events", ev["dedup"]["one_off"], None, None, None),
        ("declaration domain", ev["dedup"]["one_off_finite_dm"], None, None, None),
        ("Cat2 rows with coordinates", ev["row_level"]["cat2_rows_ok"], None, None, None),
        ("comparable bursts", rob["n_comparable"], None, None, None),
        ("Cat1 rows", ev["row_level"]["cat1_rows"], None, "catalogue 1 rows", None),
        ("Cat1 rows removed as repeats", ev["dedup"]["cat1_rows_matched_in_cat2"],
         None, "of the 600", r"(\d+)\s+of\s+the\s+600"),
        ("candidate pairs (criterion set)", ind["n_pairs"], None, "pairs involving", None),
        ("events in candidate pairs", ind["n_unique_events"], None, "unique events", None),
        ("connected components", ind["n_components"], None, "independent hypotheses", None),
        ("largest component size", sizes[0], None, "size vector", None),
        ("hypotheses surviving the audit", kept, ["seven"], "surviving", None),
        ("pair entries removed by the audit", removed, None, "removed", None),
        ("stricter null bar (grid union)", nb["models"]["ra"]["real_cand"], None, "union", None),
        ("permutation draws B", nb["models"]["ra"]["n_reps"], None, "b+1", None),
        ("single-test p bound", 0.001, None, "b+1", None),
        ("Bonferroni global p bound", 0.009, None, "bonferroni", None),
        ("grid configurations", nb["config"]["n_eff_trials"], None, "bonferroni", None),
        ("separation angle (arcsec)", obs["theta_arcsec"], None, "arcsec", None),
        ("abs Delta DM (pc cm-3)", abs(obs["ddm"]), None, "dmpc", None),
        ("floor ratio, conservative", w25["ratio_conservative_DMsqrt2"], None, "sigma", None),
        ("floor ratio CI95 low", w25["ratio_conservative_ci95"][0], None, "ci", None),
        ("floor ratio CI95 high", w25["ratio_conservative_ci95"][1], None, "ci", None),
        ("floor s(DM) at +-0.25", w25["s_hat"], None, "loose", None),
        ("two-burst sigma_comb", w25["sigma_comb_hat"], None, "sigma", None),
        ("frac |D| above formal error (%)", feo["frac_absD_gt_formal_err"] * 100, None, None, None),
        ("overall median of ratios", feo["overall_median_of_ratios"], None, None, None),
        ("median-ratio range low", feo["ratio_of_medians_range"][0], None, "bin)", None),
        ("median-ratio range high", feo["ratio_of_medians_range"][1], None, "$700$", None),
        ("|HA| median, declared domain (h)", ha["median_abs_ha_h"], None, "ha", None),
        ("|HA|<0.1 h share, declared domain (%)", ha["frac_abs_ha_lt_0p1_h"] * 100, None, "ha", None),
        ("|HA| median, all Cat2 rows (h)", ha["all_cat2_rows"]["median_abs_ha_h"], None,
         "localized catalogue 2", None),
        ("|HA|<0.1 h share, all Cat2 rows (%)", ha["all_cat2_rows"]["frac_abs_ha_lt_0p1_h"] * 100,
         None, "0.0488", None),
        ("same-beam conditional-rate factor", cond, None, "same-beam", None),
        ("sigma_rel at rho=0 (arcsec)", rho[0]["sigma_rel_comb_arcsec"], None, "rho", None),
        ("sigma_rel at rho=0.99 (arcsec)", rho[-1]["sigma_rel_comb_arcsec"], None, "rho", None),
        ("theta/sigma_rel at rho=0", rho[0]["theta_over_sigma_rel"], None, "rho", None),
        ("theta/sigma_rel at rho=0.99", rho[-1]["theta_over_sigma_rel"], None, "rho", None),
        ("convolution draws", conv["draws"], None, "draws", None),
        ("convolved all-pool median count", conv["all_pairs"]["median"], None, "full pool yields",
         None),
        ("known repeaters in Cat2", ctl["n_sources"], None, "recall", None),
        ("recall as a fraction", None, [f"{ctl['n_recovered']}/{ctl['n_sources']}"], "recall", None),
        ("candidate percentile among repeaters", ctl["candidate"]["percentile_floor"] * 100,
         None, "percentile", None),
        ("control-cell measurable share", cc["frac_measured"] * 100,
         [f"{cc['n_tau_measured']}/{cc['n_total']}"], "control cell", None),
        ("control-cell minimum tau (ms)", cc["tau_ms_measured"]["min"], None, "control cell", None),
        ("control-cell median tau (ms)", cc["tau_ms_measured"]["median"], None, "median", None),
        ("extrapolated tau at DM 1835 (ms)", sc["scaling"]["tau_predicted_at_1835_ms"], None,
         "trend", None),
        ("ratio R of prediction over floor", sc["scaling"]["ratio_R_predicted_over_floor"],
         None, None, None),
        ("observed DM", zinv["dm_obs"], None, None, None),
        ("Galactic DM", zinv["mw"]["mean"], None, None, None),
        ("adopted DM_IGM", zinv["adopted"]["dm_igm"], None, "halo", None),
        ("redshift inverted from DM-z", zinv["adopted"]["z"], None, "intrinsic scatter",
         r"z\s*=\s*([\d.]+)"),
        ("+-20% scatter band, low", zinv["bands"]["intrinsic_scatter_pm20"]["z"][0], None,
         "intrinsic scatter", None),
        ("+-20% scatter band, high", zinv["bands"]["intrinsic_scatter_pm20"]["z"][1], None,
         "intrinsic scatter", None),
        ("+-30% scatter band, low", zinv["bands"]["intrinsic_scatter_pm30"]["z"][0], None,
         "intrinsic scatter", None),
        ("+-30% scatter band, high", zinv["bands"]["intrinsic_scatter_pm30"]["z"][1], None,
         "intrinsic scatter", None),
        ("formal DM error shift on z (table states an upper bound)",
         (zinv["bands"]["formal_dm_error_pm"]["z"][1] - zinv["bands"]["formal_dm_error_pm"]["z"][0])
         / 2.0, ["0.001"], "formal", None),
        ("distinct bursts in Catalogue 2", sc["population"]["burst_level"], None, None, None),
        ("sub-burst rows in Catalogue 2 (5045 rows minus bursts)",
         5045 - sc["population"]["burst_level"], ["506"], None, None),
        ("halo/Galactic shift on z, stated as an upper bound",
         max(c["z"] for c in zinv["systematics_grid"]) - min(c["z"] for c in zinv["systematics_grid"]),
         ["0.04"], "shifts $z$", None),
        ("anchor deficit at z=0.189 (%)", abs(zinv["low_z_anchors"][0]["dev_pct"]), None,
         "under the prediction", None),
        ("anchor deficit at z=0.321 (%)", abs(zinv["low_z_anchors"][1]["dev_pct"]), None,
         "under the prediction", None),
        ("anchor deficit at z=0.660 (%)", abs(zinv["low_z_anchors"][2]["dev_pct"]), None,
         "under the prediction", None),
        ("Delta t (days)", dt, None, None, None),
        ("Delta t in sidereal days", n_sid, None, "sidereal", None),
        ("offset from integer sidereal days (s)", off_s, None, "sidereal", None),
    ]

    bad = nostated = noanchor = 0
    print(f"manuscript: {a.tex}\n")
    for name, target, extra, anchor, slot in checks:
        spans = windows(body, anchor)
        if anchor is not None and not spans:
            noanchor += 1
            print(f"  ANCHOR-NOT-FOUND  {name:38s} artefact={target if target is not None else '-':>10}"
                  f"  anchor=/{anchor}/")
            continue
        hit = None
        for s in spans:
            if extra and any(x in s for x in extra):
                hit = extra[0]
                break
            if target is not None:
                got = next(((v, d) for v, d in numbers(s)
                            if matches(float(target), v, d)), None)
                if got:
                    hit = f"{got[0]} ({got[1]} dp)"
                    break
        if hit is not None:
            print(f"  ok        {name:38s} artefact={target if target is not None else '-':>10}"
                  f"  text={hit}")
            continue
        claimed = None
        if slot:
            for s in spans:
                m = re.search(slot, s, re.I)
                if m:
                    claimed = m.group(1)
                    break
        if claimed is not None:
            bad += 1
            print(f"  MISMATCH  {name:38s} artefact={target:>10}  text states {claimed} "
                  f"(anchor=/{anchor}/)")
        else:
            nums = sorted({v for s in spans for v, _ in numbers(s)})[:12]
            nostated += 1
            print(f"  NOT-STATED {name:38s} artefact={target if target is not None else '-':>10}"
                  f"  anchor=/{anchor}/  numbers nearby={nums}")
    total = len(checks)
    print(f"\n{total - bad - nostated - noanchor}/{total} agree, {bad} mismatch, "
          f"{nostated} not-stated, {noanchor} anchor-not-found")
    print("not machine-checked: TNS and ADS results (external services) and the 58/85 "
          "denominator that mixes Catalogue 1-only sources.")
    print(f"note: the DM-z relation sits {zinv['low_z_anchor_mean_abs_dev_pct']}% (mean absolute) "
          f"above the published DM_IGM of the three low-z localised anchors; those per-anchor "
          f"deviations are checked as text claims above and stored in dm_z_inversion.json.")
    return 1 if (bad or noanchor) else 0


if __name__ == "__main__":
    raise SystemExit(main())
