"""Experiment: rescale figure type sizes to the printed size without changing content."""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.figure as mfigure
import matplotlib.text as mtext

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

TRIAL = ROOT / "docs" / "figures" / "trial"
TRIAL.mkdir(parents=True, exist_ok=True)

TEXT_W_CM = 16.2
CM_PER_IN = 2.54
TARGET_PT = 7.0
MARGIN = 1.08

FIGURES = {
    "fig2": ("make_figs", "fig2", 7.4),
    "fig3": ("make_figs", "fig3", 7.6),
    "fig4": ("make_figs", "fig4", 7.4),
    "fig5": ("make_figs", "fig5_null_battery", 7.6),
    "fig6": ("make_figs", "fig6", 4.5),
    "fig8": ("make_fast_advantage", "main", 7.4),
}
STEM = {"fig2": "fig2_event_dedup_mechanism", "fig3": "fig3_joint_isolation",
        "fig4": "fig4_sidereal_quantization", "fig5": "fig5_null_battery",
        "fig6": "fig6_scattering_discriminant", "fig8": "fig8_fast_advantage"}


def capture_and_rescale(mod, fname: str, out_stem: str, factor: float,
                        target_w_in: float, layout: str = "asis",
                        out_suffix: str = "") -> dict:
    grabbed: dict = {}
    orig_savefig = mfigure.Figure.savefig

    def patched(self, fname_, *a, **k):
        grabbed["fig"] = self
        return None

    mfigure.Figure.savefig = patched
    try:
        getattr(mod, fname)()
    finally:
        mfigure.Figure.savefig = orig_savefig

    fig = grabbed.get("fig")
    if fig is None:
        return {"error": "未截获 Figure"}

    fonts_before, fonts_after = [], []
    for t in fig.findobj(mtext.Text):
        s = t.get_fontsize()
        fonts_before.append(s)
        t.set_fontsize(s * factor)
        fonts_after.append(s * factor)

    w0, h0 = fig.get_size_inches()
    if layout == "2x2":
        axes = [a for a in fig.axes]
        fig.set_size_inches(target_w_in, target_w_in * 0.85)
        left, right, bottom, top = 0.085, 0.985, 0.105, 0.900
        gap_w, gap_h = 0.30, 0.30
        gw = (right - left - gap_w) / 2
        gh = (top - bottom - gap_h) / 2
        for i, ax in enumerate(axes[:4]):
            r, c = divmod(i, 2)
            ax.set_position([left + c * (gw + gap_w),
                             top - gh - r * (gh + gap_h), gw, gh])
        for t in fig.texts:
            if t.get_position()[1] < 0.06:
                t.set_position((t.get_position()[0], 0.008))
    else:
        fig.set_size_inches(target_w_in, h0 * target_w_in / w0)

    pdf = TRIAL / f"{out_stem}{out_suffix}_trial.pdf"
    png = TRIAL / f"{out_stem}{out_suffix}_trial_900dpi.png"
    fig.savefig(pdf, format="pdf")
    fig.savefig(png, dpi=900)
    matplotlib.pyplot.close(fig)

    import re
    head = pdf.read_bytes()[:4000]
    m = re.search(rb"/MediaBox\s*\[([^\]]+)\]", head)
    page_w_cm = None
    eff_min_pt = None
    if m:
        v = [float(x) for x in m.group(1).split()]
        page_w_cm = (v[2] - v[0]) / 72 * CM_PER_IN
        eff_min_pt = (min(fonts_after) if fonts_after else 0) * TEXT_W_CM / page_w_cm

    return {"n_text": len(fonts_before),
            "min_before": min(fonts_before) if fonts_before else None,
            "min_after": round(min(fonts_after), 2) if fonts_after else None,
            "size_in": tuple(round(float(x), 3) for x in fig.get_size_inches()),
            "pdf_page_cm": round(page_w_cm, 2) if page_w_cm else None,
            "eff_min_pt_at_16.2cm": round(eff_min_pt, 2) if eff_min_pt else None,
            "pdf_kb": round(pdf.stat().st_size / 1024)}


def main() -> int:
    keys = sys.argv[1:] or ["fig5"]
    layout = "2x2" if "--2x2" in sys.argv else "asis"
    keys = [k for k in keys if not k.startswith("--")]
    for k in keys:
        modname, fname, w0 = FIGURES[k]
        mod = __import__(modname)
        MIN_DESIGN = {"fig2": 5.3, "fig3": 5.4, "fig4": 5.0,
                      "fig5": 4.9, "fig6": 5.4, "fig8": 5.7}
        factor = TARGET_PT * MARGIN / MIN_DESIGN[k]
        print(f"== {k} ({modname}.{fname})  最小设计字号 {MIN_DESIGN[k]} pt "
              f"-> 乘 {factor:.3f}，目标最小 {TARGET_PT*MARGIN:.1f} pt  布局={layout}",
              flush=True)
        for it in range(4):
            info = capture_and_rescale(mod, fname, STEM[k], factor,
                                       TEXT_W_CM / CM_PER_IN, layout=layout,
                                       out_suffix="" if layout == "asis" else "_2x2")
            eff = info.get("eff_min_pt_at_16.2cm")
            print(f"   iter{it}: factor={factor:.3f}  page={info['pdf_page_cm']} cm  "
                  f"eff_min={eff} pt", flush=True)
            if eff and eff >= TARGET_PT * MARGIN:
                break
            if eff:
                factor *= (TARGET_PT * MARGIN) / eff
        print("   ", info, flush=True)
    print(f"\n产出目录：{TRIAL}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
