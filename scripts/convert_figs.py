"""Convert manuscript figures into submission formats.

Usage: python scripts/convert_figs.py
Output: docs/submission/figures_submission/<name>.pdf and .tif
"""
from __future__ import annotations

import shutil
from pathlib import Path

from PIL import Image

# RAA final layout textwidth is 145 mm (raa.cls); [referee] is 200 mm.
TEXT_W_CM = 14.5

SRC: list[tuple[str, str, str]] = [
    ("hires/fig2_event_dedup_mechanism.png", "hires/fig2_event_dedup_mechanism.png",
     "fig2_event_dedup_mechanism"),
    ("hires/fig3_joint_isolation.png", "hires/fig3_joint_isolation.png",
     "fig3_joint_isolation"),
    ("hires/fig4_sidereal_quantization.png", "hires/fig4_sidereal_quantization.png",
     "fig4_sidereal_quantization"),
    ("hires/fig5_null_battery.png", "hires/fig5_null_battery.png",
     "fig5_null_battery"),
    ("hires/fig6_scattering_discriminant.png", "hires/fig6_scattering_discriminant.png",
     "fig6_scattering_discriminant"),
    ("hires/fig8_fast_advantage.png", "hires/fig8_fast_advantage.png",
     "fig8_fast_advantage"),
    ("vector/flow.pdf", "vector/flow_900dpi.png", "fig10_flow"),
    ("vector/geometry.pdf", "vector/geometry_900dpi.png", "fig11_geometry"),
    ("vector/dual.pdf", "vector/dual_900dpi.png", "fig7_dual"),
    ("vector/scatter.pdf", "vector/scatter_900dpi.png", "fig12_scatter"),
    ("vector/roadmap.pdf", "vector/roadmap_900dpi.png", "fig13_roadmap"),
]

FIGDIR = Path(__file__).resolve().parents[1] / "docs/figures"
OUT = Path(__file__).resolve().parents[1] / "docs/submission/figures_submission"


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    ok = vector = 0
    for pdf_rel, tif_rel, stem in SRC:
        pdf_src = FIGDIR / pdf_rel
        tif_src = FIGDIR / tif_rel
        if not pdf_src.exists() or not tif_src.exists():
            print(f"MISSING {pdf_src if not pdf_src.exists() else tif_src}")
            continue

        pdf_dst = OUT / f"{stem}.pdf"
        if pdf_src.suffix.lower() == ".pdf":
            shutil.copy2(pdf_src, pdf_dst)
            vector += 1
            kind = f"vector PDF {pdf_dst.stat().st_size / 1024:.0f} KB"
        else:
            im = Image.open(pdf_src)
            dpi = im.info.get("dpi", (900, 900))
            im.convert("RGB").save(pdf_dst, "PDF", resolution=float(dpi[0]))
            kind = f"raster PDF @{float(dpi[0]):.0f} dpi"

        im2 = Image.open(tif_src)
        w, h = im2.size
        eff = w / (TEXT_W_CM / 2.54)
        im2.convert("RGB").save(OUT / f"{stem}.tif", "TIFF",
                                compression="tiff_lzw", dpi=(900, 900))
        print(f"{stem:32s} {kind:26s} | TIFF {w}x{h} (全幅等效 {eff:.0f} dpi)")
        ok += 1

    print(f"\n{ok}/{len(SRC)} figures converted into {OUT}  "
          f"（其中真矢量 PDF {vector} 张）")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
