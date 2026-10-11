"""Check every figure against the journal resolution and font-size requirements.

Usage: python scripts/audit_figure_specs.py
"""
from __future__ import annotations

import io
import re
import struct
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
TEX = ROOT / "docs" / "submission" / "raa_manuscript.tex"
# raa.cls sets textwidth to 145mm for the final layout (200mm under [referee]); the
# generators design at 15.00 cm, so print scales figures by 14.5/15.0.
TEXT_W_CM = 14.5
CM_PER_IN = 2.54

DPI_COLOUR = 300.0
DPI_GREY = 500.0
DPI_LINE = 1000.0
MIN_PT = 7.0

VECTOR_EXT = {".pdf", ".eps", ".ai", ".ps", ".svg"}


def is_colour(p: Path, max_px: int = 400_000) -> bool:
    try:
        from PIL import Image
        im = Image.open(p)
        im.thumbnail((1200, 1200))
        rgb = im.convert("RGB")
    except Exception:
        return True
    a = rgb.tobytes()
    n = len(a) // 3
    if n == 0:
        return False
    diff = 0
    step = max(1, n // 60_000)
    seen = 0
    for i in range(0, n, step):
        r, g, b = a[3 * i], a[3 * i + 1], a[3 * i + 2]
        seen += 1
        if max(r, g, b) - min(r, g, b) > 6:
            diff += 1
    return (diff / max(seen, 1)) > 0.002


def png_size_dpi(p: Path) -> tuple[tuple[int, int], tuple[float, float] | None]:
    b = p.read_bytes()
    w, h = struct.unpack(">II", b[16:24])
    dpi = None
    i = 8
    while i + 8 <= len(b):
        ln = struct.unpack(">I", b[i:i + 4])[0]
        typ = b[i + 4:i + 8]
        if typ == b"pHYs" and ln >= 9:
            px, py, unit = struct.unpack(">IIB", b[i + 8:i + 17])
            if unit == 1:
                dpi = (px * 0.0254, py * 0.0254)
            break
        if typ == b"IDAT":
            break
        i += 12 + ln
    return (w, h), dpi


def tex_figures() -> list[tuple[int, str, bool]]:
    out = []
    for n, line in enumerate(TEX.read_text(encoding="utf-8").splitlines(), 1):
        m = re.search(r"\\includegraphics(\[[^\]]*\])?\{([^}]+)\}", line)
        if m:
            out.append((n, m.group(2), "textwidth" in (m.group(1) or "")))
    return out


def design_specs():
    specs: dict[str, tuple[float, float, float]] = {}
    mg = (ROOT / "scripts" / "make_figs.py").read_text(encoding="utf-8")
    for m in re.finditer(r"\ndef (\w+)\(", mg):
        blk = mg[m.start(): (mg.find("\ndef ", m.start() + 1) % len(mg)) or len(mg)]
        f = re.search(r"figsize=\(([0-9.]+),\s*[0-9.]+\)", blk)
        fs = [float(x) for x in re.findall(r"fontsize=([0-9.]+)", blk)]
        base = re.search(r'"font\.size":\s*([0-9.]+)', blk)
        base = float(base.group(1)) if base else 7.0
        if f and fs:
            specs[m.group(1)] = (float(f.group(1)), min(fs), base)
    mf = (ROOT / "scripts" / "make_fast_advantage.py").read_text(encoding="utf-8")
    f = re.search(r"figsize=\(([0-9.]+),\s*[0-9.]+\)", mf)
    fs = [float(x) for x in re.findall(r"fontsize=([0-9.]+)", mf)]
    base = re.search(r'"font\.size"\]\s*=\s*([0-9.]+)', mf)
    if f and fs:
        specs["make_fast_advantage"] = (float(f.group(1)), min(fs),
                                        float(base.group(1)) if base else 7.0)
    return specs


FIG_FN = {
    "fig2_event_dedup_mechanism": "fig2",
    "fig3_joint_isolation": "fig3",
    "fig4_sidereal_quantization": "fig4",
    "fig5_null_battery": "fig5_null_battery",
    "fig6_scattering_discriminant": "fig6",
    "fig8_fast_advantage": "make_fast_advantage",
}


RETEXT_SPEC = {"flow": "SPEC_FLOW", "geometry": "SPEC_GEOM", "dual": "SPEC_DUAL",
               "scatter": "SPEC_SCATTER", "roadmap": "SPEC_ROADMAP",
               "composite": "SPEC_COMPOSITE"}


def retext_min_px(key: str) -> int:
    txt = (ROOT / "scripts" / "build_retext.py").read_text(encoding="utf-8")
    var = RETEXT_SPEC.get(key)
    if not var:
        return 0
    s = txt.find(f"{var} = [")
    if s < 0:
        return 0
    e = txt.find("\n]", s)
    blk = txt[s: e if e > 0 else len(txt)]
    sizes = [int(x) for x in re.findall(r'S\("lat",\s*(\d+)', blk)]
    return min(sizes) if sizes else 0



def vector_min_label_pt(pdf_path: Path) -> float | None:
    """Smallest line-level font size of a vector figure, read from the PDF.

    Each text line contributes its largest run: a mathtext subscript or a log-axis exponent is drawn
    at 0.7x its base size by design and shares the body font family, so a per-span minimum would
    flag correct typography. Lines made only of numerals or signs are tick labels and are skipped.
    """
    try:
        import pymupdf
    except Exception:
        return None
    sizes = []
    numerals = set("0123456789.−-+eEx")
    with pymupdf.open(str(pdf_path)) as doc:
        for page in doc:
            for blk in page.get_text("dict")["blocks"]:
                for line in blk.get("lines", []):
                    spans = [(s["text"].strip(), float(s["size"])) for s in line.get("spans", [])]
                    spans = [(t, z) for t, z in spans if t]
                    if not spans:
                        continue
                    if all(all(c in numerals for c in t) for t, _ in spans):
                        continue
                    sizes.append(max(z for _, z in spans))
    return min(sizes) if sizes else None
def main() -> int:
    specs = design_specs()
    hard, soft = [], []
    print(f"版心宽 {TEXT_W_CM} cm = {TEXT_W_CM/CM_PER_IN:.3f} in"
          f"   |   RAA 分档：彩色 {DPI_COLOUR:.0f} / 黑白 {DPI_GREY:.0f} / 线条 {DPI_LINE:.0f} dpi"
          f"   |   软线 {MIN_PT:.0f} pt（字号，内部要求）")
    print("=" * 108)
    print(f"{'行':>4}  {'文件':44s} {'类型':5s} {'像素/设计宽':>14s} "
          f"{'有效dpi':>9s} {'最小字号pt':>10s}  判定")
    print("-" * 108)

    for line_no, rel, full in tex_figures():
        p = (TEX.parent / rel).resolve()
        stem = Path(rel).stem
        name = Path(rel).name
        ext = Path(rel).suffix.lower()
        width_cm = TEXT_W_CM if full else None
        if not p.exists():
            print(f"{line_no:>4}  {name:44s} {'-':5s} {'缺失':>14s} {'-':>9s} {'-':>10s}  "
                  f"!! 文件不存在")
            hard.append((line_no, name, "缺失"))
            continue

        if ext in VECTOR_EXT:
            min_pt = vector_min_label_pt(p)
            if min_pt is not None and min_pt < MIN_PT:
                soft.append((line_no, name, f"{min_pt:.1f} pt"))
                print(f"{line_no:>4}  {name:44s} {'矢量':5s} {'-':>14s} {'∞':>9s} "
                      f"{min_pt:>10.1f}  OK 矢量（RAA 位图分档不适用）  ~ 字号 {min_pt:.1f} pt < {MIN_PT:.0f}")
            else:
                print(f"{line_no:>4}  {name:44s} {'矢量':5s} {'-':>14s} {'∞':>9s} "
                      f"{'-':>10s}  OK 矢量，RAA 位图分档不适用")
            continue

        size, dpi = png_size_dpi(p)
        w_px = size[0]
        eff_dpi = w_px / (width_cm / CM_PER_IN) if width_cm else None

        min_pt = None
        if stem in FIG_FN:
            spec = specs.get(FIG_FN[stem])
            if spec:
                fw_in, fs_min, base = spec
                scale = (width_cm / CM_PER_IN) / fw_in
                min_pt = min(fs_min, base) * scale
        elif retext_min_px(stem):
            min_pt = retext_min_px(stem) / w_px * width_cm / CM_PER_IN * 72.0

        colour = is_colour(p)
        tier = DPI_COLOUR if colour else DPI_GREY
        kind = "位图·彩色" if colour else "位图·灰度"

        notes = []
        if eff_dpi is not None and eff_dpi < tier:
            notes.append(f"!! {eff_dpi:.0f} dpi < {tier:.0f}（{'彩色' if colour else '黑白'}档）")
            hard.append((line_no, name, f"{eff_dpi:.0f} dpi"))
        if min_pt is not None and min_pt < MIN_PT:
            notes.append(f"~ 字号 {min_pt:.1f} pt < {MIN_PT:.0f}")
            soft.append((line_no, name, f"{min_pt:.1f} pt"))
        tag = "  ".join(notes) if notes else "OK"
        print(f"{line_no:>4}  {name:44s} {kind:5s} {str(size):>14s} "
              f"{(f'{eff_dpi:.0f}' if eff_dpi else '-'):>9s} "
              f"{(f'{min_pt:.1f}' if min_pt else '-'):>10s}  {tag}")

    print("=" * 108)
    print(f"硬违规（RAA 分档：彩色 {DPI_COLOUR:.0f} / 黑白 {DPI_GREY:.0f} / 线条 {DPI_LINE:.0f} dpi）: "
          f"{len(hard)} 项" +
          ("" if not hard else "  ->  " + "; ".join(f"L{n} {f}" for n, f, _ in hard)))
    print(f"软提醒（{MIN_PT:.0f} pt 字号，内部要求非 RAA 条文）: {len(soft)} 项" +
          ("" if not soft else f"  ->  最小 {min(float(v.split()[0]) for _, _, v in soft):.1f} pt"))
    return 1 if hard else 0


if __name__ == "__main__":
    sys.exit(main())
