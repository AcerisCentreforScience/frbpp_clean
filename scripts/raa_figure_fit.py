"""Check that replacement labels fit their erase boxes and that the glyphs exist."""
import sys
from pathlib import Path

import numpy as np
from PIL import Image, ImageDraw, ImageFont

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
sys.stdout.reconfigure(encoding="utf-8", errors="replace")
import build_retext as B  # noqa: E402
import retext_figures as R  # noqa: E402

FONT_PX = 24


def glyph_map(path: str, ch: str) -> np.ndarray:
    f = ImageFont.truetype(path, FONT_PX)
    im = Image.new("L", (46, 36), 255)
    ImageDraw.Draw(im).text((4, 4), ch, font=f, fill=0)
    return np.asarray(im)


def covered(path: str, ch: str) -> bool:
    try:
        return not np.array_equal(glyph_map(path, ch), glyph_map(path, chr(0xE0A0)))
    except Exception:
        return False


COV: dict[tuple[str, str], bool] = {}


def has(path: str, ch: str) -> bool:
    key = (path, ch)
    if key not in COV:
        COV[key] = covered(path, ch)
    return COV[key]


SPECS = {
    "flow": B.SPEC_FLOW, "dual": B.SPEC_DUAL, "geometry": B.SPEC_GEOM,
    "roadmap": B.SPEC_ROADMAP, "scatter": B.SPEC_SCATTER, "composite": B.SPEC_COMPOSITE,
}

n_over = n_missing = n_shrunk = 0
for key, spec in SPECS.items():
    print(f"\n===== {key}  ({len(spec)} 条) =====")
    for i, item in enumerate(spec):
        box, zh, en, st = item[0], item[1], item[2], item[3]
        eb = item[4] if len(item) > 4 else box
        fn = R._find(R.LAT_BOLD if st.get("bold") else R.LAT_REG)
        if st.get("key") == "cjk":
            fn = R._find(R.CJK_BOLD if st.get("bold") else R.CJK_REG)
        f = ImageFont.truetype(fn, st["size"])
        d = ImageDraw.Draw(Image.new("L", (8, 8)))
        text2, size2 = B.fit_text(d, en, st, eb)
        lines = text2.split("\n")
        wmax = max(d.textlength(ln, font=ImageFont.truetype(fn, size2)) for ln in lines)
        bw = eb[2] - eb[0]
        miss = sorted({c for c in en if ord(c) > 127 and not has(fn, c)})
        flags = []
        if wmax > bw + 2:
            flags.append(f"适配后仍溢出 {wmax:.0f}px > 框宽 {bw}px")
            n_over += 1
        elif size2 < st["size"]:
            flags.append(f"已适配：字号 {st['size']}->{size2}"
                         + ("（触到下限，文案应改短）" if size2 <= max(11, round(st["size"] * 0.72)) else ""))
            n_shrunk += 1
        if miss:
            flags.append("缺字形 " + " ".join(f"U+{ord(c):04X} {c!r}" for c in miss))
            n_missing += 1
        if flags:
            print(f"  [{i:2d}] {en[:46]!r}")
            print(f"       框 {box} 擦 {eb} 字号 {st['size']} -> " + "；".join(flags))
print(f"\n合计：适配后仍溢出 {n_over} 条，缺字形 {n_missing} 条，靠缩字号/折行才放下 {n_shrunk} 条")
sys.exit(1 if (n_over or n_missing) else 0)
