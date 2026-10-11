"""Generate the re-typeset reference figures with their label specifications.

Usage: python scripts/build_retext.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
from PIL import Image, ImageDraw, ImageFont

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.stdout.reconfigure(encoding="utf-8", errors="replace")
from retext_figures import (  # noqa: E402
    REFS, OUT, REF_DIR, font, erase, put, erase_pill, sample_bg, flush_rotated,
)

WM = {
    "composite": [1130, 640, 1376, 762],
    "flow": [1770, 1018, 2030, 1130],
    "geometry": [1770, 1018, 2030, 1130],
    "scatter": [1770, 1018, 2035, 1130],
    "roadmap": [1770, 1018, 2030, 1130],
    "dual": [1770, 1018, 2030, 1130],
}

NAVY = (31, 58, 95)
BLUE = (43, 91, 152)
LBLUE = (91, 141, 191)
SLATE = (95, 112, 130)
GREY = (120, 120, 120)
DARK = (40, 40, 40)
RED = (200, 40, 46)
WHITE = (255, 255, 255)


def S(key="lat", size=22, color=DARK, align="center", valign="center",
      bold=False, cjk=False, rotate=0):
    return dict(key=key, size=size, color=color, align=align,
                valign=valign, bold=bold, cjk=cjk, rotate=rotate)


def _is_cjk(text: str) -> bool:
    return any("\u4e00" <= c <= "\u9fff"
               or "\u3000" <= c <= "\u303f"
               or "\uff00" <= c <= "\uffef" for c in text)


def _fnt(st, line, size):
    cjk = st["cjk"] or _is_cjk(line)
    return font(size, cjk=cjk, bold=st["bold"])


def _ink_map(fpath: str, ch: str, size: int = 24):
    from PIL import Image as _I
    f = ImageFont.truetype(fpath, size)
    im = _I.new("L", (46, 36), 255)
    ImageDraw.Draw(im).text((4, 4), ch, font=f, fill=0)
    return np.asarray(im)


def missing_glyphs(text: str, fpath: str) -> list[str]:
    out = []
    for ch in set(text):
        if ord(ch) < 128:
            continue
        try:
            if np.array_equal(_ink_map(fpath, ch), _ink_map(fpath, chr(0xE0A0))):
                out.append(ch)
        except Exception:
            out.append(ch)
    return sorted(out)


def fit_text(d, text: str, st, box) -> tuple[str, int]:
    bw = box[2] - box[0]
    floor = max(11, int(round(st["size"] * 0.72)))
    size = st["size"]
    lines = text.split("\n")

    def widest(sz):
        return max(d.textlength(ln, font=_fnt(st, ln, sz)) for ln in lines)

    while size > floor and widest(size) > bw - 2:
        size -= 1
    if widest(size) > bw - 2 and len(lines) == 1:
        idx = [i for i, c in enumerate(text) if c == " "]
        if idx:
            best = min(idx, key=lambda i: abs(i - len(text) // 2))
            cand = [text[:best], text[best + 1:]]
            w2 = max(d.textlength(ln, font=_fnt(st, ln, st["size"])) for ln in cand)
            if w2 <= bw - 2:
                return "\n".join(cand), st["size"]
    return "\n".join(lines), size


def apply(im: Image.Image, spec, lang: int):
    """spec: list of (box, zh, en, style) / (box, zh, en, style, erasebox)"""
    for item in spec:
        box = item[0]
        text = item[1] if lang == 0 else item[2]
        eb = item[4] if len(item) > 4 else box
        bg = item[5] if len(item) > 5 else None
        erase(im, eb, pad=2, bg=bg)
    d = ImageDraw.Draw(im)
    fitted = []
    for i, item in enumerate(spec):
        box = item[0]
        text = item[1] if lang == 0 else item[2]
        st = item[3]
        if not text:
            continue
        probe = _fnt(st, text, st["size"])
        miss = missing_glyphs(text, getattr(probe, "path", ""))
        if miss:
            raise SystemExit(
                f"[{lang and 'en' or 'zh'} spec 第 {i} 条] 字体缺字形 "
                + " ".join(f"U+{ord(c):04X} {c!r}" for c in miss)
                + f"  文本={text[:60]!r}  字体={getattr(probe, 'path', '?')}"
                + "\n（PIL 会把缺字形画成空心方框 = 图里的\"豆腐块\"。请换字或改写法。）")
        new_text, size = fit_text(d, text, st, box)
        if new_text != text or size != st["size"]:
            fitted.append((i, text[:34], f"折行={new_text != text} 字号 {st['size']}->{size}"))
        st2 = dict(st)
        st2["size"] = size
        put(d, box, new_text, _fnt(st2, new_text, size), color=st["color"],
            align=st["align"], valign=st["valign"],
            rotate=st.get("rotate", 0))
    flush_rotated(im)
    im.info["retext_fitted"] = fitted
    return im


def load(key: str) -> Image.Image:
    p = REF_DIR / REFS[key]
    im = Image.open(p).convert("RGB")
    erase_pill(im, WM[key], pad=3)
    return im


SPEC_SCATTER = [
    ([1706, 30, 2018, 74], "p ≤ 0.009（4 种零模型，Bonferroni）",
     "p ≤ 0.009 (4 null models, Bonferroni)", S("lat", 26, DARK, "right", "center")),
    ([1050, 343, 1492, 378], "FRB20200511C × FRB20211128E",
     "FRB20200511C × FRB20211128E", S("lat", 25, RED)),
    ([1050, 381, 1345, 417], "(θ = 23.9 arcsec, ΔDM = 0.18 pc cm⁻³)",
     "(θ = 23.9 arcsec, ΔDM = 0.18 pc cm⁻³)", S("lat", 22, RED)),
    ([1281, 692, 1564, 726], "极端孤立区", "extreme isolated region",
     S("lat", 22, RED, "left", "center")),
    ([243, 961, 615, 995], "3665 个表观单次事件（背景）",
     "3665 apparent single events (background)", S("lat", 22, DARK, "left", "center")),
]

SPEC_FLOW = [
    ([60, 196, 495, 240], "CHIME/FRB 目录 1 & 2", "CHIME/FRB Catalog 1 & 2",
     S("lat", 30, NAVY, "center", "center", bold=True)),
    ([179, 243, 375, 280], "（2018 – 2023）", "(2018 – 2023)",
     S("lat", 24, SLATE)),
    ([60, 296, 400, 334], "5608 目录行 → 5014 独立事件",
     "5608 rows → 5014 events", S("lat", 22, WHITE)),
    ([216, 335, 339, 389], "5014", "5014", S("lat", 40, NAVY, bold=True)),
    ([140, 395, 412, 430], "独立事件", "independent events", S("lat", 24, WHITE)),
    ([217, 591, 338, 644], "3665", "3665", S("lat", 40, NAVY, bold=True)),
    ([148, 645, 407, 687], "表观单次事件子空间", "apparent single-event subspace",
     S("lat", 24, DARK)),
    ([221, 687, 334, 718], "（有限 DM + 坐标）", "(finite DM + position)", S("lat", 20, SLATE)),
    ([722, 257, 1411, 298], "跨名称同源关联框架", "Cross-Name Association Framework",
     S("lat", 32, NAVY, bold=True)),
    ([675, 350, 780, 398], "步骤 1", "Step 1", S("lat", 22, WHITE, bold=True)),
    ([1017, 354, 1122, 397], "步骤 2", "Step 2", S("lat", 22, WHITE, bold=True)),
    ([1357, 354, 1460, 397], "步骤 3", "Step 3", S("lat", 22, WHITE, bold=True)),
    ([615, 417, 838, 453], "时空锥 +", "Space-Time Cone +", S("lat", 21, DARK)),
    ([611, 448, 842, 486], "DM 窗口配对判据", "DM Window Pairing Criterion", S("lat", 21, DARK)),
    ([676, 483, 779, 515], "判据", "", S("lat", 20, DARK)),
    ([1007, 422, 1133, 454], "并查集", "Union-Find", S("lat", 21, DARK)),
    ([952, 450, 1184, 491], "假设归并算法", "Hypothesis Merging", S("lat", 21, DARK)),
    ([1010, 484, 1130, 522], "算法", "Algorithm", S("lat", 20, DARK)),
    ([1345, 427, 1472, 459], "零模型", "Null Model", S("lat", 21, DARK)),
    ([1314, 460, 1507, 491], "置换检验", "Permutation Test", S("lat", 21, DARK)),
    ([1408, 546, 1524, 584], "p ≤ 0.009", "p ≤ 0.009", S("lat", 21, RED, bold=True)),
    ([1343, 727, 1464, 755], "检验统计量", "Test Statistic", S("lat", 20, DARK)),
    ([700, 754, 740, 785], "DM", "DM", S("lat", 20, DARK)),
    ([764, 656, 850, 688], "Δt、ΔDM", "Δt, ΔDM", S("lat", 20, DARK)),
    ([1676, 383, 1982, 408], "最强候选", "STRONGEST CANDIDATE", S("lat", 20, (92,112,135), bold=True)),
    ([1709, 462, 1967, 502], "FRB20200511C", "FRB20200511C", S("lat", 24, WHITE, bold=True)),
    ([1710, 564, 1968, 602], "FRB20211128E", "FRB20211128E", S("lat", 24, WHITE, bold=True)),
    ([150, 505, 420, 545], "剔除 594 个重复行", "594 duplicate rows removed",
     S("lat", 19, SLATE)),
]


SPEC_DUAL = [
    ([368, 31, 814, 80], "散射时标 τ 与 DM", "Scattering Time vs. DM",
     S("lat", 30, DARK, "center", "center", bold=True)),
    ([1315, 13, 1820, 60], "DM 预算分解", "DM Budget Decomposition",
     S("lat", 30, DARK, "center", "center", bold=True)),
    ([242, 146, 514, 180], "已知 FRB / 类 FRB 样本", "Known FRB / FRB-like data",
     S("lat", 22, DARK, "left", "center")),
    ([242, 190, 663, 224], "湍流屏 τ–DM 趋势（银河系标度）",
     "Turbulent screen τ–DM trend (Galactic)", S("lat", 22, DARK, "left", "center")),
    ([240, 229, 476, 269], "±1σ（预期范围）", "±1σ (expected range)",
     S("lat", 22, DARK, "left", "center")),
    ([226, 910, 600, 944], "主候选（未测得散射）", "Candidate (no scattering detection)",
     S("lat", 22, DARK, "left", "center")),
    ([653, 646, 894, 682], "排除为", "Excluded as", S("lat", 23, RED, "center", "center", bold=True)),
    ([653, 675, 999, 710], "湍流屏主导", "turbulence-screen dominant",
     S("lat", 23, RED, "center", "center", bold=True)),
    ([1236, 240, 1534, 282], "总 DM", "Total DM", S("lat", 23, DARK, "center", "center", bold=True),
     [1228, 238, 1542, 285]),
    ([1236, 278, 1534, 318], "1835 pc cm⁻³", "1835 pc cm⁻³", S("lat", 23, DARK, "center", "center", bold=True),
     [1222, 276, 1548, 322]),
    ([1337, 484, 1442, 524], "IGM DM", "IGM DM", S("lat", 23, WHITE, bold=True)),
    ([1353, 525, 1428, 565], "1671", "1671", S("lat", 24, WHITE, bold=True)),
    ([1348, 558, 1431, 597], "(91.1%)", "(91.1%)", S("lat", 21, WHITE, bold=True)),
    ([1341, 764, 1441, 799], "银晕 DM", "Halo DM", S("lat", 22, WHITE, bold=True)),
    ([1361, 796, 1418, 834], "40", "40", S("lat", 23, WHITE, bold=True)),
    ([1349, 826, 1431, 865], "(2.2%)", "(2.2%)", S("lat", 21, WHITE, bold=True)),
    ([1330, 886, 1451, 918], "银河系 DM", "Galactic DM", S("lat", 22, WHITE, bold=True)),
    ([1361, 916, 1420, 954], "124", "124", S("lat", 23, WHITE, bold=True)),
    ([1349, 948, 1430, 984], "(6.7%)", "(6.7%)", S("lat", 21, WHITE, bold=True)),
    ([1700, 824, 2020, 872], "IGM DM（1671 pc cm-3）", "IGM DM (1671 pc cm⁻³)",
     S("lat", 21, DARK, "left", "center"),
     [1698, 822, 2036, 874]),
    ([1700, 887, 2020, 931], "银晕 DM（40 pc cm-3）", "Halo DM (40 pc cm⁻³)",
     S("lat", 21, DARK, "left", "center"),
     [1698, 885, 2036, 933]),
    ([1700, 947, 2020, 995], "银河系 DM（124 pc cm-3）", "Galactic DM (124 pc cm⁻³)",
     S("lat", 21, DARK, "left", "center"),
     [1698, 945, 2036, 997]),
    ([1424, 129, 1599, 165], "H1：z ≈ 1.5", "H1: z ≈ 1.5", S("lat", 21, DARK, "left", "center", bold=True)),
    ([1424, 155, 1599, 195], "高红移重复暴候选", "high-redshift repeater",
     S("lat", 20, DARK, "left", "center")),
    ([1801, 142, 1981, 179], "H2：低红移", "H2: low-redshift", S("lat", 21, DARK, "left", "center", bold=True)),
    ([1808, 175, 1974, 206], "致密本地介质", "dense local medium",
     S("lat", 20, DARK, "left", "center")),
]


SPEC_GEOM = [
    ([633, 52, 1184, 100], "CHIME 固定天顶望远镜", "CHIME Fixed Zenith Telescope",
     S("lat", 32, NAVY, "center", "center", bold=True)),
    ([668, 102, 1150, 141], "源赤纬 δ 下的观测几何",
     "Geometry with Source Declination δ", S("lat", 25, NAVY)),
    ([33, 124, 296, 168], "时间关系", "Time Relation",
     S("lat", 30, NAVY, "left", "center", bold=True)),
    ([69, 238, 409, 295], "565 太阳日", "565 solar days",
     S("lat", 32, NAVY, "center", "center", bold=True)),
    ([69, 356, 409, 409], "= 567 恒星日", "= 567 sidereal days",
     S("lat", 32, NAVY, "center", "center", bold=True)),
    ([158, 422, 320, 464], "（恰为整数）", "(integerized)", S("lat", 23, SLATE)),
    ([166, 626, 311, 668], "几何锁定", "Geometric Lock",
     S("lat", 27, NAVY, "center", "center", bold=True)),
    ([72, 846, 410, 882], "地球自转与恒星日周期", "Earth rotation and sidereal",
     S("lat", 22, SLATE)),
    ([53, 885, 427, 924], "由几何锁定（565 = 567）", "cycle are locked by geometry",
     S("lat", 22, SLATE)),
    ([959, 211, 1160, 242], "地球自转轴", "Earth's rotation axis",
     S("lat", 21, SLATE, "left", "center")),
    ([968, 301, 1048, 336], "天顶", "Zenith", S("lat", 21, SLATE, "center", "center")),
    ([1178, 276, 1256, 308], "源", "Source", S("lat", 21, SLATE, "center", "center")),
    ([1177, 303, 1316, 334], "（赤纬 δ）", "(declination δ)", S("lat", 21, SLATE, "left", "center")),
    ([1057, 501, 1230, 564], "CHIME 固定\n天顶望远镜", "CHIME fixed\nzenith telescope",
     S("lat", 21, SLATE, "left", "center")),
    ([529, 749, 705, 783], "地球自转", "Earth's rotation", S("lat", 21, SLATE, "left", "center")),
    ([528, 781, 704, 821], "（1 恒星日）", "(1 sidereal day)", S("lat", 21, SLATE, "left", "center")),
    ([1364, 100, 1691, 149], "关键判别项", "Key Discriminants",
     S("lat", 28, NAVY, "left", "center", bold=True)),
    ([1369, 209, 1451, 241], "特征", "Feature", S("lat", 21, NAVY, "center", "center", bold=True)),
    ([1614, 192, 1799, 257], "周期假设\n（物理）", "Physical\nPeriod Hypothesis",
     S("lat", 20, NAVY, "center", "center", bold=True)),
    ([1819, 192, 2019, 257], "几何产物\n假设", "Geometric\nProduct Hypothesis",
     S("lat", 20, NAVY, "center", "center", bold=True)),
    ([1369, 306, 1567, 370], "长期相位稳定性", "Phase stability over\nlong timescales",
     S("lat", 20, DARK, "left", "center")),
    ([1616, 326, 1784, 357], "漂移 / 不稳定", "Drifts / not stable", S("lat", 20, DARK)),
    ([1886, 324, 1951, 357], "稳定", "Stable", S("lat", 20, DARK, "left", "center")),
    ([1367, 415, 1492, 480], "多源普适性", "Multi-source\nuniversality",
     S("lat", 20, DARK, "left", "center")),
    ([1615, 434, 1796, 468], "依赖源", "Source-dependent", S("lat", 20, DARK)),
    ([1886, 434, 1979, 466], "普适", "Universal", S("lat", 20, DARK, "left", "center")),
    ([1369, 532, 1492, 594], "周期值\n（565 vs 567）", "Period value\n(565 vs 567)",
     S("lat", 20, DARK, "left", "center")),
    ([1613, 533, 1743, 570], "非天然", "Not naturally", S("lat", 20, DARK)),
    ([1614, 565, 1727, 600], "整数化", "integerized", S("lat", 20, DARK)),
    ([1884, 529, 1997, 567], "由几何", "Integerized", S("lat", 20, DARK)),
    ([1884, 561, 2009, 598], "整数化", "by geometry", S("lat", 20, DARK)),
    ([1370, 640, 1526, 671], "与地球自转", "Correlation with", S("lat", 20, DARK, "left", "center")),
    ([1370, 671, 1518, 702], "的关联", "Earth's rotation", S("lat", 20, DARK, "left", "center")),
    ([1615, 641, 1747, 676], "无必然", "No necessary", S("lat", 20, DARK)),
    ([1613, 669, 1695, 704], "关系", "relation", S("lat", 20, DARK)),
    ([1887, 657, 2021, 689], "直接关联", "Directly linked", S("lat", 20, DARK, "left", "center")),
    ([1370, 746, 1593, 806], "对未观测源的\n预测力", "Predictive power for\nunobserved sources",
     S("lat", 20, DARK, "left", "center")),
    ([1614, 761, 1692, 793], "有限", "Limited", S("lat", 20, DARK)),
    ([1883, 758, 1937, 797], "高", "High", S("lat", 20, DARK, "left", "center")),
    ([1369, 852, 1561, 916], "对仪器变更的\n稳健性", "Robustness to\ninstrument changes",
     S("lat", 20, DARK, "left", "center")),
    ([1614, 865, 1705, 900], "敏感", "Sensitive", S("lat", 20, DARK)),
    ([1886, 867, 1959, 899], "稳健", "Robust", S("lat", 20, DARK, "left", "center")),
    ([1368, 958, 1586, 1019], "解释恒星日/\n太阳日关系", "Explains sidereal/solar\nrelation",
     S("lat", 20, DARK, "left", "center")),
    ([1611, 969, 1651, 1008], "否", "No", S("lat", 20, DARK)),
    ([1885, 972, 1925, 1006], "是", "Yes", S("lat", 20, DARK, "left", "center")),
]


SPEC_ROADMAP = [
    ([241, 485, 479, 528], "FAST 望远镜", "FAST Telescope",
     S("lat", 25, WHITE, bold=True)),
    ([174, 525, 551, 565], "第三次爆发确认", "Third Burst Confirmation",
     S("lat", 25, WHITE, bold=True)),
    ([800, 487, 1250, 526], "偏振 RM 测量", "Polarization RM Measurement",
     S("lat", 25, WHITE, bold=True)),
    ([781, 526, 1271, 566], "约束磁场环境", "Constrain Magnetic Environment",
     S("lat", 25, WHITE, bold=True)),
    ([1508, 465, 1882, 504], "VLBI / 光学对应体", "VLBI / Optical Counterpart",
     S("lat", 25, WHITE, bold=True)),
    ([1548, 503, 1841, 540], "高精度定位", "Precise Localization",
     S("lat", 25, WHITE, bold=True)),
    ([1559, 539, 1830, 582], "区分 H1 / H2", "Distinguish H1 / H2",
     S("lat", 25, WHITE, bold=True)),
    ([287, 574, 335, 622], "1", "1", S("lat", 30, WHITE, bold=True)),
    ([976, 574, 1024, 622], "2", "2", S("lat", 30, WHITE, bold=True)),
    ([1655, 574, 1703, 622], "3", "3", S("lat", 30, WHITE, bold=True)),
    ([233, 857, 478, 891], "关键科学产出", "Key Scientific Output",
     S("lat", 20, SLATE)),
    ([903, 857, 1147, 891], "关键科学产出", "Key Scientific Output", S("lat", 20, SLATE)),
    ([1577, 864, 1829, 898], "关键科学产出", "Key Scientific Output", S("lat", 20, SLATE)),
    ([245, 895, 464, 943], "可重复性确认", "Repeatability", S("lat", 28, BLUE, "center", "center", bold=True)),
    ([245, 940, 464, 981], "确认", "Confirmation", S("lat", 28, BLUE, "center", "center", bold=True)),
    ([857, 898, 1192, 939], "介质磁化诊断", "Magnetized Medium", S("lat", 28, BLUE, "center", "center", bold=True)),
    ([941, 939, 1108, 986], "诊断", "Diagnosis", S("lat", 28, BLUE, "center", "center", bold=True)),
    ([1575, 901, 1828, 941], "红移锚点", "Redshift Anchor", S("lat", 28, BLUE, "center", "center", bold=True)),
    ([1231, 159, 1282, 197], "RM", "RM", S("lat", 22, DARK)),
    ([1935, 223, 1960, 251], "H1", "H1", S("lat", 20, WHITE, "left", "center", bold=True)),
    ([1934, 277, 1961, 306], "H2", "H2", S("lat", 20, WHITE, "left", "center", bold=True)),
    ([1768, 747, 1808, 783], "H1", "H1", S("lat", 20, DARK, "left", "center")),
    ([1767, 791, 1812, 830], "H2", "H2", S("lat", 20, DARK, "left", "center")),
    ([1063, 763, 1115, 805], "RM", "RM", S("lat", 22, DARK)),
    ([769, 1051, 1276, 1089], "当前红移不确定窗：z ∈ [0.3, 2.5]（对应 H1；H2 为 z < 0.1）",
     "Current uncertainty window: z in [0.3, 2.5] (H1); H2 lies at z < 0.1",
     S("lat", 21, SLATE)),
]


SPEC_COMPOSITE = [
    ([15, 10, 240, 46], "1. 数据与框架", "1. Data and Framework",
     S("lat", 21, WHITE, "left", "center", bold=True),
     [12, 8, 320, 48]),
    ([414, 13, 720, 45], "2. 最强同源候选证据", "2. Strongest Candidate Evidence",
     S("lat", 21, WHITE, "left", "center", bold=True),
     [410, 11, 762, 47]),
    ([1006, 15, 1322, 45], "3. 物理机制与红移", "3. Physical Mechanism & Redshift",
     S("lat", 21, WHITE, "left", "center", bold=True),
     [1002, 13, 1368, 47]),
    ([18, 508, 400, 540], "4. 方法学守卫与验证", "4. Methodology Safeguard & Verification",
     S("lat", 21, WHITE, "left", "center", bold=True),
     [14, 506, 600, 542]),
    ([975, 546, 1256, 574], "FAST 后续与未来验证", "FAST Follow-up & Future Verification",
     S("lat", 20, NAVY, "left", "center", bold=True),
     [970, 542, 1260, 578]),
    ([20, 78, 112, 124], "CHIME/FRB\n目录 1 & 2", "CHIME/FRB\nCatalog 1 & 2",
     S("lat", 20, NAVY, "left", "center", bold=True)),
    ([259, 78, 314, 108], "5014", "5014", S("lat", 25, NAVY, bold=True)),
    ([258, 101, 365, 125], "独立事件", "independent events", S("lat", 18, SLATE)),
    ([258, 148, 313, 178], "3665", "3665", S("lat", 25, NAVY, bold=True)),
    ([257, 171, 379, 197], "表观单次事件", "apparent single events", S("lat", 18, SLATE)),
    ([20, 221, 280, 246], "跨名称同源关联框架", "Cross-name association framework",
     S("lat", 20, NAVY, "left", "center", bold=True)),
    ([33, 258, 162, 283], "时空锥", "Spatiotemporal cone", S("lat", 18, DARK, "left", "center", bold=True)),
    ([251, 258, 323, 283], "并查集", "Union-Find", S("lat", 18, DARK, "left", "center", bold=True)),
    ([222, 272, 350, 302], "假设归并", "hypothesis merging", S("lat", 17, DARK, "left", "center")),
    ([34, 276, 118, 298], "+ 时空窗判据", "+ space-time window", S("lat", 16, DARK, "left", "center")),
    ([25, 407, 115, 432], "· Δθ（位置）", "· Δθ (positional)", S("lat", 17, DARK, "left", "center")),
    ([25, 427, 108, 452], "· Δt（时间）", "· Δt (temporal)", S("lat", 17, DARK, "left", "center")),
    ([25, 445, 177, 471], "· ΔDM（色散量）", "· ΔDM (dispersion measure)", S("lat", 17, DARK, "left", "center")),
    ([533, 58, 844, 83], "FRB20200511C × FRB20211128E", "FRB20200511C × FRB20211128E",
     S("lat", 21, RED, "center", "center", bold=True)),
    ([436, 98, 614, 123], "位置–DM 联合空间", "Position–DM Joint Space",
     S("lat", 18, NAVY, "center", "center", bold=True)),
    ([714, 98, 938, 123], "散射时标 τ 对比", "Scattering Time (τ) Comparison",
     S("lat", 18, NAVY, "center", "center", bold=True)),
    ([714, 121, 938, 141], "(1.0–1.5 GHz)", "(1.0–1.5 GHz)", S("lat", 15, SLATE, "center", "center")),
    ([826, 122, 946, 148], "观测", "Observed", S("lat", 14, DARK, "left", "center")),
    ([826, 144, 946, 170], "DM 趋势预期", "DM trend prediction", S("lat", 14, DARK, "left", "center")),
    ([550, 186, 612, 215], "θ = 23.9″", "θ = 23.9″", S("lat", 18, RED, "left", "center", bold=True)),
    ([549, 202, 677, 232], "ΔDM = 0.18 pc cm⁻³", "ΔDM = 0.18 pc cm⁻³",
     S("lat", 18, RED, "left", "center", bold=True)),
    ([752, 269, 827, 291], "未测得", "No scattering", S("lat", 17, DARK, "left", "center")),
    ([764, 284, 815, 306], "散射", "detected", S("lat", 17, DARK, "left", "center")),
    ([750, 364, 844, 386], "FRB20200511C ×", "FRB20200511C ×", S("lat", 15, DARK)),
    ([755, 380, 840, 402], "FRB20211128E", "FRB20211128E", S("lat", 15, DARK)),
    ([864, 362, 945, 389], "典型 FRB", "Typical FRB", S("lat", 15, DARK)),
    ([864, 379, 945, 404], "（DM 相近）", "(at similar DM)", S("lat", 15, DARK)),
    ([483, 379, 651, 403], "角距 θ（角秒）", "Angular Separation (arcsec)", S("lat", 17, DARK)),
    ([457, 412, 641, 436], "● FRB20200511C × FRB20211128E", "● FRB20200511C × FRB20211128E",
     S("lat", 15, RED, "left", "center")),
    ([461, 434, 536, 456], "○ 其它事件", "○ Other events", S("lat", 15, DARK, "left", "center")),
    ([718, 433, 926, 458], "排除湍流屏模型", "Excludes turbulence screen model",
     S("lat", 17, NAVY, "center", "center", bold=True)),
    ([714, 451, 929, 475], "（低于 DM 趋势 > 2 个数量级）",
     "(>2 orders of magnitude below DM trend)", S("lat", 15, NAVY, "center", "center")),
    ([1010, 63, 1206, 88], "DM 预算分解", "DM Budget Decomposition",
     S("lat", 18, NAVY, "center", "center", bold=True)),
    ([1029, 112, 1059, 138], "124", "124", S("lat", 16, WHITE, bold=True)),
    ([1024, 126, 1064, 155], "(6.7%)", "(6.7%)", S("lat", 14, WHITE)),
    ([1093, 113, 1122, 136], "40", "40", S("lat", 16, WHITE, bold=True)),
    ([1089, 126, 1128, 155], "(2.2%)", "(2.2%)", S("lat", 14, WHITE)),
    ([1206, 110, 1242, 137], "1671", "1671", S("lat", 16, WHITE, bold=True)),
    ([1206, 128, 1244, 155], "(91.1%)", "(91.1%)", S("lat", 14, WHITE)),
    ([1316, 105, 1351, 128], "总 DM", "Total", S("lat", 15, DARK)),
    ([1316, 119, 1354, 145], "1835", "1835", S("lat", 18, NAVY, bold=True)),
    ([1311, 132, 1362, 163], "pc cm⁻³", "pc cm⁻³", S("lat", 14, DARK)),
    ([1031, 180, 1108, 202], "银河系 DM", "Galactic DM", S("lat", 15, DARK, "left", "center")),
    ([1137, 178, 1211, 203], "银晕 DM", "Halo DM", S("lat", 15, DARK, "left", "center")),
    ([1256, 177, 1310, 202], "IGM DM", "IGM DM", S("lat", 15, DARK, "left", "center")),
    ([1053, 329, 1125, 360], "z ≈ 1.5", "z ≈ 1.5", S("lat", 17, NAVY, "left", "center", bold=True)),
    ([1049, 360, 1132, 385], "高红移重复暴", "High-redshift", S("lat", 15, NAVY, "left", "center")),
    ([1002, 375, 1171, 404], "候选", "repeating burst candidate", S("lat", 15, DARK)),
    ([1219, 332, 1324, 360], "低红移", "Low-redshift", S("lat", 17, NAVY, "left", "center", bold=True)),
    ([1206, 359, 1336, 384], "致密本地介质", "Dense local medium", S("lat", 15, DARK, "left", "center")),
    ([1127, 417, 1233, 442], "(待验证)", "(TO BE VERIFIED)", S("lat", 15, SLATE)),
    ([107, 556, 326, 584], "565 太阳日 = 567 恒星日", "565 days = 567 sidereal days",
     S("lat", 20, NAVY, "left", "center")),
    ([89, 670, 345, 699], "CHIME 过境观测几何产物，", "CHIME transit observation geometry product,",
     S("lat", 17, SLATE, "left", "center")),
    ([89, 691, 212, 715], "非物理周期。", "not a physical period.", S("lat", 17, SLATE, "left", "center")),
    ([453, 557, 741, 582], "零模型检验（4 种置换模型）", "Zero-model Test (4 permutation models)",
     S("lat", 20, NAVY, "left", "center")),
    ([479, 598, 539, 623], "时间置换", "Time shift", S("lat", 16, DARK)),
    ([595, 598, 661, 623], "DM 置换", "DM shuffle", S("lat", 16, DARK)),
    ([702, 597, 794, 622], "位置置换", "Position shuffle", S("lat", 16, DARK)),
    ([833, 596, 903, 624], "标签置换", "Label swap", S("lat", 16, DARK)),
    ([477, 625, 539, 652], "p ≤ 0.009", "p ≤ 0.009", S("lat", 16, NAVY, bold=True)),
    ([596, 625, 659, 652], "p ≤ 0.009", "p ≤ 0.009", S("lat", 16, NAVY, bold=True)),
    ([717, 624, 779, 652], "p ≤ 0.009", "p ≤ 0.009", S("lat", 16, NAVY, bold=True)),
    ([836, 625, 900, 652], "p ≤ 0.009", "p ≤ 0.009", S("lat", 16, NAVY, bold=True)),
    ([557, 680, 812, 705], "显著超出随机（p ≤ 0.009，Bonferroni 校正）",
     "Significant beyond random (p ≤ 0.009, Bonferroni)", S("lat", 17, WHITE, bold=True)),
    ([1221, 584, 1285, 609], "第三次", "Third burst", S("lat", 16, DARK, "left", "center")),
    ([1218, 598, 1295, 625], "爆发确认", "confirmation", S("lat", 16, DARK, "left", "center")),
    ([1219, 640, 1323, 664], "RM 测量", "RM measurement", S("lat", 16, DARK, "left", "center")),
    ([1218, 686, 1331, 716], "高精度定位", "Precise localization", S("lat", 16, DARK, "left", "center")),
    ([1027, 706, 1073, 732], "FAST", "FAST", S("lat", 18, NAVY, bold=True)),
]


JOBS = [
    ("composite", SPEC_COMPOSITE),
    ("flow", SPEC_FLOW),
    ("geometry", SPEC_GEOM),
    ("scatter", SPEC_SCATTER),
    ("roadmap", SPEC_ROADMAP),
    ("dual", SPEC_DUAL),
]


def main():
    (OUT / "zh").mkdir(parents=True, exist_ok=True)
    (OUT / "en").mkdir(parents=True, exist_ok=True)
    for key, spec in JOBS:
        for lang, tag in ((0, "zh"), (1, "en")):
            im = load(key)
            apply(im, spec, lang)
            dst = OUT / tag / f"{key}.png"
            im.save(dst)
            fitted = im.info.get("retext_fitted") or []
            print(f"saved -> {dst}   自动适配 {len(fitted)} 条"
                  + ("" if fitted else "（无）"))
            for i, txt, how in fitted:
                print(f"      [{i:2d}] {txt!r:38s} {how}")


if __name__ == "__main__":
    main()
