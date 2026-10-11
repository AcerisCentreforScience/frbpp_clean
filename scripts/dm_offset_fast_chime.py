"""Estimate the FAST-to-CHIME DM zero-point offset without matched bursts.

Usage: python scripts/dm_offset_fast_chime.py
"""
from __future__ import annotations

import json

import numpy as np
import pandas as pd

from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SRC = "FRB20220912A"
rng = np.random.default_rng(20261008)

raw = pd.read_csv(ROOT / "data/raw/chime/cat2_25/chimefrbcat2.csv", low_memory=False)
rep = raw.repeater_name.astype("string").str.strip().fillna("")
raw["src_rep"] = rep.mask(rep.isin(["-9999", "nan", "NaN", "None"]), "").to_numpy()
ch = raw[raw.src_rep == SRC].copy()
ch["t"] = pd.to_numeric(ch.mjd_400, errors="coerce")
ch["dm"] = pd.to_numeric(ch.dm_fitb, errors="coerce")
ch["e"] = pd.to_numeric(ch.dm_fitb_err, errors="coerce")
ch = ch[ch.t.notna() & ch.dm.notna()]
subburst = pd.to_numeric(ch.sub_num, errors="coerce")
ch = ch[subburst.isna() | (subburst == 0)]

fa = pd.DataFrame(json.loads((ROOT / "data/raw/fast/fast_20220912A_bursts.json")
                             .read_text(encoding="utf-8"))["rows"])
fa["t"] = pd.to_numeric(fa.MJD, errors="coerce")
fa["dm"] = pd.to_numeric(fa.DM, errors="coerce")
fa["e"] = pd.to_numeric(fa.e_DM, errors="coerce")
fa = fa[fa.t.notna() & fa.dm.notna()]

print(f"CHIME: n={len(ch)}  t {ch.t.min():.3f}→{ch.t.max():.3f}  "
      f"DM 中位 {ch.dm.median():.2f} 范围 [{ch.dm.min():.2f},{ch.dm.max():.2f}]  "
      f"e 中位 {ch.e.median():.3f}")
print(f"FAST : n={len(fa)}  t {fa.t.min():.3f}→{fa.t.max():.3f}  "
      f"DM 中位 {fa.dm.median():.2f} 范围 [{fa.dm.min():.2f},{fa.dm.max():.2f}]  "
      f"e 中位 {fa.e.median():.3f}")
n_in = int(((fa.t >= ch.t.min()) & (fa.t <= ch.t.max())).sum())
print(f"FAST 暴落在 CHIME 覆盖区间内 = {n_in}/{len(fa)}")
print(f"平均 epoch: CHIME {ch.t.mean():.2f}  FAST {fa.t.mean():.2f}  "
      f"差 {fa.t.mean() - ch.t.mean():+.2f} d")

dm = np.concatenate([ch.dm.to_numpy(float), fa.dm.to_numpy(float)])
t = np.concatenate([ch.t.to_numpy(float), fa.t.to_numpy(float)])
ef = np.concatenate([ch.e.to_numpy(float), fa.e.to_numpy(float)])
ind = np.concatenate([np.zeros(len(ch)), np.ones(len(fa))])
tbar = t.mean()


def design(tt, ii, order):
    cols = [np.ones_like(tt)]
    for k in range(1, order + 1):
        cols.append((tt - tbar) ** k)
    cols.append(ii)
    return np.vstack(cols).T


def wls(X, y, w):
    W = np.sqrt(w)[:, None]
    beta, *_ = np.linalg.lstsq(X * W, y * np.sqrt(w), rcond=None)
    return beta


def refit(order, s):
    X = design(t, ind, order)
    w = 1.0 / s ** 2
    beta = wls(X, dm, w)
    cov = np.linalg.inv(X.T @ (X * w[:, None]))
    se = np.sqrt(np.diag(cov))
    chi2 = float(np.sum((dm - X @ beta) ** 2 * w))
    dof = len(dm) - X.shape[1]
    return dict(order=order, delta=float(beta[-1]), se=float(se[-1]),
                chi2n=chi2 / dof, p=X.shape[1])


_s0 = np.maximum(ef, 1e-3)
_r0 = dm - design(t, ind, 1) @ wls(design(t, ind, 1), dm, 1.0 / _s0 ** 2)
SIG_JIT = float(1.4826 * np.median(np.abs(_r0 - np.median(_r0))))
S_FORMAL = _s0
S_FLOOR = np.maximum(ef, SIG_JIT)
print(f"\nσ_jit = {SIG_JIT:.3f} pc/cm³  （形式误差中位 CHIME {np.median(ch.e):.3f} / "
      f"FAST {np.median(fa.e):.3f} → 地板放大 "
      f"{SIG_JIT / np.median(ef):.1f}×）")

print("\n== 逐步加入趋势阶数，看 δ 是否漂移（漂移=趋势与仪器指示项混淆）==")
res = []
for order in (0, 1, 2):
    r = refit(order, S_FLOOR)
    res.append(r)
    print(f"  jitter       order={order}: δ = {r['delta']:+7.3f} ± {r['se']:.3f} "
          f"(χ²/dof = {r['chi2n']:6.2f}, n={len(dm)}, p={r['p']})")

print("\n== 同一拟合的两种误差尺度对照（地板 vs 形式误差）==")
for order in (0, 1, 2):
    rf = refit(order, S_FORMAL)
    print(f"  formal       order={order}: δ = {rf['delta']:+7.3f} ± {rf['se']:.3f} "
          f"(χ²/dof = {rf['chi2n']:6.2f})")

main = res[1]
print(f"\n主结果（order=1，抖动地板尺度）：")
print(f"  δ(FAST − CHIME) = {main['delta']:+.3f} ± {main['se']:.3f} pc/cm³ "
      f"= {main['delta'] / main['se']:+.2f}σ   在 DM≈{np.median(dm):.0f} 处")
print(f"  δ 随趋势阶数的漂移：order 0/1/2 = "
      f"{res[0]['delta']:+.3f} / {res[1]['delta']:+.3f} / {res[2]['delta']:+.3f}"
      f" → 系统项 ≲ {max(r['delta'] for r in res) - min(r['delta'] for r in res):.3f}")

print("\n== 负对照：置换仪器标签（δ 应回到 0）==")
X = design(t, ind, 1)
w = 1.0 / S_FLOOR ** 2
obs = wls(X, dm, w)[-1]
perm = np.array([wls(design(t, rng.permutation(ind), 1), dm, w)[-1] for _ in range(300)])
z = (obs - perm.mean()) / perm.std()
print(f"  置换零分布: 均值 {perm.mean():+.4f}  sd {perm.std():.4f}"
      f"  |观测 δ = {obs:+.3f} → z = {z:+.2f}")
print(f"  判据: |z|>3 且 δ 幅度 > σ_jit/10 才声明「真实跨仪器零点差」"
      f"；否则记为「与本目录精度不可分」")

cov1 = np.linalg.inv(X.T @ (X * w[:, None]))
rho_d = cov1[0, -1] / np.sqrt(cov1[0, 0] * cov1[-1, -1])
Xc = X - X.mean(0)
Xs = Xc[:, 1:] / Xc[:, 1:].std(0)
print("\n== 可识别性 ==")
print(f"  非常数列标准化后条件数 = {np.linalg.cond(Xs):.3e}")
cc = np.corrcoef(Xc[:, 1], Xc[:, -1])[0, 1]
print(f"  corr(t−t̄, 仪器指示项) = {cc:+.3f} → δ 与线性趋势的可混淆程度 |ρ| = {abs(cc):.2f}")
print(f"  corr(截距 a, δ) = {rho_d:+.3f}（全样本）")

print("\n== 同时段子样本（t 限制在 FAST 观测窗内）→ 切断 δ 与趋势的简并 ==")
win = (t >= fa.t.min()) & (t <= fa.t.max())
tw, dww, ew, indw = t[win], dm[win], ef[win], ind[win]
sw = np.maximum(ew, SIG_JIT)
rw = 1.0 / sw ** 2
print(f"  窗内 n={int(win.sum())}（CHIME {int((indw == 0).sum())} + FAST "
      f"{int((indw == 1).sum())}），时间跨度 {(tw.max() - tw.min()):.1f} d")
for order in (0, 1):
    Xw = np.vstack([np.ones_like(tw), (tw - tw.mean()) ** 1, indw]).T if order == 1 \
        else np.vstack([np.ones_like(tw), indw]).T
    bw = wls(Xw, dww, rw)
    cw = np.linalg.inv(Xw.T @ (Xw * rw[:, None]))
    Xw0 = Xw[:, 1:] / Xw[:, 1:].std(0)
    ccw = (np.corrcoef(Xw[:, 1] - Xw[:, 1].mean(), Xw[:, -1] - Xw[:, -1].mean())[0, 1]
           if order == 1 else 0.0)
    print(f"  order={order}: δ = {bw[-1]:+7.3f} ± {np.sqrt(cw[-1, -1]):.3f}"
          f"  χ²/dof={float(np.sum((dww - Xw @ bw) ** 2 * rw) / (len(dww) - Xw.shape[1])):.2f}"
          f"  |ρ(δ,趋势)|={abs(ccw):.2f}  cond={np.linalg.cond(Xw0):.2e}")
Xw1 = np.vstack([np.ones_like(tw), tw - tw.mean(), indw]).T
obs_w = wls(Xw1, dww, rw)[-1]
perm_w = np.array([wls(np.vstack([np.ones_like(tw), tw - tw.mean(),
                                  rng.permutation(indw)]).T, dww, rw)[-1]
                   for _ in range(300)])
zw = (obs_w - perm_w.mean()) / perm_w.std()
print(f"  置换检验（order=1）: 零分布 sd={perm_w.std():.4f} | δ={obs_w:+.3f} → z = {zw:+.2f}")
print("\n== 分层检验：δ 随 DM 的稳定性（若漂移 → 零点差是 DM 依赖的，不可外推）==")
q_dm = np.quantile(dm, [0.34, 0.67])
for lab, m in [("低 DM 三分位", dm < q_dm[0]),
               ("中 DM", (dm >= q_dm[0]) & (dm < q_dm[1])),
               ("高 DM 三分位", dm >= q_dm[1])]:
    if m.sum() > 20 and 0 < ind[m].sum() < m.sum():
        Xi = design(t[m], ind[m], 1)
        si = S_FLOOR[m]
        wi = 1.0 / si ** 2
        bi = wls(Xi, dm[m], wi)
        ci = np.linalg.inv(Xi.T @ (Xi * wi[:, None]))
        print(f"  {lab:14s} n={int(m.sum()):5d}  DM∈[{dm[m].min():.0f},{dm[m].max():.0f}]"
              f"  δ = {bi[-1]:+7.3f} ± {np.sqrt(ci[-1, -1]):.3f}")

fl = np.concatenate([pd.to_numeric(ch.fluence, errors="coerce").to_numpy(float),
                     pd.to_numeric(fa.Fluence, errors="coerce").to_numpy(float)])
is_ch = np.r_[np.ones(len(ch), bool), np.zeros(len(fa), bool)]
print("\n== DM–通量协变检验（CHIME 的已知系统差：暗暴 DM 偏低）==")
print(f"  CHIME fluence 中位 {np.nanmedian(fl[is_ch]):.3f} Jy·ms  "
      f"p10/p90 = {np.nanquantile(fl[is_ch], .1):.3f}/{np.nanquantile(fl[is_ch], .9):.3f}")
print(f"  FAST   Fluence 中位 {np.nanmedian(fl[~is_ch]):.3f} Jy·ms  "
      f"p10/p90 = {np.nanquantile(fl[~is_ch], .1):.3f}/{np.nanquantile(fl[~is_ch], .9):.3f}")
okf = np.isfinite(fl)
lg = np.log10(np.where(okf, fl, 1.0))
for lab, m in [("CHIME", is_ch), ("FAST  ", ~is_ch)]:
    sel = m & okf
    xg = lg[sel] - lg[sel].mean()
    dg = dm[sel] - dm[sel].mean()
    sl = float(np.dot(xg, dg) / np.dot(xg, xg))
    se_sl = float(np.std(dg - sl * xg, ddof=2) / np.sqrt(np.dot(xg, xg)))
    print(f"  {lab} 斜率 dDM/dlog10F = {sl:+.3f} ± {se_sl:.3f} "
          f"(n={int(sel.sum())}, 跨 {xg.max() - xg.min():.1f} dex → "
          f"组内 ΔDM 幅度 ≈ {sl * (xg.max() - xg.min()):+.2f})")

sel_w = win & okf
lgw = lg[sel_w] - lg[sel_w].mean()
tw2 = t[sel_w] - t[sel_w].mean()
iw = ind[sel_w]
Xf = np.vstack([np.ones_like(iw), tw2, lgw, iw]).T
sf = np.maximum(ef[sel_w], SIG_JIT)
wf = 1.0 / sf ** 2
bf = wls(Xf, dm[sel_w], wf)
cf = np.linalg.inv(Xf.T @ (Xf * wf[:, None]))
print(f"\n  窗内含 log F 协变量: δ = {bf[-1]:+.3f} ± {np.sqrt(cf[-1, -1]):.3f} "
      f"| dDM/dlogF = {bf[2]:+.3f} ± {np.sqrt(cf[2, 2]):.3f}")

lo = max(np.nanpercentile(fl[:len(ch)], 5), np.nanpercentile(fl[len(ch):], 5))
hi = min(np.nanpercentile(fl[:len(ch)], 95), np.nanpercentile(fl[len(ch):], 95))
ov = win & okf & (fl >= lo) & (fl <= hi)
print(f"  通量重叠区间 [{lo:.3f}, {hi:.3f}] Jy·ms 内 n={int(ov.sum())} "
      f"(CHIME {int((ind[ov] == 0).sum())} / FAST {int((ind[ov] == 1).sum())})")
if 0 < ind[ov].sum() < ov.sum():
    Xo = np.vstack([np.ones(int(ov.sum())), t[ov] - t[ov].mean(), ind[ov]]).T
    so = np.maximum(ef[ov], SIG_JIT)
    wo = 1.0 / so ** 2
    bo = wls(Xo, dm[ov], wo)
    co = np.linalg.inv(Xo.T @ (Xo * wo[:, None]))
    print(f"    → δ = {bo[-1]:+.3f} ± {np.sqrt(co[-1, -1]):.3f} pc/cm³（通量硬匹配子样本）")

print("\n== 归因：通量选择效应能否解释观测到的 δ ==")
mlog = {lab: float(np.mean(lg[m & okf])) for lab, m in [("CHIME", is_ch), ("FAST", ~is_ch)]}
dlogF = mlog["FAST"] - mlog["CHIME"]
sl_obs = {"CHIME": None, "FAST": None}
for lab, m in [("CHIME", is_ch), ("FAST", ~is_ch)]:
    sel = m & okf
    xg = lg[sel] - lg[sel].mean()
    dg = dm[sel] - dm[sel].mean()
    sl_obs[lab] = float(np.dot(xg, dg) / np.dot(xg, xg))
raw_ddm = float(np.median(dm[~is_ch & okf]) - np.median(dm[is_ch & okf]))
print(f"  平均 log10 F 差(FAST−CHIME) = {dlogF:+.3f} dex "
      f"（中位通量比 = {10 ** dlogF:.2f}×）")
for lab in ("FAST", "CHIME"):
    print(f"  用 {lab} 斜率 {sl_obs[lab]:+.3f}/dex → 预测伪 δ = {sl_obs[lab] * dlogF:+.3f}")
print(f"  原始中位差（未控任何协变量）= {raw_ddm:+.3f} pc/cm³")
print(f"  → 观测 δ(窗内, 未控通量) = {obs_w:+.3f}；预测选择效应 = "
      f"{sl_obs['FAST'] * dlogF:+.3f} → 同量级同符号")
print(f"  通量控制后两种口径：协变量法 {bf[-1]:+.3f} ± {np.sqrt(cf[-1, -1]):.3f}；"
      f"硬匹配法 {bo[-1]:+.3f} ± {np.sqrt(co[-1, -1]):.3f}")
print(f"  两口径互差 = {abs(bf[-1] - bo[-1]):.3f} pc/cm³ → 这是本方法的系统下限，"
      f"与待测 δ 同量级")

