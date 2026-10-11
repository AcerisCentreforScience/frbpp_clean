"""Data contracts: column sets, dtype checks, routing invariants, coordinate identity assertions."""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

FLAG_RFI = 1 << 0
FLAG_BASEBAND = 1 << 1
FLAG_BEAM_FAR = 1 << 2
FLAG_EXCLUDED = 1 << 3

CATEGORIES = ("repeater", "suspected", "one_off")

IDENTITY_TOL_DEG = 10.0 / 3600.0


@dataclass
class BurstCatalog:

    t: np.ndarray
    source: np.ndarray
    dm: np.ndarray
    dm_err: np.ndarray
    snr: np.ndarray
    fluence: np.ndarray
    band_lo: np.ndarray
    band_hi: np.ndarray
    ra: np.ndarray
    dec: np.ndarray
    flags: np.ndarray
    meta: dict = field(default_factory=dict)

    _ARRAY_FIELDS = ("source", "dm", "dm_err", "snr", "fluence",
                     "band_lo", "band_hi", "ra", "dec", "flags")

    def __post_init__(self) -> None:
        n = len(self.t)
        for name in self._ARRAY_FIELDS:
            if len(getattr(self, name)) != n:
                raise ValueError(
                    f"BurstCatalog.{name} 长度 {len(getattr(self, name))} != t 长度 {n}")
        if n and not np.all(np.isfinite(self.t)):
            raise ValueError("BurstCatalog.t 含非有限值（时间系统必须显式可用）")

    @property
    def n_bursts(self) -> int:
        return int(len(self.t))


@dataclass
class SourceRecord:

    name: str
    ra: float
    dec: float
    category: str
    n_bursts: int
    n_episodes: int
    min_episodes: int = 5
    cross_match_status: str = "unmatched"
    multi_component: bool = False
    meta: dict = field(default_factory=dict)

    def __post_init__(self) -> None:
        if self.category not in CATEGORIES:
            raise ValueError(f"未知 category: {self.category!r}（可选 {CATEGORIES}）")


@dataclass
class ObservationLog:

    t_start: np.ndarray
    t_end: np.ndarray
    sens: np.ndarray
    tgt_ra: float
    tgt_dec: float
    beam_ra: np.ndarray
    beam_dec: np.ndarray
    meta: dict = field(default_factory=dict)

    def __post_init__(self) -> None:
        m = len(self.t_start)
        for name in ("t_end", "sens", "beam_ra", "beam_dec"):
            if len(getattr(self, name)) != m:
                raise ValueError(f"ObservationLog.{name} 长度 != t_start 长度")
        if m and np.any(self.t_end <= self.t_start):
            raise ValueError("ObservationLog 存在非正长度观测段")


@dataclass
class EfficiencyCurve:

    snr_grid: np.ndarray
    dm_grid: np.ndarray
    eff: np.ndarray
    meta: dict = field(default_factory=dict)

    def __post_init__(self) -> None:
        if self.eff.shape != (len(self.snr_grid), len(self.dm_grid)):
            raise ValueError(
                f"EfficiencyCurve.eff 形状 {self.eff.shape} 与网格不匹配")


_PIPELINE_BY_CATEGORY: dict[str, tuple[str, ...]] = {
    "repeater": ("periodicity", "sr", "dm"),
    "suspected": ("census",),
    "one_off": ("upgrade_side_channel", "population_baseline"),
}


def pipeline_for(category: str) -> tuple[str, ...]:
    try:
        return _PIPELINE_BY_CATEGORY[category]
    except KeyError:
        raise ValueError(f"未知 category: {category!r}") from None


def enters_periodicity_likelihood(category: str) -> bool:
    return category == "repeater"


def angular_separation_deg(ra1, dec1, ra2, dec2):
    ra1, dec1, ra2, dec2 = (np.radians(np.asarray(x, dtype=float))
                            for x in (ra1, dec1, ra2, dec2))
    d_ra = ra2 - ra1
    d_dec = dec2 - dec1
    a = (np.sin(d_dec / 2.0) ** 2
         + np.cos(dec1) * np.cos(dec2) * np.sin(d_ra / 2.0) ** 2)
    return np.degrees(2.0 * np.arcsin(np.sqrt(np.clip(a, 0.0, 1.0))))


def match_coordinates(ra1, dec1, ra2, dec2,
                      tol_deg: float = IDENTITY_TOL_DEG) -> bool:
    return bool(angular_separation_deg(ra1, dec1, ra2, dec2) <= tol_deg)