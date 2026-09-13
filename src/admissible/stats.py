"""Small statistics helpers.  Stdlib only, all closed-form."""

from __future__ import annotations

import math


def median(xs: list[float]) -> float:
    if not xs:
        return float("nan")
    s = sorted(xs)
    n = len(s)
    mid = n // 2
    return s[mid] if n % 2 else (s[mid - 1] + s[mid]) / 2.0


def mad(xs: list[float]) -> float:
    """Median absolute deviation, scaled to be comparable to a standard deviation."""
    if len(xs) < 2:
        return float("nan")
    m = median(xs)
    return 1.4826 * median([abs(x - m) for x in xs])


def wilson_interval(successes: int, trials: int, z: float = 1.96) -> tuple[float, float]:
    """Wilson score interval: behaves at small n, where the normal interval does not."""
    if trials <= 0:
        return (float("nan"), float("nan"))
    p = successes / trials
    denom = 1.0 + z * z / trials
    centre = (p + z * z / (2 * trials)) / denom
    half = (z / denom) * math.sqrt(p * (1 - p) / trials + z * z / (4 * trials * trials))
    return (max(0.0, centre - half), min(1.0, centre + half))


def binom_two_sided_p(k: int, n: int, p: float = 0.5) -> float:
    """Exact two-sided binomial p-value.

    Used instead of a fixed allele-balance window because a fixed window
    mis-scales with depth: 3/10 is unremarkable for a true heterozygote,
    60/200 is not, and both sit at an allele balance of 0.30.
    """
    if n <= 0 or not (0.0 < p < 1.0):
        return float("nan")
    obs = _binom_pmf(k, n, p)
    tol = obs * (1 + 1e-7)
    return min(1.0, sum(pmf for i in range(n + 1) if (pmf := _binom_pmf(i, n, p)) <= tol))


def _binom_pmf(k: int, n: int, p: float) -> float:
    return math.comb(n, k) * (p**k) * ((1 - p) ** (n - k))


def largest_gap(values: list[float], lo: float, hi: float) -> tuple[float, float] | None:
    """Widest empty interval within [lo, hi] in a sorted 1-D sample.

    Used to calibrate a boundary on the cohort itself rather than hard-coding
    one, which is the only defensible way to threshold a depth-dependent metric.
    """
    pts = sorted(v for v in values if not math.isnan(v))
    if len(pts) < 2:
        return None
    best: tuple[float, float] | None = None
    best_w = 0.0
    for a, b in zip(pts, pts[1:], strict=False):
        left, right = max(a, lo), min(b, hi)
        if right - left > best_w:
            best_w, best = right - left, (a, b)
    return best if best_w > 0 else None
