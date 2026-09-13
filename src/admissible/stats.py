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
    """Two-sided binomial p-value, by exact summation in log space.

    Used instead of a fixed allele-balance window because a fixed window
    mis-scales with depth: 3/10 is unremarkable for a true heterozygote,
    60/200 is not, and both sit at an allele balance of 0.30.

    The arithmetic is done on log probabilities rather than on
    ``math.comb(n, k) * p**k * ...`` because the binomial coefficient overflows
    a float long before real data runs out of depth: at n = 7166 - an ordinary
    pileup depth in panel or amplicon data - the direct form raises
    OverflowError and takes the whole check down with it.  Above
    ``_EXACT_MAX`` trials the exact sum is replaced by a normal approximation
    with a continuity correction, which keeps the cost constant instead of
    linear in depth.  That approximation tracks the exact value to within a
    fraction of a percent near any threshold anyone sets here (0.03% at
    p = 5e-3, 0.4% at p = 2e-8 for n = 20000) and drifts to order-of-magnitude
    accuracy far out in the tail - which changes no classification, because a
    p-value of 1e-45 and one of 1e-44 are the same verdict.
    """
    if n <= 0 or not (0.0 < p < 1.0):
        return float("nan")
    k = min(max(k, 0), n)
    if n > _EXACT_MAX:
        mean = n * p
        sd = math.sqrt(n * p * (1.0 - p))
        if sd == 0.0:
            return 1.0
        z = (abs(k - mean) - 0.5) / sd
        if z <= 0.0:
            return 1.0
        return min(1.0, math.erfc(z / math.sqrt(2.0)))

    log_obs = _binom_logpmf(k, n, p)
    tol = log_obs + 1e-7
    total = 0.0
    for i in range(n + 1):
        lp = _binom_logpmf(i, n, p)
        if lp <= tol:
            total += math.exp(lp)
    return min(1.0, total)


# Above this many trials the exact sum is both unnecessary and slow: the normal
# approximation is accurate to well beyond the precision any threshold here uses.
_EXACT_MAX = 20_000


def _binom_logpmf(k: int, n: int, p: float) -> float:
    return (
        math.lgamma(n + 1)
        - math.lgamma(k + 1)
        - math.lgamma(n - k + 1)
        + k * math.log(p)
        + (n - k) * math.log1p(-p)
    )


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
