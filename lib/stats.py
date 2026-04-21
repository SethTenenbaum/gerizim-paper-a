"""
stats.py — Statistical test wrappers.

Thin wrappers around scipy.stats that make the test logic explicit
and testable. Each function documents what it computes and why.
"""

from dataclasses import dataclass
from typing import List, Tuple

from scipy.stats import binomtest, chisquare, fisher_exact, spearmanr, norm
from scipy import stats as scipy_stats


# ---------------------------------------------------------------------------
# Binomial enrichment test
# ---------------------------------------------------------------------------

@dataclass
class BinomialResult:
    """Result of a one-tailed binomial test."""
    observed: int
    total: int
    null_rate: float
    expected: float
    ratio: float
    p_value: float

    @property
    def observed_rate(self) -> float:
        return self.observed / self.total if self.total else 0.0

    def summary(self) -> str:
        return (
            f"{self.observed}/{self.total} = {self.observed_rate:.1%} "
            f"(expected {self.expected:.1f}, ratio {self.ratio:.2f}x, "
            f"p = {self.p_value:.4f})"
        )


def binomial_enrichment(k: int, n: int, null_rate: float) -> BinomialResult:
    """
    One-tailed binomial test: is the observed count higher than expected?

    This is the core statistical test in the project. It asks:
    "If each site had a `null_rate` probability of being A+, what is the
    chance of seeing `k` or more A+ sites out of `n`?"

    Parameters
    ----------
    k : int
        Number of successes (e.g. A+ sites).
    n : int
        Total trials (e.g. total sites in population).
    null_rate : float
        Expected probability under H0 (e.g. P_NULL_AP from lib.beru for A+ tier).

    Returns
    -------
    BinomialResult
    """
    expected = n * null_rate
    ratio = (k / expected) if expected > 0 else float("inf")
    result = binomtest(k, n, null_rate, alternative="greater")
    return BinomialResult(
        observed=k,
        total=n,
        null_rate=null_rate,
        expected=expected,
        ratio=ratio,
        p_value=result.pvalue,
    )


def binomial_depletion(k: int, n: int, null_rate: float) -> BinomialResult:
    """
    One-tailed binomial test: is the observed count LOWER than expected?

    Mirror of binomial_enrichment for the depletion direction.
    Asks: "If each site had a `null_rate` probability of being A+, what
    is the chance of seeing `k` or fewer A+ sites out of `n`?"

    A significant result (p < 0.05) means the population is significantly
    *depleted* relative to the geometric null — the signature of a region
    where sites are structurally avoiding beru harmonics.

    Parameters
    ----------
    k : int
        Number of successes (e.g. A+ sites observed).
    n : int
        Total trials (e.g. total sites in the corpus).
    null_rate : float
        Expected probability under H0 (e.g. P_NULL_AP from lib.beru for A+ tier).

    Returns
    -------
    BinomialResult
        ratio < 1.0 indicates depletion; p_value is one-tailed (less).
    """
    expected = n * null_rate
    ratio = (k / expected) if expected > 0 else 0.0
    result = binomtest(k, n, null_rate, alternative="less")
    return BinomialResult(
        observed=k,
        total=n,
        null_rate=null_rate,
        expected=expected,
        ratio=ratio,
        p_value=result.pvalue,
    )


# ---------------------------------------------------------------------------
# Chi-square uniformity test
# ---------------------------------------------------------------------------

@dataclass
class ChiSquareResult:
    """Result of a chi-square goodness-of-fit test."""
    observed_bins: List[int]
    expected_per_bin: float
    statistic: float
    p_value: float
    n_bins: int


def chi_square_uniform(deviations: List[float], n_bins: int = 5,
                       max_dev: float = 0.05) -> ChiSquareResult:
    """
    Chi-square test for uniformity of beru deviations.

    Divides the deviation range [0, max_dev] into equal bins and tests
    whether sites are distributed uniformly across them.

    A low p-value means the deviations are NOT uniformly distributed,
    which supports the hypothesis of clustering near harmonics.

    Parameters
    ----------
    deviations : list of float
        Beru deviations (one per site).
    n_bins : int
        Number of equal-width bins (default 5).
    max_dev : float
        Maximum deviation to consider (0.05 beru = half a step).

    Returns
    -------
    ChiSquareResult
    """
    bin_width = max_dev / n_bins
    observed = [0] * n_bins
    for d in deviations:
        idx = min(int(d / bin_width), n_bins - 1)
        observed[idx] += 1

    n = len(deviations)
    expected_per_bin = n / n_bins
    stat, p = chisquare(observed, f_exp=[expected_per_bin] * n_bins)
    return ChiSquareResult(
        observed_bins=observed,
        expected_per_bin=expected_per_bin,
        statistic=stat,
        p_value=p,
        n_bins=n_bins,
    )


# ---------------------------------------------------------------------------
# Fisher exact test (2×2 contingency)
# ---------------------------------------------------------------------------

@dataclass
class FisherResult:
    """Result of a one-tailed Fisher exact test."""
    table: List[List[int]]
    odds_ratio: float
    p_value: float


def fisher_enrichment(a: int, b: int, c: int, d: int,
                      alternative: str = "greater") -> FisherResult:
    """
    Fisher exact test on a 2×2 contingency table.

    Table layout:
        [[a, b],    # group 1: hits, misses
         [c, d]]    # group 2: hits, misses

    Used for comparing A+ rates between sub-populations
    (e.g. founding vs non-founding, early vs late).

    Parameters
    ----------
    a, b, c, d : int
        Cell counts of the 2×2 table.
    alternative : str
        'greater', 'less', or 'two-sided'.

    Returns
    -------
    FisherResult
    """
    table = [[a, b], [c, d]]
    odds, p = fisher_exact(table, alternative=alternative)
    return FisherResult(table=table, odds_ratio=odds, p_value=p)


# ---------------------------------------------------------------------------
# Cochran-Armitage trend test
# ---------------------------------------------------------------------------

@dataclass
class CochranArmitageResult:
    """Result of a Cochran-Armitage trend test."""
    z_statistic: float
    p_value: float
    direction: str  # "decreasing" or "increasing" or "none"


def cochran_armitage(ns: List[int], successes: List[int],
                     scores: List[float] = None) -> CochranArmitageResult:
    """
    Cochran-Armitage trend test for ordered proportions.

    Tests whether the success rate changes monotonically across
    ordered groups (e.g. inscription-year cohorts).

    Parameters
    ----------
    ns : list of int
        Sample sizes for each group.
    successes : list of int
        Number of successes in each group.
    scores : list of float, optional
        Ordinal scores for each group (default: 1, 2, ..., k).

    Returns
    -------
    CochranArmitageResult
    """
    import numpy as np
    ns = np.asarray(ns, dtype=float)
    aps = np.asarray(successes, dtype=float)
    if scores is None:
        scores = np.arange(1, len(ns) + 1, dtype=float)
    else:
        scores = np.asarray(scores, dtype=float)

    N_total = ns.sum()
    R_total = aps.sum()
    p_bar = R_total / N_total
    s_bar = np.sum(scores * ns) / N_total

    T_num = np.sum(scores * aps) - R_total * s_bar
    T_den_sq = p_bar * (1 - p_bar) * (np.sum(scores**2 * ns) - N_total * s_bar**2)
    T_den = np.sqrt(T_den_sq) if T_den_sq > 0 else 0
    Z = T_num / T_den if T_den > 0 else 0.0

    # One-sided: Z < 0 means decreasing trend
    p = float(norm.cdf(Z))

    if Z < -0.01:
        direction = "decreasing"
    elif Z > 0.01:
        direction = "increasing"
    else:
        direction = "none"

    return CochranArmitageResult(z_statistic=Z, p_value=p, direction=direction)


# ---------------------------------------------------------------------------
# Spearman correlation
# ---------------------------------------------------------------------------

@dataclass
class SpearmanResult:
    """Result of a Spearman rank correlation."""
    rho: float
    p_value: float


def spearman_correlation(x: List[float], y: List[float]) -> SpearmanResult:
    """
    Spearman rank correlation between two variables.

    Used to test whether inscription year correlates with beru deviation
    (i.e. earlier sites tend to be closer to harmonics).

    Returns
    -------
    SpearmanResult
    """
    rho, p = spearmanr(x, y)
    return SpearmanResult(rho=float(rho), p_value=float(p))


# ---------------------------------------------------------------------------
# Two-proportion z-test
# ---------------------------------------------------------------------------

@dataclass
class TwoPropZResult:
    """Result of a one-sided two-proportion z-test."""
    p1: float
    p2: float
    z_statistic: float
    p_value: float


def two_prop_z(k1: int, n1: int, k2: int, n2: int) -> TwoPropZResult:
    """
    One-sided z-test: is proportion 1 greater than proportion 2?

    Uses a pooled standard error under H0: p1 = p2.

    Parameters
    ----------
    k1, n1 : int
        Successes and trials in group 1.
    k2, n2 : int
        Successes and trials in group 2.

    Returns
    -------
    TwoPropZResult
    """
    import numpy as np
    p1     = k1 / n1
    p2     = k2 / n2
    p_pool = (k1 + k2) / (n1 + n2)
    se     = np.sqrt(p_pool * (1 - p_pool) * (1 / n1 + 1 / n2))
    z      = (p1 - p2) / se if se > 0 else 0.0
    p      = float(norm.sf(z))
    return TwoPropZResult(p1=p1, p2=p2, z_statistic=float(z), p_value=p)


# ---------------------------------------------------------------------------
# Mantel-Haenszel stratified test
# ---------------------------------------------------------------------------

@dataclass
class MantelHaenszelResult:
    """Result of a Cochran-Mantel-Haenszel test across strata."""
    common_odds_ratio: float
    chi2_statistic: float
    p_value: float
    n_strata: int


def mantel_haenszel(
    strata: List[Tuple[int, int, int, int]]
) -> MantelHaenszelResult:
    """
    Cochran-Mantel-Haenszel test for a common odds ratio across strata.

    Each stratum is a 2x2 table represented as (k1, n1, k2, n2):
        group 1: k1 successes out of n1 trials
        group 2: k2 successes out of n2 trials

    Used to test whether group-1 enrichment persists after stratifying
    by a confound (e.g. regional composition).

    Parameters
    ----------
    strata : list of (k1, n1, k2, n2) tuples

    Returns
    -------
    MantelHaenszelResult
    """
    import numpy as np
    num = den = 0.0
    chi_num = chi_den = 0.0
    for (k1, n1, k2, n2) in strata:
        N = n1 + n2
        if N < 2:
            continue
        a = k1; b = n1 - k1; c = k2; d = n2 - k2
        num     += a * d / N
        den     += b * c / N
        E_a      = n1 * (k1 + k2) / N
        V_a      = (n1 * n2 * (k1 + k2) * (N - k1 - k2) / (N ** 2 * (N - 1))
                    if N > 1 else 0)
        chi_num += (a - E_a)
        chi_den += V_a
    common_or = num / den if den else float("nan")
    chi2      = chi_num ** 2 / chi_den if chi_den else float("nan")
    p         = float(scipy_stats.chi2.sf(chi2, df=1)) if not np.isnan(chi2) else float("nan")
    return MantelHaenszelResult(
        common_odds_ratio=common_or,
        chi2_statistic=chi2,
        p_value=p,
        n_strata=len(strata),
    )


# ---------------------------------------------------------------------------
# Significance helpers
# ---------------------------------------------------------------------------

def significance_label(p: float) -> str:
    """
    Human-readable significance label for a p-value.

    Returns:
        '***' if p < 0.001
        '**'  if p < 0.01
        '*'   if p < 0.05
        '~'   if p < 0.10
        'ns'  otherwise
    """
    if p < 0.001:
        return "***"
    if p < 0.01:
        return "**"
    if p < 0.05:
        return "*"
    if p < 0.10:
        return "~"
    return "ns"


def bonferroni(p: float, k: int) -> float:
    """
    Bonferroni correction: multiply p-value by the number of tests.

    Parameters
    ----------
    p : float
        Raw p-value from a single test.
    k : int
        Number of independent tests performed.

    Returns
    -------
    float
        Adjusted p-value, capped at 1.0.
    """
    return min(p * k, 1.0)


def null_rate_at_spacing(threshold: float, spacing: float) -> float:
    """
    Geometric null rate for a given threshold and harmonic spacing.

    Each harmonic has a window of 2 × threshold beru. Harmonics repeat
    every `spacing` beru. The null fraction is (2 × threshold) / spacing.

    Parameters
    ----------
    threshold : float
        Tier threshold in beru (e.g. 0.002 for A+).
    spacing : float
        Harmonic spacing in beru (e.g. 0.10).

    Returns
    -------
    float
        Null probability under uniform assumption.
    """
    return 2 * threshold / spacing
