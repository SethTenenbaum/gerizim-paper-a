"""
Unit tests for the statistical test wrappers.
"""

import pytest
from lib.stats import (
    binomial_enrichment, BinomialResult,
    chi_square_uniform, ChiSquareResult,
    fisher_enrichment, FisherResult,
    cochran_armitage, CochranArmitageResult,
    spearman_correlation, SpearmanResult,
    significance_label, bonferroni, null_rate_at_spacing,
)


# ---------------------------------------------------------------------------
# Binomial enrichment
# ---------------------------------------------------------------------------

class TestBinomialEnrichment:

    def test_perfect_match(self):
        """If observed == expected, p should be around 0.5 (not significant)."""
        result = binomial_enrichment(4, 100, 0.04)
        assert result.observed == 4
        assert result.total == 100
        assert result.expected == pytest.approx(4.0)
        assert result.ratio == pytest.approx(1.0)
        assert result.p_value > 0.3  # not significant

    def test_strong_enrichment(self):
        """Observing many more than expected should give small p."""
        result = binomial_enrichment(20, 100, 0.04)
        assert result.ratio > 4.0
        assert result.p_value < 0.001

    def test_zero_observed(self):
        """Zero observed should have p close to 1."""
        result = binomial_enrichment(0, 100, 0.04)
        assert result.p_value > 0.99

    def test_observed_rate_property(self):
        result = binomial_enrichment(10, 200, 0.04)
        assert result.observed_rate == pytest.approx(0.05)

    def test_summary_string(self):
        result = binomial_enrichment(10, 100, 0.04)
        s = result.summary()
        assert "10/100" in s
        assert "p =" in s


# ---------------------------------------------------------------------------
# Chi-square uniformity
# ---------------------------------------------------------------------------

class TestChiSquareUniform:

    def test_uniform_distribution(self):
        """Perfectly uniform deviations should give a high p-value."""
        # 100 sites evenly spaced from 0 to 0.05
        deviations = [i * 0.05 / 100 for i in range(100)]
        result = chi_square_uniform(deviations)
        assert result.p_value > 0.05

    def test_clustered_distribution(self):
        """All sites near zero should give a low p-value."""
        deviations = [0.001] * 50 + [0.04] * 10
        result = chi_square_uniform(deviations)
        assert result.p_value < 0.01

    def test_bin_count_sum(self):
        """Sum of bin counts should equal total sites."""
        deviations = [0.005, 0.015, 0.025, 0.035, 0.045]
        result = chi_square_uniform(deviations)
        assert sum(result.observed_bins) == 5

    def test_custom_bins(self):
        result = chi_square_uniform([0.01] * 10, n_bins=3)
        assert result.n_bins == 3
        assert len(result.observed_bins) == 3


# ---------------------------------------------------------------------------
# Fisher exact test
# ---------------------------------------------------------------------------

class TestFisherEnrichment:

    def test_strong_enrichment(self):
        """Group 1 has much higher hit rate than group 2."""
        result = fisher_enrichment(a=20, b=5, c=5, d=20)
        assert result.odds_ratio > 1.0
        assert result.p_value < 0.001

    def test_no_enrichment(self):
        """Equal hit rates should give OR ≈ 1."""
        result = fisher_enrichment(a=10, b=10, c=10, d=10)
        assert result.odds_ratio == pytest.approx(1.0)
        assert result.p_value > 0.3

    def test_table_stored(self):
        result = fisher_enrichment(a=5, b=10, c=3, d=12)
        assert result.table == [[5, 10], [3, 12]]


# ---------------------------------------------------------------------------
# Cochran-Armitage trend test
# ---------------------------------------------------------------------------

class TestCochranArmitage:

    def test_decreasing_trend(self):
        """Clear decreasing rates should give negative Z."""
        ns = [100, 100, 100]
        successes = [20, 10, 2]
        result = cochran_armitage(ns, successes)
        assert result.z_statistic < 0
        assert result.p_value < 0.05
        assert result.direction == "decreasing"

    def test_no_trend(self):
        """Equal rates across groups should give Z near 0."""
        ns = [100, 100, 100]
        successes = [10, 10, 10]
        result = cochran_armitage(ns, successes)
        assert abs(result.z_statistic) < 1.0
        assert result.p_value > 0.1

    def test_custom_scores(self):
        result = cochran_armitage([50, 50], [10, 5], scores=[1.0, 2.0])
        assert isinstance(result, CochranArmitageResult)


# ---------------------------------------------------------------------------
# Spearman correlation
# ---------------------------------------------------------------------------

class TestSpearmanCorrelation:

    def test_perfect_positive(self):
        result = spearman_correlation([1, 2, 3, 4, 5], [10, 20, 30, 40, 50])
        assert result.rho == pytest.approx(1.0)
        assert result.p_value < 0.01

    def test_perfect_negative(self):
        result = spearman_correlation([1, 2, 3, 4, 5], [50, 40, 30, 20, 10])
        assert result.rho == pytest.approx(-1.0)
        assert result.p_value < 0.01

    def test_no_correlation(self):
        # Alternating pattern — weak or no monotonic trend
        result = spearman_correlation([1, 2, 3, 4], [10, 30, 20, 40])
        assert abs(result.rho) < 1.0


# ---------------------------------------------------------------------------
# Significance helpers
# ---------------------------------------------------------------------------

class TestSignificanceLabel:

    def test_three_stars(self):
        assert significance_label(0.0001) == "***"

    def test_two_stars(self):
        assert significance_label(0.005) == "**"

    def test_one_star(self):
        assert significance_label(0.03) == "*"

    def test_marginal(self):
        assert significance_label(0.07) == "~"

    def test_not_significant(self):
        assert significance_label(0.5) == "ns"

    def test_boundary_001(self):
        """Exactly 0.001 is NOT < 0.001, so it falls into the ** bucket."""
        assert significance_label(0.001) == "**"

    def test_boundary_01(self):
        """Exactly 0.01 is NOT < 0.01, so it falls into the * bucket."""
        assert significance_label(0.01) == "*"

    def test_boundary_05(self):
        """Exactly 0.05 is NOT < 0.05, so it falls into the ~ bucket."""
        assert significance_label(0.05) == "~"


class TestBonferroni:

    def test_basic(self):
        assert bonferroni(0.01, 5) == pytest.approx(0.05)

    def test_capped_at_one(self):
        assert bonferroni(0.5, 5) == 1.0

    def test_single_test(self):
        assert bonferroni(0.03, 1) == pytest.approx(0.03)


class TestNullRateAtSpacing:

    def test_standard_aplus(self):
        """At 0.10 spacing, A+ null rate = 2×0.002/0.1 = 0.04."""
        rate = null_rate_at_spacing(0.002, 0.10)
        assert rate == pytest.approx(0.04)

    def test_wider_spacing(self):
        """At 0.20 spacing, A+ null rate = 2×0.002/0.2 = 0.02."""
        rate = null_rate_at_spacing(0.002, 0.20)
        assert rate == pytest.approx(0.02)
