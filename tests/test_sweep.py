"""
Unit tests for the anchor sweep module.
"""

import pytest
import numpy as np
from lib.sweep import count_hits, run_sweep, percentile_rank, summarize_anchor
from lib.beru import GERIZIM, TIER_APLUS, TIER_A_MAX


class TestCountHits:

    def test_exact_harmonic_sites(self):
        """Sites exactly on harmonics should all be counted."""
        # Sites at exactly 0.0, 1.0, 2.0 beru from Gerizim
        lons = [35.272, 65.272, 95.272]
        count = count_hits(35.272, lons, TIER_APLUS)
        assert count == 3

    def test_no_hits(self):
        """Sites far from any harmonic should not be counted."""
        # 1.5° from anchor = 0.05 beru = exactly at midpoint → dev=0.05
        lons = [35.272 + 1.5]
        count = count_hits(35.272, lons, TIER_APLUS)
        assert count == 0

    def test_empty_list(self):
        count = count_hits(35.272, [], TIER_APLUS)
        assert count == 0

    def test_different_thresholds(self):
        """A+ threshold is stricter than A threshold."""
        lons = [35.272 + 0.15]  # 0.15° = 0.005 beru → A but not A+
        count_ap = count_hits(35.272, lons, TIER_APLUS)
        count_a = count_hits(35.272, lons, TIER_A_MAX)
        assert count_ap == 0  # dev=0.005, threshold=0.002
        assert count_a == 1   # dev=0.005, threshold=0.010


class TestRunSweep:

    def test_returns_expected_keys(self):
        lons = [35.272, 65.272]
        result = run_sweep(lons, start=35.0, end=36.0, step=0.1)
        assert "sweep_anchors" in result
        assert "counts_aplus" in result
        assert "counts_a" in result

    def test_sweep_length(self):
        result = run_sweep([35.272], start=35.0, end=36.0, step=0.1)
        assert len(result["sweep_anchors"]) == len(result["counts_aplus"])

    def test_gerizim_anchor_in_sweep(self):
        """The Gerizim longitude should produce max hits for harmonic sites."""
        lons = [35.272 + 30.0, 35.272 + 60.0]  # exact harmonics
        result = run_sweep(lons, start=35.0, end=35.5, step=0.01)
        # At anchor=35.272, both sites should be hits
        idx = int(round((35.272 - 35.0) / 0.01))
        assert result["counts_aplus"][idx] >= 2


class TestPercentileRank:

    def test_all_below(self):
        """If value is the maximum, percentile should be 100."""
        counts = np.array([1, 2, 3, 4, 5])
        assert percentile_rank(counts, 5) == 100.0

    def test_all_above(self):
        """If value is below the minimum, percentile should be 0."""
        counts = np.array([5, 6, 7, 8])
        assert percentile_rank(counts, 2) == 0.0

    def test_median(self):
        """Value at the median should give ~50th percentile."""
        counts = np.array([1, 2, 3, 4, 5, 6, 7, 8, 9, 10])
        pctile = percentile_rank(counts, 5)
        assert 40 < pctile < 60


class TestSummariseAnchor:

    def test_returns_expected_keys(self):
        lons = [35.272, 65.272]
        sweep_result = run_sweep(lons, start=35.0, end=36.0, step=0.1)
        summary = summarize_anchor(lons, 35.272, sweep_result)
        assert "count_aplus" in summary
        assert "count_a" in summary
        assert "pctile_aplus" in summary
        assert "pctile_a" in summary
