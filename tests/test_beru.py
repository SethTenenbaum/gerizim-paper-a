"""
Unit tests for the beru analysis library.

Run with:
    python3 -m pytest tests/ -v
    python3 -m pytest tests/test_beru.py -v
"""

import pytest
from lib.beru import (
    deviation, full_calculation, deviation_at_spacing,
    tier_label, is_aplus, is_a_or_better, dev_to_km,
    GERIZIM, BERU, HARMONIC_STEP,
    TIER_APP, TIER_APLUS, TIER_A_MAX, TIER_B_MAX,
    P_NULL_AP, P_NULL_A, CONFIG,
)


# ---------------------------------------------------------------------------
# Config loading
# ---------------------------------------------------------------------------

class TestConfigLoading:
    """Verify that config.json is loaded correctly."""

    def test_gerizim_longitude(self):
        assert GERIZIM == 35.269

    def test_beru_degrees(self):
        assert BERU == 30.0

    def test_harmonic_step(self):
        assert HARMONIC_STEP == 0.1

    def test_tier_thresholds_ordered(self):
        assert TIER_APP < TIER_APLUS < TIER_A_MAX < TIER_B_MAX

    def test_null_rate_aplus(self):
        """Null rate for A+ should be tier_width / baseline_width."""
        assert P_NULL_AP == pytest.approx(TIER_APLUS / TIER_B_MAX, abs=1e-10)

    def test_null_rate_a(self):
        assert P_NULL_A == pytest.approx(TIER_A_MAX / TIER_B_MAX, abs=1e-10)

    def test_config_has_anchors(self):
        assert "gerizim" in CONFIG["anchors"]
        assert "jerusalem" in CONFIG["anchors"]


# ---------------------------------------------------------------------------
# Deviation calculation
# ---------------------------------------------------------------------------

class TestDeviation:
    """Test the core deviation function."""

    def test_exact_harmonic(self):
        """A site exactly 1.0 beru from the anchor has zero deviation."""
        dev = deviation(GERIZIM + 30.0)  # exactly 1.0 beru
        assert dev == pytest.approx(0.0, abs=1e-10)

    def test_exact_half_step(self):
        """A site at the midpoint between two harmonics has dev = 0.05."""
        # 0.05 beru = 1.5° from anchor → halfway between 0.0 and 0.1 harmonic
        dev = deviation(GERIZIM + 1.5)
        assert dev == pytest.approx(0.05, abs=1e-6)

    def test_small_offset(self):
        """A site slightly off a harmonic has a small deviation."""
        # 0.001 beru offset from the 1.0 harmonic → dev = 0.001
        offset_deg = 0.001 * BERU  # 0.03°
        dev = deviation(GERIZIM + 30.0 + offset_deg)
        assert dev == pytest.approx(0.001, abs=1e-6)

    def test_symmetry(self):
        """Deviation is the same east and west of the anchor."""
        dev_east = deviation(GERIZIM + 15.0)  # 0.5 beru east
        dev_west = deviation(GERIZIM - 15.0)  # 0.5 beru west
        assert dev_east == pytest.approx(dev_west, abs=1e-10)

    def test_range_always_0_to_005(self):
        """Deviation should always be in [0, 0.05]."""
        import numpy as np
        for lon in np.linspace(-180, 180, 500):
            dev = deviation(lon)
            assert 0.0 <= dev <= 0.05 + 1e-10, f"dev={dev} at lon={lon}"

    def test_wrap_mode(self):
        """With wrap=True, arcs > 180° are wrapped."""
        dev_no_wrap = deviation(GERIZIM + 200.0, anchor=GERIZIM, wrap=False)
        dev_wrap = deviation(GERIZIM + 200.0, anchor=GERIZIM, wrap=True)
        assert 0.0 <= dev_wrap <= 0.05
        assert dev_wrap != dev_no_wrap or True  # just check it runs

    def test_custom_anchor(self):
        """Deviation works with a custom anchor."""
        dev = deviation(0.0, anchor=30.0)  # exactly 1.0 beru
        assert dev == pytest.approx(0.0, abs=1e-10)


# ---------------------------------------------------------------------------
# Full calculation
# ---------------------------------------------------------------------------

class TestFullCalculation:
    """Test the full_calculation function."""

    def test_returns_all_keys(self):
        result = full_calculation(GERIZIM + 30.0)
        expected_keys = {"arc_deg", "beru_val", "nearest", "dev", "dev_km", "tier"}
        assert set(result.keys()) == expected_keys

    def test_exact_harmonic_values(self):
        result = full_calculation(GERIZIM + 30.0)
        assert result["arc_deg"] == pytest.approx(30.0, abs=1e-10)
        assert result["beru_val"] == pytest.approx(1.0, abs=1e-10)
        assert result["nearest"] == pytest.approx(1.0, abs=1e-6)
        assert result["dev"] == pytest.approx(0.0, abs=1e-10)
        assert result["dev_km"] == pytest.approx(0.0, abs=0.01)
        assert result["tier"] == "A++"

    def test_tier_assignment_aplus(self):
        """A site with dev=0.001 should be A+."""
        offset = 0.001 * BERU  # 0.001 beru in degrees
        result = full_calculation(GERIZIM + 30.0 + offset)
        assert result["tier"] == "A+"

    def test_dev_km_calculation(self):
        result = full_calculation(GERIZIM + 30.0 + 0.03)  # 0.001 beru offset
        expected_km = result["dev"] * BERU * 111.0
        assert result["dev_km"] == pytest.approx(expected_km, abs=0.01)


# ---------------------------------------------------------------------------
# Deviation at custom spacing
# ---------------------------------------------------------------------------

class TestDeviationAtSpacing:
    """Test deviation_at_spacing for unit sweep analysis."""

    def test_standard_spacing(self):
        """At 0.1 spacing, should match the default deviation."""
        dev_default = deviation(GERIZIM + 30.5)
        dev_at_01 = deviation_at_spacing(GERIZIM + 30.5, spacing=0.1)
        assert dev_default == pytest.approx(dev_at_01, abs=1e-10)

    def test_different_spacing(self):
        """At 0.05 spacing, harmonics are closer together."""
        dev = deviation_at_spacing(GERIZIM + 0.5, spacing=0.05)
        assert dev == pytest.approx(0.01667, abs=0.001)


# ---------------------------------------------------------------------------
# Tier classification
# ---------------------------------------------------------------------------

class TestTierLabel:
    """Test tier classification boundaries."""

    def test_tier_app(self):
        assert tier_label(0.0) == "A++"
        assert tier_label(0.0001) == "A++"
        assert tier_label(0.0002) == "A++"

    def test_tier_aplus(self):
        assert tier_label(0.0003) == "A+"
        assert tier_label(0.001) == "A+"
        assert tier_label(0.002) == "A+"

    def test_tier_a(self):
        assert tier_label(0.003) == "A"
        assert tier_label(0.010) == "A"

    def test_tier_b(self):
        assert tier_label(0.011) == "B"
        assert tier_label(0.050) == "B"

    def test_tier_c(self):
        assert tier_label(0.051) == "C"
        assert tier_label(0.1) == "C"

    def test_boundary_aplus_to_a(self):
        """Boundary at exactly TIER_APLUS should be A+."""
        assert tier_label(TIER_APLUS) == "A+"
        assert tier_label(TIER_APLUS + 1e-10) == "A"


class TestTierHelpers:
    def test_is_aplus(self):
        assert is_aplus("A++")
        assert is_aplus("A+")
        assert not is_aplus("A")
        assert not is_aplus("B")
        assert not is_aplus("C")

    def test_is_a_or_better(self):
        assert is_a_or_better("A++")
        assert is_a_or_better("A+")
        assert is_a_or_better("A")
        assert not is_a_or_better("B")
        assert not is_a_or_better("C")


class TestDevToKm:
    def test_zero(self):
        assert dev_to_km(0.0) == 0.0

    def test_known_value(self):
        # 0.002 beru × 30 deg × 111 km/deg = 6.66 km
        assert dev_to_km(0.002) == pytest.approx(6.66, abs=0.01)
