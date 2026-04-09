"""
Integration tests that verify the config.json ↔ lib module pipeline.
"""

import json
import pytest
from pathlib import Path


class TestConfigIntegrity:
    """Verify config.json structure and that lib modules read it correctly."""

    @pytest.fixture
    def config(self):
        config_path = Path(__file__).parent.parent / "config.json"
        with open(config_path) as f:
            return json.load(f)

    def test_config_has_required_sections(self, config):
        for section in ["anchors", "units", "tiers", "null_rates", "anchor_sweep"]:
            assert section in config, f"Missing section: {section}"

    def test_gerizim_anchor(self, config):
        assert "gerizim" in config["anchors"]
        assert config["anchors"]["gerizim"]["longitude"] == 35.269

    def test_jerusalem_anchor(self, config):
        assert "jerusalem" in config["anchors"]
        assert isinstance(config["anchors"]["jerusalem"]["longitude"], (int, float))

    def test_beru_unit(self, config):
        assert config["units"]["beru"]["degrees"] == 30.0
        assert config["units"]["harmonic_step"] == 0.1

    def test_tiers_ordered(self, config):
        app = config["tiers"]["A++"]["max_deviation_beru"]
        ap = config["tiers"]["A+"]["max_deviation_beru"]
        a = config["tiers"]["A"]["max_deviation_beru"]
        b = config["tiers"]["B"]["max_deviation_beru"]
        assert app < ap < a < b

    def test_null_rates_consistent(self, config):
        """Null rates should equal tier_width / baseline_width."""
        ap_threshold = config["tiers"]["A+"]["max_deviation_beru"]
        a_threshold = config["tiers"]["A"]["max_deviation_beru"]
        b_threshold = config["tiers"]["B"]["max_deviation_beru"]
        assert config["null_rates"]["tier_aplus"] == pytest.approx(
            ap_threshold / b_threshold, abs=1e-10
        )
        assert config["null_rates"]["tier_a"] == pytest.approx(
            a_threshold / b_threshold, abs=1e-10
        )

    def test_lib_beru_matches_config(self, config):
        """lib.beru constants should match config.json values."""
        from lib.beru import GERIZIM, BERU, TIER_APLUS, TIER_A_MAX, TIER_B_MAX
        assert GERIZIM == config["anchors"]["gerizim"]["longitude"]
        assert BERU == config["units"]["beru"]["degrees"]
        assert TIER_APLUS == config["tiers"]["A+"]["max_deviation_beru"]
        assert TIER_A_MAX == config["tiers"]["A"]["max_deviation_beru"]
        assert TIER_B_MAX == config["tiers"]["B"]["max_deviation_beru"]

    def test_sweep_configs_present(self, config):
        assert "levant" in config["anchor_sweep"]
        assert "global" in config["anchor_sweep"]


class TestBackwardCompatibility:
    """Ensure data.scripts.unesco_corpus backward-compatible aliases still work."""

    def test_imports(self):
        from data.scripts.unesco_corpus import (
            GERIZIM, BERU, TIER_APLUS, P_NULL_AP,
            beru_deviation, tier_label, is_aplus, sig_label,
        )
        assert GERIZIM == 35.269
        assert callable(beru_deviation)
        assert callable(tier_label)

    def test_beru_deviation_alias(self):
        from data.scripts.unesco_corpus import beru_deviation
        from lib.beru import deviation
        # Both should give the same result
        lon = 65.274
        assert beru_deviation(lon) == pytest.approx(deviation(lon), abs=1e-10)

    def test_tier_label_alias(self):
        from data.scripts.unesco_corpus import tier_label as tl1
        from lib.beru import tier_label as tl2
        assert tl1(0.001) == tl2(0.001)
        assert tl1(0.05) == tl2(0.05)
