"""
Unit tests for the dome keyword filter module.
"""

import pytest
from lib.dome_filter import (
    UNAMBIGUOUS_KEYWORDS, AMBIGUOUS_KEYWORDS, FORM_KEYWORDS,
    FORM_KEYWORD_RES,
    validate_ambiguous_sentence, validate_keyword_match,
)


class TestKeywordSets:

    def test_form_keywords_is_union(self):
        """FORM_KEYWORDS should be unambiguous + ambiguous."""
        assert set(FORM_KEYWORDS) == set(UNAMBIGUOUS_KEYWORDS + AMBIGUOUS_KEYWORDS)

    def test_unambiguous_keywords(self):
        assert "stupa" in UNAMBIGUOUS_KEYWORDS
        assert "tholos" in UNAMBIGUOUS_KEYWORDS

    def test_ambiguous_keywords(self):
        assert "dome" in AMBIGUOUS_KEYWORDS
        assert "domed" in AMBIGUOUS_KEYWORDS
        assert "spherical" in AMBIGUOUS_KEYWORDS

    def test_regex_patterns_exist(self):
        for kw in FORM_KEYWORDS:
            assert kw in FORM_KEYWORD_RES


class TestKeywordRegex:

    def test_word_boundary_match(self):
        """Keywords should match at word boundaries."""
        assert FORM_KEYWORD_RES["stupa"].search("the great stupa of Sanchi")
        assert FORM_KEYWORD_RES["dome"].search("a magnificent dome")

    def test_no_partial_match(self):
        """Keywords should NOT match partial words."""
        assert not FORM_KEYWORD_RES["dome"].search("the kingdom was powerful")
        assert not FORM_KEYWORD_RES["stupa"].search("the stupambitious plan")

    def test_case_insensitive(self):
        assert FORM_KEYWORD_RES["dome"].search("The Dome of the Rock")
        assert FORM_KEYWORD_RES["stupa"].search("STUPA architecture")


class TestAmbiguousSentenceValidation:

    def test_architectural_context_accepted(self):
        """A sentence with 'dome' + 'mosque' should be accepted."""
        sent = "The mosque features a large dome covered in tiles."
        assert validate_ambiguous_sentence(sent) is True

    def test_geological_context_rejected(self):
        """A sentence with 'dome' + 'karst' should be rejected."""
        sent = "The karst dome rises above the limestone valley."
        assert validate_ambiguous_sentence(sent) is False

    def test_vernacular_context_rejected(self):
        """A sentence with 'dome' + 'hut' should be rejected."""
        sent = "The traditional hut has a dome-shaped roof."
        assert validate_ambiguous_sentence(sent) is False

    def test_atomic_bomb_context_rejected(self):
        """The Hiroshima Peace Memorial dome should be rejected."""
        sent = "The atomic bomb dome is a peace memorial."
        assert validate_ambiguous_sentence(sent) is False

    def test_no_context_rejected(self):
        """A sentence with 'dome' but no context words should be rejected."""
        sent = "The structure has a dome."
        assert validate_ambiguous_sentence(sent) is False

    def test_church_context(self):
        sent = "The church is crowned by a dome decorated with mosaics."
        assert validate_ambiguous_sentence(sent) is True

    def test_byzantine_context(self):
        sent = "The domed Byzantine basilica dates from the 6th century."
        assert validate_ambiguous_sentence(sent) is True


class TestValidateKeywordMatch:

    def test_unambiguous_always_valid(self):
        """Unambiguous keywords are always valid, no context needed."""
        is_valid, _, note = validate_keyword_match("Some random text", "stupa")
        assert is_valid is True
        assert note == "unambiguous"

    def test_unambiguous_tholos(self):
        is_valid, _, note = validate_keyword_match("No architecture here", "tholos")
        assert is_valid is True

    def test_ambiguous_with_good_context(self):
        text = "The mosque features a magnificent dome. It was built in 1453."
        is_valid, sentences, note = validate_keyword_match(text, "dome")
        assert is_valid is True
        assert len(sentences) > 0

    def test_ambiguous_with_bad_context(self):
        text = "The karst dome rises above the geological formation."
        is_valid, sentences, note = validate_keyword_match(text, "dome")
        assert is_valid is False

    def test_ambiguous_no_sentence(self):
        """If the keyword doesn't actually appear, it should not validate."""
        is_valid, _, note = validate_keyword_match("No keywords here.", "dome")
        assert is_valid is False
        assert note == "no-sentence-found"
