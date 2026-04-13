"""
dome_filter.py
==============
Shared context-aware dome/stupa/tholos keyword matching logic.

Used by:
  - spherical_monument_test.py (primary analysis)
  - simulation_null_model.py (permutation tests)
  - unit_sweep_fill.py (unit sweep tables)

The matching system has two stages:
  1. KEYWORD MATCH: Find morphological keywords in full UNESCO text
  2. CONTEXT VALIDATION: For ambiguous keywords (dome/domed/domes/spherical),
     validate that the keyword appears in a sentence describing built
     monumental architecture, not natural geology, vernacular housing,
     or non-architectural objects.

Unambiguous keywords (stupa, stupas, tholos) are accepted without
context validation — these terms have no meaning outside monumental
architecture.

All keyword lists and context-validation patterns are loaded at import
time from keywords.json (project root). No strings or regexes are
hardcoded in this file — edit keywords.json to add, remove, or change
any keyword or pattern.
"""

import re
import json
from pathlib import Path

# ── Load ALL config from keywords.json ───────────────────────────────────────
_KW_PATH = Path(__file__).parent.parent / "keywords.json"
with open(_KW_PATH) as _f:
    _KW = json.load(_f)

_dome_cfg = _KW["dome_forms"]

UNAMBIGUOUS_KEYWORDS = _dome_cfg["unambiguous"]
AMBIGUOUS_KEYWORDS   = _dome_cfg["ambiguous"]
FORM_KEYWORDS        = UNAMBIGUOUS_KEYWORDS + AMBIGUOUS_KEYWORDS

FORM_KEYWORD_RES = {
    kw: re.compile(r"\b" + re.escape(kw) + r"\b", re.IGNORECASE)
    for kw in FORM_KEYWORDS
}

# ── Negative-context patterns ────────────────────────────────────────────────
# If any of these appear in the SAME SENTENCE as an ambiguous keyword,
# the sentence is disqualified.
NEGATIVE_CONTEXT     = _dome_cfg["negative_context"]
NEGATIVE_CONTEXT_RES = [re.compile(pat, re.IGNORECASE) for pat in NEGATIVE_CONTEXT]

# ── Positive-context patterns ────────────────────────────────────────────────
# For an ambiguous keyword to count, the sentence must contain at least
# one architectural-context term.
POSITIVE_CONTEXT     = _dome_cfg["positive_context"]
POSITIVE_CONTEXT_RES = [re.compile(pat, re.IGNORECASE) for pat in POSITIVE_CONTEXT]


def validate_ambiguous_sentence(sentence: str) -> bool:
    """Check if a sentence with an ambiguous keyword describes
    monumental architecture (True) or something else (False).

    Rules:
      1. If any NEGATIVE_CONTEXT pattern matches → False
      2. If any POSITIVE_CONTEXT pattern matches → True
      3. Otherwise → False (conservative: reject unclear contexts)
    """
    for pat in NEGATIVE_CONTEXT_RES:
        if pat.search(sentence):
            return False
    for pat in POSITIVE_CONTEXT_RES:
        if pat.search(sentence):
            return True
    return False


def rejection_reason(sentence: str) -> str:
    """Return a human-readable reason why a sentence was rejected.

    Returns the first matching negative-context pattern, or
    'no-positive-context' if no negative pattern matched but no
    positive pattern matched either.
    """
    for pat in NEGATIVE_CONTEXT_RES:
        m = pat.search(sentence)
        if m:
            return f"negative-context: «{m.group(0)}» matched /{pat.pattern}/"
    return "no-positive-context (no architectural term found)"

# Keep private alias for internal use
_validate_ambiguous_sentence = validate_ambiguous_sentence


def validate_keyword_match(full_text: str, keyword: str) -> tuple:
    """Validate a keyword match in context.

    For unambiguous keywords (stupa, stupas, tholos): always valid.
    For ambiguous keywords (dome, domed, domes, spherical): valid only if
    at least one matching sentence passes context validation.

    Returns (is_valid, matched_sentences, validation_notes)
    """
    if keyword in UNAMBIGUOUS_KEYWORDS:
        return True, [], "unambiguous"

    sentences = re.split(r"(?<=[.!?])\s+", full_text)
    kw_re = FORM_KEYWORD_RES[keyword]

    valid_sentences = []
    rejected_sentences = []
    for sent in sentences:
        if kw_re.search(sent):
            if validate_ambiguous_sentence(sent):
                valid_sentences.append(sent.strip())
            else:
                rejected_sentences.append(sent.strip())

    if valid_sentences:
        return True, valid_sentences, "context-validated"
    elif rejected_sentences:
        return False, rejected_sentences, "rejected-context"
    else:
        return False, [], "no-sentence-found"


def is_dome_site(site_obj) -> bool:
    """Determine whether a UNESCO site qualifies as a domed/stupa/tholos
    monument using context-aware keyword matching.

    Parameters
    ----------
    site_obj : object with .full_text, .short_description, .extended_description

    Returns
    -------
    bool : True if the site qualifies
    """
    full_text = site_obj.full_text  # lowercased already
    # Stage 1: raw keyword match
    raw_matched = [k for k in FORM_KEYWORDS if FORM_KEYWORD_RES[k].search(full_text)]
    if not raw_matched:
        return False

    # If any unambiguous keyword matched, accept immediately
    if any(k in UNAMBIGUOUS_KEYWORDS for k in raw_matched):
        return True

    # Stage 2: context-validate ambiguous keywords
    all_text = (site_obj.short_description or "")
    if site_obj.extended_description:
        all_text += " " + site_obj.extended_description

    sentences = re.split(r"(?<=[.!?])\s+", all_text)
    for kw in raw_matched:
        kw_re = FORM_KEYWORD_RES[kw]
        for sent in sentences:
            if kw_re.search(sent) and _validate_ambiguous_sentence(sent):
                return True

    return False
