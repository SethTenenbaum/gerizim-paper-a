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
"""

import re
import json
from pathlib import Path

# ── Load keyword sets from config ────────────────────────────────────────────
_CONFIG_PATH = Path(__file__).parent.parent / "config.json"
with open(_CONFIG_PATH) as _f:
    _CONFIG = json.load(_f)

_dome_cfg = _CONFIG["keywords"]["dome_forms"]
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
NEGATIVE_CONTEXT = [
    r"\bkarst\b", r"\brock island", r"\bgeological\b",
    r"\bconglomerate\b", r"\bsandstone\b", r"\bbasalt\b", r"\berosion\b",
    r"\blimestone formation", r"\brock dome", r"\brock formation",
    r"\bnatural formation", r"\bnatural monument",
    # Vernacular / domestic housing
    r"\bhut\b", r"\bhuts\b", r"\bhutment\b",
    r"\bpastoral\b", r"\bpastoralist\b", r"\btranshumance\b",
    r"\bvernacular house", r"\btraditional .{0,15}house",
    r"\bmat\b.{0,20}\bbraid",
    # Non-architectural objects
    r"\bstone sphere", r"\bstone ball",
    # War memorials (dome as incidental name, not architectural form)
    r"\batomic bomb\b", r"\bfirst atomic\b", r"\bpeace memorial\b",
    r"\bgenbaku\b", r"\bhypocenter\b",
]
NEGATIVE_CONTEXT_RES = [re.compile(pat, re.IGNORECASE) for pat in NEGATIVE_CONTEXT]

# ── Positive-context patterns ────────────────────────────────────────────────
# For an ambiguous keyword to count, the sentence must contain at least
# one architectural-context term.
POSITIVE_CONTEXT = [
    # Building types (with plural/inflected forms)
    r"\bchurch\w*\b", r"\bmosque\w*\b", r"\bcathedral\w*\b", r"\btemple\w*\b",
    r"\bbasilica\w*\b", r"\bmausoleum\w*\b", r"\btomb\w*\b",
    r"\bpalace\w*\b", r"\bmonument\w*\b", r"\bbuilding\w*\b", r"\bhall\b",
    r"\bchapel\w*\b", r"\bobservator\w*\b", r"\blighthouse\w*\b", r"\bshrine\w*\b",
    r"\bmonaster\w*\b", r"\bmedresseh\w*\b", r"\bmadrasa\w*\b",
    r"\bhospice\w*\b", r"\bchamber\w*\b", r"\btower\w*\b",
    r"\bfortress\w*\b", r"\bcitadel\w*\b",
    # Architectural elements (with plural/inflected forms)
    r"\bminaret\w*\b", r"\bportal\w*\b", r"\bvault\w*\b",
    r"\bapse\w*\b", r"\bnave\w*\b", r"\bcupola\w*\b", r"\bdrum\b",
    r"\boctagonal\b", r"\brotunda\w*\b", r"\blantern\w*\b",
    r"\bmihrab\w*\b", r"\biwan\w*\b", r"\bsquinch\w*\b",
    r"\bbaptister\w*\b", r"\bsanctuar\w*\b", r"\bcorbel\w*\b",
    r"\bbelfr\w*\b", r"\breliqu\w*\b", r"\bportico\w*\b", r"\bcolumn\w*\b",
    r"\bcross.{0,5}plan\b", r"\bcross.{0,5}hall\b",
    r"\bribbed\b", r"\bfaience\b",
    r"\bspire\w*\b", r"\bceiling\w*\b", r"\broof\w*\b",
    # Construction / architecture terms
    r"\barchitect\w*\b", r"\bconstruct\w*\b", r"\bbuilt\b",
    r"\bmasonry\b", r"\bbrick\w*\b", r"\bstone\b",
    r"\bfounded\b", r"\berected\b", r"\bsurmounted\b",
    r"\broofed\b", r"\bcovered\b", r"\bcrowned\b",
    r"\bdecorat\w*\b", r"\bmosaic\w*\b", r"\btile\w*\b",
    r"\bfresco\w*\b", r"\bstatuar\w*\b", r"\bsculptur\w*\b",
    # Specific dome-architecture terms
    r"\bdouble.shelled\b", r"\bhemispheric\w*\b",
    r"\bonion dome\b", r"\bcross.{0,3}dome\b",
    # Period/style markers for architecture
    r"\bbyzantin\w*\b", r"\bromanesque\b", r"\bgothic\b",
    r"\brenaissance\b", r"\bbaroque\b", r"\bmamluk\b",
    r"\btimurid\b", r"\bseljuk\b", r"\bsassanid\b", r"\bottoman\b",
    r"\bislamic\b", r"\bmedieval\b",
    r"\bcaldarium\b", r"\bmuqarnas\b",
]
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
