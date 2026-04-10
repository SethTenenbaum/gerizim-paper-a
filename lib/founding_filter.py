"""
founding_filter.py
==================
Context-aware keyword matching for founding/sacred site classification.

Modelled on dome_filter.py: every category has unambiguous keywords
(accepted on sight) and ambiguous keywords (require sentence-level
context validation to reduce false positives).

CATEGORIES
──────────
  F  Founding Capital / Seat of Power
  S  Sacred Origin / Birth of a Tradition
  M  Founding Monument / Prototype
  X  Founding Axis / Imperial Infrastructure
  L  Ancient Continuous Landscape

USAGE
─────
  from lib.founding_filter import classify_site, CATEGORY_LABELS

  cats = classify_site(site_obj)   # → {'F', 'S'}  or set()

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

CATEGORY_LABELS = {
    "F": "FOUNDING CAPITAL / SEAT OF POWER",
    "S": "SACRED ORIGIN / BIRTH OF A TRADITION",
    "M": "FOUNDING MONUMENT / PROTOTYPE",
    "X": "FOUNDING AXIS / IMPERIAL INFRASTRUCTURE",
    "L": "ANCIENT CONTINUOUS LANDSCAPE",
    "?": "UNCLASSIFIED (no matching keyword)",
}

PRIORITY = ["F", "S", "M", "X", "L", "?"]


# ══════════════════════════════════════════════════════════════════════════════
#  F — FOUNDING CAPITAL
# ══════════════════════════════════════════════════════════════════════════════

_F_cfg = _KW["founding_capital"]
F_UNAMBIGUOUS = _F_cfg["unambiguous"]
F_AMBIGUOUS   = _F_cfg["ambiguous"]
F_ALL         = F_UNAMBIGUOUS + F_AMBIGUOUS

F_KEYWORD_RES = {
    kw: re.compile(r"\b" + re.escape(kw) + r"\b", re.IGNORECASE)
    if " " in kw else re.compile(re.escape(kw), re.IGNORECASE)
    for kw in F_ALL
}

F_POSITIVE_CONTEXT = _F_cfg["positive_context"]
F_POSITIVE_RES     = [re.compile(p, re.IGNORECASE) for p in F_POSITIVE_CONTEXT]

F_NEGATIVE_CONTEXT = _F_cfg["negative_context"]
F_NEGATIVE_RES     = [re.compile(p, re.IGNORECASE) for p in F_NEGATIVE_CONTEXT]


def _validate_F_sentence(sentence: str) -> bool:
    """Validate that an ambiguous F keyword appears in a political-capital context."""
    for pat in F_NEGATIVE_RES:
        if pat.search(sentence):
            return False
    for pat in F_POSITIVE_RES:
        if pat.search(sentence):
            return True
    return False


# ══════════════════════════════════════════════════════════════════════════════
#  S — SACRED ORIGIN
# ══════════════════════════════════════════════════════════════════════════════

_S_cfg = _KW["sacred_origin"]
S_UNAMBIGUOUS = _S_cfg["unambiguous"]
S_AMBIGUOUS   = _S_cfg["ambiguous"]
S_ALL         = S_UNAMBIGUOUS + S_AMBIGUOUS

S_KEYWORD_RES = {
    kw: re.compile(r"\b" + re.escape(kw) + r"\b", re.IGNORECASE)
    if " " in kw else re.compile(re.escape(kw), re.IGNORECASE)
    for kw in S_ALL
}

S_POSITIVE_CONTEXT = _S_cfg["positive_context"]
S_POSITIVE_RES     = [re.compile(p, re.IGNORECASE) for p in S_POSITIVE_CONTEXT]

S_NEGATIVE_CONTEXT = _S_cfg["negative_context"]
S_NEGATIVE_RES     = [re.compile(p, re.IGNORECASE) for p in S_NEGATIVE_CONTEXT]


def _validate_S_sentence(sentence: str) -> bool:
    """Validate that an ambiguous S keyword appears in a religious/sacred context."""
    for pat in S_NEGATIVE_RES:
        if pat.search(sentence):
            return False
    for pat in S_POSITIVE_RES:
        if pat.search(sentence):
            return True
    return False


# ══════════════════════════════════════════════════════════════════════════════
#  M — FOUNDING MONUMENT
# ══════════════════════════════════════════════════════════════════════════════

_M_cfg = _KW["founding_monument"]
M_UNAMBIGUOUS = _M_cfg["unambiguous"]
M_AMBIGUOUS   = _M_cfg["ambiguous"]
M_ALL         = M_UNAMBIGUOUS + M_AMBIGUOUS

M_KEYWORD_RES = {
    kw: re.compile(r"\b" + re.escape(kw) + r"\b", re.IGNORECASE)
    if " " in kw else re.compile(re.escape(kw), re.IGNORECASE)
    for kw in M_ALL
}

M_POSITIVE_CONTEXT = _M_cfg["positive_context"]
M_POSITIVE_RES     = [re.compile(p, re.IGNORECASE) for p in M_POSITIVE_CONTEXT]

M_NEGATIVE_CONTEXT = _M_cfg["negative_context"]
M_NEGATIVE_RES     = [re.compile(p, re.IGNORECASE) for p in M_NEGATIVE_CONTEXT]


def _validate_M_sentence(sentence: str) -> bool:
    """Validate that an ambiguous M keyword describes a founding/prototype monument."""
    for pat in M_NEGATIVE_RES:
        if pat.search(sentence):
            return False
    for pat in M_POSITIVE_RES:
        if pat.search(sentence):
            return True
    return False


# ══════════════════════════════════════════════════════════════════════════════
#  X — FOUNDING AXIS
# ══════════════════════════════════════════════════════════════════════════════

_X_cfg = _KW["founding_axis"]
X_UNAMBIGUOUS = _X_cfg["unambiguous"]
X_AMBIGUOUS   = _X_cfg["ambiguous"]
X_ALL         = X_UNAMBIGUOUS + X_AMBIGUOUS

# Regex-based unambiguous patterns (match against full text, no context needed)
X_REGEX_UNAMBIGUOUS = [
    re.compile(p, re.IGNORECASE)
    for p in _X_cfg.get("regex_unambiguous", [])
]

X_KEYWORD_RES = {
    kw: re.compile(r"\b" + re.escape(kw) + r"\b", re.IGNORECASE)
    if " " in kw else re.compile(re.escape(kw), re.IGNORECASE)
    for kw in X_ALL
}

X_POSITIVE_CONTEXT = _X_cfg["positive_context"]
X_POSITIVE_RES     = [re.compile(p, re.IGNORECASE) for p in X_POSITIVE_CONTEXT]

X_NEGATIVE_CONTEXT = _X_cfg["negative_context"]
X_NEGATIVE_RES     = [re.compile(p, re.IGNORECASE) for p in X_NEGATIVE_CONTEXT]


def _validate_X_sentence(sentence: str) -> bool:
    """Validate that an ambiguous X keyword describes a major route/fortification."""
    for pat in X_NEGATIVE_RES:
        if pat.search(sentence):
            return False
    for pat in X_POSITIVE_RES:
        if pat.search(sentence):
            return True
    return False


# ══════════════════════════════════════════════════════════════════════════════
#  L — ANCIENT CONTINUOUS LANDSCAPE
# ══════════════════════════════════════════════════════════════════════════════

_L_cfg = _KW["ancient_landscape"]
L_UNAMBIGUOUS = _L_cfg["unambiguous"]
L_AMBIGUOUS   = _L_cfg["ambiguous"]
L_ALL         = L_UNAMBIGUOUS + L_AMBIGUOUS

L_KEYWORD_RES = {
    kw: re.compile(r"\b" + re.escape(kw) + r"\b", re.IGNORECASE)
    if " " in kw else re.compile(re.escape(kw), re.IGNORECASE)
    for kw in L_ALL
}

L_POSITIVE_CONTEXT = _L_cfg["positive_context"]
L_POSITIVE_RES     = [re.compile(p, re.IGNORECASE) for p in L_POSITIVE_CONTEXT]

L_NEGATIVE_CONTEXT = _L_cfg["negative_context"]
L_NEGATIVE_RES     = [re.compile(p, re.IGNORECASE) for p in L_NEGATIVE_CONTEXT]


def _validate_L_sentence(sentence: str) -> bool:
    """Validate that an ambiguous L keyword describes an ancient continuous landscape."""
    for pat in L_NEGATIVE_RES:
        if pat.search(sentence):
            return False
    for pat in L_POSITIVE_RES:
        if pat.search(sentence):
            return True
    return False


# ══════════════════════════════════════════════════════════════════════════════
#  GENERIC VALIDATION DISPATCHER
# ══════════════════════════════════════════════════════════════════════════════

_VALIDATORS = {
    "F": (_validate_F_sentence, F_UNAMBIGUOUS, F_AMBIGUOUS, F_KEYWORD_RES),
    "S": (_validate_S_sentence, S_UNAMBIGUOUS, S_AMBIGUOUS, S_KEYWORD_RES),
    "M": (_validate_M_sentence, M_UNAMBIGUOUS, M_AMBIGUOUS, M_KEYWORD_RES),
    "X": (_validate_X_sentence, X_UNAMBIGUOUS, X_AMBIGUOUS, X_KEYWORD_RES),
    "L": (_validate_L_sentence, L_UNAMBIGUOUS, L_AMBIGUOUS, L_KEYWORD_RES),
}


def _check_category(cat: str, full_text: str, sentences: list[str]) -> tuple[bool, list[str]]:
    """
    Check if a text qualifies for category `cat`.

    Returns (matched: bool, reasons: list[str]).
    """
    validate_fn, unambiguous, ambiguous, kw_res = _VALIDATORS[cat]
    reasons = []

    # Stage 1: Check for unambiguous keywords anywhere in full_text
    for kw in unambiguous:
        if kw_res[kw].search(full_text):
            reasons.append(f"unambiguous:{kw}")
            return True, reasons

    # Stage 1b: Check regex-based unambiguous patterns (X category)
    if cat == "X":
        for pat in X_REGEX_UNAMBIGUOUS:
            m = pat.search(full_text)
            if m:
                reasons.append(f"regex-unambiguous:{pat.pattern}")
                return True, reasons

    # Stage 2: Check ambiguous keywords with sentence-level context validation
    for kw in ambiguous:
        kw_re = kw_res[kw]
        for sent in sentences:
            if kw_re.search(sent) and validate_fn(sent):
                reasons.append(f"context-validated:{kw}")
                return True, reasons

    return False, []


def classify_site(site_obj) -> dict[str, list[str]]:
    """
    Classify a UNESCO site into founding/sacred categories.

    Parameters
    ----------
    site_obj : UNESCOSite
        Must have .full_text (site name + short_description + justification + extended_description).

    Returns
    -------
    dict
        Maps category code → list of reason strings.
        Empty dict if no category matches.
        Example: {'F': ['unambiguous:first capital'], 'S': ['context-validated:birthplace of']}
    """
    full_text = site_obj.full_text  # already lowercased (includes all available text + extended data)

    # Split full_text into sentences for context validation
    # (full_text includes site name, short_description, justification, extended_description)
    sentences = re.split(r"(?<=[.!?])\s+", full_text)

    result = {}
    for cat in PRIORITY[:-1]:  # skip "?"
        matched, reasons = _check_category(cat, full_text, sentences)
        if matched:
            result[cat] = reasons

    return result


def classify_text(text: str) -> dict[str, list[str]]:
    """
    Classify a raw text string (for use without a UNESCOSite object).

    Parameters
    ----------
    text : str
        The full text to classify (will be lowercased internally).

    Returns
    -------
    dict
        Same format as classify_site().
    """
    text_lower = text.lower()
    sentences = re.split(r"(?<=[.!?])\s+", text_lower)

    result = {}
    for cat in PRIORITY[:-1]:
        matched, reasons = _check_category(cat, text_lower, sentences)
        if matched:
            result[cat] = reasons

    return result


def primary_category(cats: dict) -> str:
    """Return the highest-priority category from a classify_site() result."""
    for c in PRIORITY:
        if c in cats:
            return c
    return "?"
