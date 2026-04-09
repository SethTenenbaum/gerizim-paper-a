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
"""

import re
import json
from pathlib import Path

# ── Load keyword sets from config ────────────────────────────────────────────
_CONFIG_PATH = Path(__file__).parent.parent / "config.json"
with open(_CONFIG_PATH) as _f:
    _CONFIG = json.load(_f)

_kw = _CONFIG["keywords"]

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

_F_cfg = _kw["founding_capital"]
F_UNAMBIGUOUS = _F_cfg["unambiguous"]
F_AMBIGUOUS   = _F_cfg["ambiguous"]
F_ALL         = F_UNAMBIGUOUS + F_AMBIGUOUS

F_KEYWORD_RES = {
    kw: re.compile(r"\b" + re.escape(kw) + r"\b", re.IGNORECASE)
    if " " in kw else re.compile(re.escape(kw), re.IGNORECASE)
    for kw in F_ALL
}

# ── Positive context: sentence must mention a polity / political entity ──────
F_POSITIVE_CONTEXT = [
    # Polity types
    r"\bkingdom\b", r"\bempire\b", r"\bsultanate\b", r"\bcaliphate\b",
    r"\bdynasty\b", r"\bkhanate\b", r"\brepublic\b", r"\bstate\b",
    r"\bprovince\b", r"\bcolony\b", r"\bviceroyalty\b",
    r"\bcaptaincy\b", r"\bprefecture\b", r"\bprincipality\b",
    r"\bconfederac\w+\b", r"\bregion\b", r"\bcounty\b",
    # Political descriptors
    r"\bpolitical\b", r"\badministrat\w+\b", r"\bimperial\b",
    r"\broyal\b", r"\bceremonial\b", r"\bseat of\b",
    # City-founding language
    r"\bcity\b", r"\btown\b", r"\burban\b", r"\bsettlement\b",
    r"\bmetropol\w+\b", r"\bgarrison\b", r"\bcitadel\b",
    r"\bport\b", r"\bfortress\b", r"\bcolonial\b",
    # "capital of [proper noun]" — the phrase "capital of" itself implies a polity
    # if followed by a proper-noun-looking word (uppercase initial in original)
    # We accept any "capital of" where no negative pattern fires, since
    # the concept of "capital of X" almost always means a political capital.
    r"\bcapital of\b",
]
F_POSITIVE_RES = [re.compile(p, re.IGNORECASE) for p in F_POSITIVE_CONTEXT]

# ── Negative context: sentence describes something other than a capital ──────
# NOTE: Only domain-level exclusions here. No specific proper nouns,
# art movements, or items identified by inspecting the data post-hoc.
F_NEGATIVE_CONTEXT = [
    # Religious institution founding (→ not a political capital)
    r"\bmonaster\w*\b", r"\bchurch\b", r"\bcathedral\b",
    r"\btemple\b", r"\bmosque\b", r"\bshrine\b",
    r"\bhermit\b",
]
F_NEGATIVE_RES = [re.compile(p, re.IGNORECASE) for p in F_NEGATIVE_CONTEXT]


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

_S_cfg = _kw["sacred_origin"]
S_UNAMBIGUOUS = _S_cfg["unambiguous"]
S_AMBIGUOUS   = _S_cfg["ambiguous"]
S_ALL         = S_UNAMBIGUOUS + S_AMBIGUOUS

S_KEYWORD_RES = {
    kw: re.compile(r"\b" + re.escape(kw) + r"\b", re.IGNORECASE)
    if " " in kw else re.compile(re.escape(kw), re.IGNORECASE)
    for kw in S_ALL
}

# ── Positive context: sentence must describe a religious/sacred tradition ────
S_POSITIVE_CONTEXT = [
    # Religious traditions
    r"\breligio\w+\b", r"\bfaith\b", r"\bspiritual\b",
    r"\bbuddhis[mt]\b", r"\bbuddha\b",
    r"\bchristian\w*\b", r"\bislam\w*\b", r"\bjudai\w+\b",
    r"\bhindu\w*\b", r"\bshinto\b", r"\btaois[mt]\b",
    r"\bconfucian\w*\b", r"\bzoroastrian\w*\b",
    r"\bfranciscan\b", r"\bdominican\b", r"\bbenedictine\b",
    r"\bjesuits?\b", r"\bmonastic\b",
    # Sacred entities
    r"\bgod\b", r"\bdeity\b", r"\bsaint\b", r"\bprophet\b",
    r"\bapostle\b", r"\brelics?\b", r"\bmartyrs?\b",
    r"\bdevoti\w+\b", r"\bworship\b", r"\bprayer\b",
    r"\bsacred\b", r"\bholy\b", r"\bdivine\b",
    # Pilgrimage-related
    r"\bpilgrim\w*\b", r"\bshrine\b", r"\bsanctuar\w+\b",
    r"\btemple\b", r"\bmosque\b", r"\bchurch\b", r"\bsynagogue\b",
    r"\bcathedral\b", r"\bbasilica\b",
]
S_POSITIVE_RES = [re.compile(p, re.IGNORECASE) for p in S_POSITIVE_CONTEXT]

# ── Negative context: sentence is clearly in a non-sacred domain ─────────────
# NOTE: Only domain-level exclusions. The positive list does the real work
# (requiring religious/sacred vocabulary in the same sentence). We do NOT
# list specific secular movements, people, or art styles seen in the data —
# that would be post-hoc Texas Sharpshooter targeting.
S_NEGATIVE_CONTEXT = [
    # Natural science domain (species origins, geological evolution)
    r"\bspecies\b", r"\bbiodiversity\b", r"\bglacial\b",
    r"\bgeological\b", r"\bgeomorpholog\w+\b",
]
S_NEGATIVE_RES = [re.compile(p, re.IGNORECASE) for p in S_NEGATIVE_CONTEXT]


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

_M_cfg = _kw["founding_monument"]
M_UNAMBIGUOUS = _M_cfg["unambiguous"]
M_AMBIGUOUS   = _M_cfg["ambiguous"]
M_ALL         = M_UNAMBIGUOUS + M_AMBIGUOUS

M_KEYWORD_RES = {
    kw: re.compile(r"\b" + re.escape(kw) + r"\b", re.IGNORECASE)
    if " " in kw else re.compile(re.escape(kw), re.IGNORECASE)
    for kw in M_ALL
}

# ── Positive context: must describe a built form that originated a tradition ─
M_POSITIVE_CONTEXT = [
    # Monumental/architectural terms
    r"\barchitect\w*\b", r"\bmonument\w*\b", r"\bbuilding\w*\b",
    r"\btemple\w*\b", r"\bchurch\w*\b", r"\bmosque\w*\b",
    r"\bpalace\w*\b", r"\btomb\w*\b", r"\bmausoleum\w*\b",
    r"\bfortress\w*\b", r"\bcathedral\w*\b", r"\bbasilica\w*\b",
    r"\bstupa\w*\b", r"\bpagoda\w*\b", r"\btower\w*\b",
    r"\bcastle\w*\b", r"\bwall\w*\b", r"\bbridge\w*\b",
    r"\bdam\w?\b", r"\bcanal\w*\b", r"\baqueduct\w*\b",
    r"\btype\b", r"\bstyle\b", r"\btradition\b", r"\btypolog\w+\b",
    r"\bconstruct\w*\b", r"\bbuilt\b", r"\berected\b",
    # Innovation/influence language
    r"\binspired\b", r"\binfluence\w*\b", r"\bprecursor\b",
    r"\bpioneer\w*\b", r"\boriginated\b", r"\bprototyp\w*\b",
    r"\binnovati\w+\b", r"\brevolution\w*\b",
]
M_POSITIVE_RES = [re.compile(p, re.IGNORECASE) for p in M_POSITIVE_CONTEXT]

# ── Negative context: generic UNESCO praise, not prototype language ──────────
M_NEGATIVE_CONTEXT = [
    # Generic superlatives applied to any site
    r"\bmasterpiece of .{0,20}(human creative genius|creative genius)",
    r"\bmasterpiece of .{0,20}landscape\b",
    r"\boutstanding example of\b",
]
M_NEGATIVE_RES = [re.compile(p, re.IGNORECASE) for p in M_NEGATIVE_CONTEXT]


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

_X_cfg = _kw["founding_axis"]
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

# ── Positive context: must be a named/major route or fortification system ────
X_POSITIVE_CONTEXT = [
    # Named routes
    r"\bsilk road\b", r"\bincense route\b", r"\bspice route\b",
    r"\bappian\b", r"\bvia appia\b", r"\broyal road\b",
    r"\bqhapaq\b", r"\binca\b",
    # Major infrastructure
    r"\bempire\b", r"\bimperial\b", r"\bfrontier\b", r"\blimes\b",
    r"\btranscontinental\b", r"\bnetwork\b",
    r"\bstrategic\b", r"\bmilitary\b", r"\bconquest\b",
    # Scale markers
    r"\bthousands? of (km|kilometres|kilometers|miles)\b",
    r"\b\d{3,}\s*(km|kilometres|kilometers|miles)\b",
    r"\bmore than \d{2,}\s*(km|kilometres|kilometers)\b",
    # Linking/connecting language
    r"\blinking\b", r"\bconnect\w+\b",
]
X_POSITIVE_RES = [re.compile(p, re.IGNORECASE) for p in X_POSITIVE_CONTEXT]

# ── Negative context: just a city with walls, not a founding axis ────────────
X_NEGATIVE_CONTEXT = [
    # Defensive walls of a single city (very common in UNESCO descriptions)
    r"\bcity wall\w*\b", r"\btown wall\w*\b",
    r"\bfortified (city|town|settlement)\b",
    r"\benclosed\b", r"\bencircl\w+\b",
    r"\brampart\w*\b",
]
X_NEGATIVE_RES = [re.compile(p, re.IGNORECASE) for p in X_NEGATIVE_CONTEXT]


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

_L_cfg = _kw["ancient_landscape"]
L_UNAMBIGUOUS = _L_cfg["unambiguous"]
L_AMBIGUOUS   = _L_cfg["ambiguous"]
L_ALL         = L_UNAMBIGUOUS + L_AMBIGUOUS

L_KEYWORD_RES = {
    kw: re.compile(r"\b" + re.escape(kw) + r"\b", re.IGNORECASE)
    if " " in kw else re.compile(re.escape(kw), re.IGNORECASE)
    for kw in L_ALL
}

# ── Positive context: must describe antiquity or long-term continuity ────────
L_POSITIVE_CONTEXT = [
    # Time-depth markers
    r"\bmillenn\w+\b", r"\bcenturies\b",
    r"\bthousands? of years\b",
    r"\bancient\b", r"\bprehistoric\b", r"\bneolithic\b",
    r"\bbronze age\b", r"\biron age\b",
    r"\boldest\b", r"\blong.?standing\b",
    r"\bcontinuous\w*\b", r"\buninterrupt\w+\b",
    r"\btraditional .{0,20}(agricultur|land use|practic|farming)",
    r"\bterracing\b", r"\birrigation\b",
    r"\bagricultural\b", r"\bpastoral\b", r"\bagrofore\w+\b",
    # Living landscape language
    r"\bliving landscape\b", r"\brelict landscape\b",
    r"\beveryday life\b", r"\blived.in\b",
    r"\bsacred .{0,10}landscape\b",
]
L_POSITIVE_RES = [re.compile(p, re.IGNORECASE) for p in L_POSITIVE_CONTEXT]

# ── Negative context: modern or purely urban "cultural landscape" ────────────
L_NEGATIVE_CONTEXT = [
    # Modern urban / designed landscapes (not "ancient")
    r"\b(18th|19th|20th|21st).century\b",
    r"\bmodern\b", r"\bcontemporary\b",
    r"\burban (planning|design|develop)\w*\b",
    r"\bparks? and (garden|boulevard)\w*\b",
    r"\bbotanical garden\b",
    r"\b(palace|royal) garden\b",
    r"\bdesigned landscape\b",
]
L_NEGATIVE_RES = [re.compile(p, re.IGNORECASE) for p in L_NEGATIVE_CONTEXT]


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
