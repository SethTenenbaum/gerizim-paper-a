"""
generate_audit_founding.py
==========================
Produce supplementary/audit/founding_keyword_audit.txt

Audit of the founding/sacred-origin site keyword classifier (Test 3 / Table 2).
Uses lib/founding_filter.py, which loads all keyword lists from keywords.json.

Five categories are reported:
  F  Founding Capital / Seat of Power
  S  Sacred Origin / Birth of a Tradition
  M  Founding Monument / Prototype
  X  Founding Axis / Imperial Infrastructure
  L  Ancient Continuous Landscape

For each UNESCO Cultural/Mixed site the script shows:
  • Whether the site was included (≥1 accepted category)
  • Which categories matched and why (unambiguous / context-validated / REJECTED)
  • The sentence(s) that triggered the match
  • The site's beru deviation and tier

Sites are sorted by tier (A++ → A+ → A → B → C) within each output section.

Run from repo root:
    python3 tools/generate_audit_founding.py
"""

import re
import sys
from pathlib import Path
from datetime import datetime, timezone

sys.path.insert(0, str(Path(__file__).parent.parent))
from data.unesco_corpus import load_corpus
from lib.beru import GERIZIM, BERU, TIER_APLUS, deviation as beru_dev, tier_label
from lib.founding_filter import (
    classify_site,
    CATEGORY_LABELS,
    PRIORITY,
    primary_category,
    # Per-category data for the keyword summary header
    F_UNAMBIGUOUS, F_AMBIGUOUS,
    S_UNAMBIGUOUS, S_AMBIGUOUS,
    M_UNAMBIGUOUS, M_AMBIGUOUS,
    X_UNAMBIGUOUS, X_AMBIGUOUS,
    L_UNAMBIGUOUS, L_AMBIGUOUS,
)

OUT = Path(__file__).parent.parent / "supplementary" / "audit" / "founding_keyword_audit.txt"
SEP = "─" * 88

TIER_ORDER = {"A++": 0, "A+": 1, "A": 2, "B": 3, "C": 4}


def sort_key(r):
    return (TIER_ORDER.get(r["tier"], 9), r["dev"])


def get_match_sentences(site, cat_reasons):
    """
    Extract the 1-2 trigger sentences for each matched keyword/reason from the
    site's full_text (already lowercased inside the property).
    We re-use site.full_text for sentence splitting.
    """
    full = site.full_text  # lowercased
    sentences = re.split(r"(?<=[.!?])\s+", full)
    result = {}
    for reason in cat_reasons:
        # reason is like "unambiguous:first capital" or "context-validated:birthplace of"
        # or "regex-unambiguous:<pattern>"
        if ":" in reason:
            kind, kw = reason.split(":", 1)
        else:
            kw = reason
        if kind == "regex-unambiguous":
            # kw is a regex pattern
            try:
                pat = re.compile(kw, re.IGNORECASE)
                hits = [s.strip() for s in sentences if pat.search(s)]
            except re.error:
                hits = []
        else:
            # Plain keyword — word-boundary match
            try:
                pat = re.compile(r"\b" + re.escape(kw) + r"\b", re.IGNORECASE)
                hits = [s.strip() for s in sentences if pat.search(s)]
            except re.error:
                hits = []
        if hits:
            result[reason] = hits[:2]
    return result


def run():
    corpus = load_corpus()
    included = []   # sites with ≥1 accepted category
    rejected = []   # sites that hit a keyword but every keyword was rejected

    # Collect sites that mention ANY keyword from any category
    # (founding_filter handles context-validation internally)
    # We also want to capture "near-misses" = sites that have keyword hits but
    # all were context-rejected.  For those we fall back to classify_text on
    # the XML-only text to confirm rejection.

    for site in corpus:
        if site.category not in ("Cultural", "Mixed"):
            continue
        if site.longitude is None:
            continue

        cats = classify_site(site)  # dict cat → [reasons] or {}

        dev    = beru_dev(site.longitude)
        tier   = tier_label(dev)
        dev_km = dev * BERU * 111.0

        primary = primary_category(cats)

        record = dict(
            name   = site.site,
            lon    = site.longitude,
            dev    = dev,
            dev_km = dev_km,
            tier   = tier,
            cats   = cats,
            primary= primary,
        )

        if cats:
            # Attach trigger sentences for each matched category
            record["sentences"] = {
                cat: get_match_sentences(site, reasons)
                for cat, reasons in cats.items()
            }
            included.append(record)
        # (sites with no keyword hit at all are not logged — they're the
        # vast majority and would balloon the file)

    included.sort(key=sort_key)

    # ── A+ subset stats ────────────────────────────────────────────────────────
    n_total    = sum(1 for s in corpus if s.category in ("Cultural", "Mixed"))
    n_included = len(included)
    n_ap_incl  = sum(1 for r in included if r["tier"] in ("A++", "A+"))

    # Category breakdown across included sites
    cat_counts = {c: 0 for c in "FSMLX"}
    for r in included:
        for c in r["cats"]:
            cat_counts[c] = cat_counts.get(c, 0) + 1

    # ── Compose output ─────────────────────────────────────────────────────────
    ts = datetime.now(timezone.utc).strftime("%a %b %d %H:%M:%S UTC %Y")
    lines = [
        "UNESCO FOUNDING / SACRED-ORIGIN SITE KEYWORD AUDIT",
        f"Generated  : {ts}",
        f"Script     : tools/generate_audit_founding.py",
        f"Keywords   : loaded from keywords.json via lib/founding_filter.py",
        "",
        "KEYWORD LISTS (from keywords.json)",
        SEP,
        f"  F (Founding Capital)    unamb : {sorted(F_UNAMBIGUOUS)}",
        f"  F                       amb   : {sorted(F_AMBIGUOUS)}",
        f"  S (Sacred Origin)       unamb : {sorted(S_UNAMBIGUOUS)}",
        f"  S                       amb   : {sorted(S_AMBIGUOUS)}",
        f"  M (Founding Monument)   unamb : {sorted(M_UNAMBIGUOUS)}",
        f"  M                       amb   : {sorted(M_AMBIGUOUS)}",
        f"  X (Founding Axis)       unamb : {sorted(X_UNAMBIGUOUS)}",
        f"  X                       amb   : {sorted(X_AMBIGUOUS)}",
        f"  L (Ancient Landscape)   unamb : {sorted(L_UNAMBIGUOUS)}",
        f"  L                       amb   : {sorted(L_AMBIGUOUS)}",
        "",
        "SUMMARY",
        SEP,
        f"  Cultural/Mixed sites in corpus : {n_total}",
        f"  Classified (≥1 category)       : {n_included}",
        f"  A+ among classified            : {n_ap_incl}",
        "",
        "  Category breakdown (sites may appear in multiple categories):",
    ]
    for c in "FSMLX":
        lines.append(f"    {c} — {CATEGORY_LABELS[c]}: {cat_counts.get(c, 0)}")
    lines.append("")

    # ── Included sites ─────────────────────────────────────────────────────────
    lines += ["INCLUDED SITES (sorted by tier, then beru deviation):", SEP]
    for r in included:
        flag = " ★" if r["tier"] in ("A++", "A+") else ""
        pcat = r["primary"]
        lines.append(f"  {r['name']}{flag}")
        lines.append(
            f"    LON     : {r['lon']:.4f}°E  |  dev {r['dev']:.4f} beru"
            f"  ({r['dev_km']:.1f} km)  |  Tier {r['tier']}"
        )
        lines.append(f"    PRIMARY : {pcat} — {CATEGORY_LABELS.get(pcat, pcat)}")
        lines.append(f"    ALL CATS: {sorted(r['cats'].keys())}")
        for cat, reasons in sorted(r["cats"].items()):
            lines.append(f"    [{cat}] reasons : {reasons}")
            sents_by_reason = r.get("sentences", {}).get(cat, {})
            for reason, sents in sents_by_reason.items():
                for s in sents:
                    lines.append(f"      ↳ [{reason}] \"{s[:160]}\"")
        lines.append("")

    lines += [
        f"Total included: {n_included}",
        f"Generated by  : tools/generate_audit_founding.py",
        f"Timestamp     : {ts}",
    ]

    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text("\n".join(lines), encoding="utf-8")
    print(f"Written → {OUT}  ({n_included} classified sites)")


if __name__ == "__main__":
    run()
