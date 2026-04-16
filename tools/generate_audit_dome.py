"""
generate_audit_dome.py
======================
Produce supplementary/audit/dome_keyword_audit.txt

Audit of the domed/spherical monument keyword sweep (Test 2).
Keywords: stupa, stupas, tholos (unambiguous) +
          dome, domed, domes, spherical (context-validated).

Run from repo root:
    python3 tools/generate_audit_dome.py
"""

import re
import sys
from pathlib import Path
from datetime import datetime, timezone

sys.path.insert(0, str(Path(__file__).parent.parent))
from data.unesco_corpus import load_corpus
from lib.beru import GERIZIM, BERU, TIER_APLUS, deviation as beru_dev, tier_label
from lib.dome_filter import (
    UNAMBIGUOUS_KEYWORDS, AMBIGUOUS_KEYWORDS, FORM_KEYWORD_RES,
    validate_keyword_match,
    rejection_reason,
)

OUT = Path(__file__).parent.parent / "supplementary" / "audit" / "dome_keyword_audit.txt"

SEP = "─" * 88

def strip_html(t):
    return re.sub(r"<[^>]+>", " ", t or "")

def form_sentences(text, keywords):
    sents = re.split(r"(?<=[.!?])\s+", text)
    return [s.strip() for s in sents
            if any(FORM_KEYWORD_RES[k].search(s) for k in keywords if k in FORM_KEYWORD_RES)]

def run():
    corpus = load_corpus()
    included = []
    rejected = []

    for site in corpus:
        if site.category not in ("Cultural", "Mixed"):
            continue
        if site.longitude is None:
            continue

        # Use full_text (lowercased) for hit detection — includes extended_description
        full = site.full_text  # already lowercased

        hit_unamb = [k for k in UNAMBIGUOUS_KEYWORDS
                     if re.search(r"\b" + re.escape(k) + r"\b", full)]
        hit_amb   = [k for k in AMBIGUOUS_KEYWORDS
                     if re.search(r"\b" + re.escape(k) + r"\b", full)]

        if not hit_unamb and not hit_amb:
            continue

        # Preserve original casing for sentence extraction / context validation
        full_orig = strip_html(" ".join(filter(None, [
            site.site, site.short_description,
            site.justification,
            getattr(site, "extended_description", ""),
        ])))

        validated = {}
        reject_reasons = {}   # keyword → reason string (for rejected keywords)
        for k in hit_unamb:
            validated[k] = "unambiguous"
        for k in hit_amb:
            ok, matched_sents, _note = validate_keyword_match(full_orig, k)
            validated[k] = "context-validated" if ok else "REJECTED"
            if not ok:
                # Find the first matching sentence and explain why it failed
                kw_re = FORM_KEYWORD_RES[k]
                sents = re.split(r"(?<=[.!?])\s+", full_orig)
                hit_sents = [s.strip() for s in sents if kw_re.search(s)]
                if hit_sents:
                    reject_reasons[k] = rejection_reason(hit_sents[0])
                else:
                    reject_reasons[k] = "no matching sentence found"

        all_keys = hit_unamb + hit_amb
        accepted_keys = [k for k in all_keys if validated.get(k) != "REJECTED"]
        is_rejected = len(accepted_keys) == 0

        dev = beru_dev(site.longitude)
        tier = tier_label(dev)
        dev_km = dev * BERU * 111.0

        record = dict(
            name=site.site,
            lon=site.longitude,
            dev=dev,
            dev_km=dev_km,
            tier=tier,
            keys=all_keys,
            validated=validated,
            reject_reasons=reject_reasons,
            sentences={k: form_sentences(full_orig, [k]) for k in all_keys},
        )
        if is_rejected:
            rejected.append(record)
        else:
            included.append(record)

    # sort included by tier then dev
    tier_order = {"A++": 0, "A+": 1, "A": 2, "B": 3, "C": 4, "C-": 5, "C--": 6}
    included.sort(key=lambda r: (tier_order.get(r["tier"], 9), r["dev"]))
    rejected.sort(key=lambda r: r["name"])

    lines = []
    ts = datetime.now(timezone.utc).strftime("%a %b %d %H:%M:%S UTC %Y")
    lines += [
        "UNESCO DOME / SPHERICAL MONUMENT KEYWORD AUDIT",
        f"Generated : {ts}",
        f"Script    : tools/generate_audit_dome.py",
        f"Keywords  : {sorted(UNAMBIGUOUS_KEYWORDS + AMBIGUOUS_KEYWORDS)}",
        f"Unambiguous (always included): {sorted(UNAMBIGUOUS_KEYWORDS)}",
        f"Ambiguous (context-validated): {sorted(AMBIGUOUS_KEYWORDS)}",
        "",
        f"INCLUDED  : {len(included)}",
        f"REJECTED  : {len(rejected)}",
        f"FP rate   : {len(rejected)/(len(included)+len(rejected))*100:.1f}%"
            if (len(included)+len(rejected)) > 0 else "FP rate: n/a",
        "",
        f"A+ among included : {sum(1 for r in included if r['tier'] in ('A++','A+'))}",
        "",
    ]

    lines += ["INCLUDED SITES:", SEP]
    for r in included:
        tier_flag = " ★" if r["tier"] in ("A++", "A+") else ""
        lines.append(f"  {r['name']}{tier_flag}")
        lines.append(f"    LON  : {r['lon']:.4f}°E  |  dev {r['dev']:.4f} beru"
                     f"  ({r['dev_km']:.1f} km)  |  Tier {r['tier']}")
        lines.append(f"    KEYS : {r['keys']}  →  {r['validated']}")
        for k in r["keys"]:
            for s in r["sentences"].get(k, [])[:2]:
                lines.append(f"    [{k}] \"{s[:140]}\"")
        lines.append("")

    lines += ["", "CONTEXT-REJECTED SITES:", SEP]
    for r in rejected:
        lines.append(f"  {r['name']}")
        lines.append(f"    KEYS : {r['keys']}  →  {r['validated']}")
        for k in r["keys"]:
            reason = r["reject_reasons"].get(k)
            if reason:
                lines.append(f"    REASON [{k}]: {reason}")
            for s in r["sentences"].get(k, [])[:2]:
                lines.append(f"    [{k}] \"{s[:140]}\"")
        lines.append("")

    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text("\n".join(lines), encoding="utf-8")
    print(f"Written → {OUT}  ({len(included)} included, {len(rejected)} rejected)")

if __name__ == "__main__":
    run()
