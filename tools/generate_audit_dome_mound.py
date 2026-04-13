"""
generate_audit_dome_mound.py
============================
Produce supplementary/audit/dome_mound_keyword_audit.txt

Audit of the dome + hemispherical mound evolution keyword sweep (Test 2b).
Combines dome keywords (Test 2) with mound-evolution keywords from
keywords.json["mound_evolution"].

Run from repo root:
    python3 tools/generate_audit_dome_mound.py
"""

import re
import sys
import json
from pathlib import Path
from datetime import datetime, timezone

sys.path.insert(0, str(Path(__file__).parent.parent))
from data.unesco_corpus import load_corpus
from lib.beru import GERIZIM, BERU, deviation as beru_dev, tier_label
from lib.dome_filter import (
    UNAMBIGUOUS_KEYWORDS as DOME_UNAMB,
    AMBIGUOUS_KEYWORDS   as DOME_AMB,
    FORM_KEYWORD_RES,
    validate_keyword_match as dome_validate,
)

OUT = Path(__file__).parent.parent / "supplementary" / "audit" / "dome_mound_keyword_audit.txt"
SEP = "─" * 88

_KW_PATH = Path(__file__).parent.parent / "keywords.json"
with open(_KW_PATH) as f:
    _KW = json.load(f)

_evo = _KW["mound_evolution"]
MOUND_UNAMB   = _evo["mound_unambiguous"]    # tumulus, tumuli, barrow, barrows, kofun
MOUND_AMB     = _evo["mound_ambiguous"]       # mound
MOUND_POS_CTX = _evo["mound_positive_context"]
MOUND_POS_RES = [re.compile(r"\b" + re.escape(p) + r"\b", re.IGNORECASE)
                 for p in MOUND_POS_CTX]

ALL_KEYWORDS = DOME_UNAMB + DOME_AMB + MOUND_UNAMB + MOUND_AMB
KW_RES = {k: re.compile(r"\b" + re.escape(k) + r"\b", re.IGNORECASE)
          for k in ALL_KEYWORDS}

def strip_html(t):
    return re.sub(r"<[^>]+>", " ", t or "")

def full_text(site):
    return strip_html(" ".join(filter(None, [
        site.site, site.short_description,
        getattr(site, "long_description", ""),
        getattr(site, "justification", ""),
    ])))

def form_sentences(text, kw):
    pat = KW_RES.get(kw) or re.compile(r"\b" + re.escape(kw) + r"\b", re.IGNORECASE)
    return [s.strip() for s in re.split(r"(?<=[.!?])\s+", text) if pat.search(s)]

def validate(text_orig, kw):
    """Return (accepted: bool, verdict: str)."""
    if kw in DOME_UNAMB or kw in MOUND_UNAMB:
        return True, "unambiguous"
    if kw in DOME_AMB:
        ok, _, _note = dome_validate(text_orig, kw)
        return ok, "context-validated" if ok else "REJECTED"
    if kw == "mound":
        sents = form_sentences(text_orig, kw)
        ok = any(any(p.search(s) for p in MOUND_POS_RES) for s in sents)
        return ok, "mound+arch.context" if ok else "REJECTED"
    return True, "unambiguous"

def run():
    corpus = load_corpus()
    included, rejected = [], []

    for site in corpus:
        if site.category not in ("Cultural", "Mixed"):
            continue
        txt = full_text(site)
        txt_lo = txt.lower()

        hits = [k for k in ALL_KEYWORDS
                if re.search(r"\b" + re.escape(k) + r"\b", txt_lo)]
        if not hits:
            continue

        validated = {k: validate(txt, k) for k in hits}
        accepted_keys = [k for k, (ok, _) in validated.items() if ok]

        dev    = beru_dev(site.longitude)
        tier   = tier_label(dev)
        dev_km = dev * BERU * 111.0

        record = dict(
            name=site.site, lon=site.longitude,
            dev=dev, dev_km=dev_km, tier=tier,
            keys=hits, validated={k: v for k, (_, v) in validated.items()},
            sentences={k: form_sentences(txt, k)[:2] for k in hits},
        )
        if accepted_keys:
            included.append(record)
        else:
            rejected.append(record)

    tier_order = {"A++": 0, "A+": 1, "A": 2, "B": 3, "C": 4}
    included.sort(key=lambda r: (tier_order.get(r["tier"], 9), r["dev"]))
    rejected.sort(key=lambda r: r["name"])

    n_ap = sum(1 for r in included if r["tier"] in ("A++", "A+"))
    ts = datetime.now(timezone.utc).strftime("%a %b %d %H:%M:%S UTC %Y")
    lines = [
        "UNESCO DOME + HEMISPHERICAL MOUND EVOLUTION KEYWORD AUDIT",
        f"Generated  : {ts}",
        f"Script     : tools/generate_audit_dome_mound.py",
        f"Dome unamb : {sorted(DOME_UNAMB)}",
        f"Dome amb   : {sorted(DOME_AMB)}  (context-validated)",
        f"Mound unamb: {sorted(MOUND_UNAMB)}",
        f"Mound amb  : {MOUND_AMB}  (requires archaeological context)",
        "",
        f"INCLUDED   : {len(included)}",
        f"REJECTED   : {len(rejected)}",
        f"FP rate    : {len(rejected)/(len(included)+len(rejected))*100:.1f}%"
            if (len(included)+len(rejected)) else "n/a",
        f"A+ included: {n_ap}",
        "",
        "INCLUDED SITES:", SEP,
    ]
    for r in included:
        flag = " ★" if r["tier"] in ("A++", "A+") else ""
        lines.append(f"  {r['name']}{flag}")
        lines.append(f"    LON  : {r['lon']:.4f}°E  |  dev {r['dev']:.4f} beru"
                     f"  ({r['dev_km']:.1f} km)  |  Tier {r['tier']}")
        lines.append(f"    KEYS : {r['keys']}")
        lines.append(f"    VALID: {r['validated']}")
        for k in r["keys"]:
            for s in r["sentences"].get(k, []):
                lines.append(f"    [{k}] \"{s[:140]}\"")
        lines.append("")

    lines += ["", "CONTEXT-REJECTED SITES:", SEP]
    for r in rejected:
        lines.append(f"  {r['name']}")
        lines.append(f"    KEYS : {r['keys']}  →  {r['validated']}")
        for k in r["keys"]:
            for s in r["sentences"].get(k, []):
                lines.append(f"    [{k}] \"{s[:140]}\"")
        lines.append("")

    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text("\n".join(lines), encoding="utf-8")
    print(f"Written → {OUT}  ({len(included)} included, {len(rejected)} rejected)")

if __name__ == "__main__":
    run()
