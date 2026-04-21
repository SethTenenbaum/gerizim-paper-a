"""
anchor_site_comparison.py
=========================
Computes beru deviations for individually named sites from BOTH the
Gerizim and Jerusalem anchors.  Produces the LaTeX macro block for
GROUP 16 (Lumbini, Takht-e Soleyman, Khoja Ahmed Yasawi) and
GROUP 22 (World Peace Pagoda) in paper_a_primary_unesco.tex.

Every macro printed by this script maps 1-to-1 to a \\newcommand in
the manuscript.  If any value changes, re-run this script and update
the .tex file.

USAGE
-----
  python3 analysis/global/anchor_site_comparison.py

OUTPUT
------
  Prints LaTeX macro definitions to stdout.
  GROUP 16: Anchor-comparison macros (Lumbini, Takht, Khoja, NappTotal)
  GROUP 22: World Peace Pagoda macros (WPPDeltaKm, WPPDeltaBeru)
"""

from __future__ import annotations

import sys
from pathlib import Path
import json as _json
import numpy as np

# ── Repo root ─────────────────────────────────────────────────────────────────
ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from data.unesco_corpus import load_corpus, cultural_sites_with_coords
from lib.beru import GERIZIM, BERU, TIER_APP, TIER_APLUS, TIER_A_MAX, P_NULL_AP

# ── Constants from config.json ────────────────────────────────────────────────
_CFG      = _json.loads((ROOT / "config.json").read_text())
JERUSALEM = _CFG["anchors"]["jerusalem"]["longitude"]    # 35.2317 — from config
KM_PER_DEG = _CFG["units"]["km_per_degree"]             # 111.0

# World Peace Pagoda, Lumbini — Wikidata Q6540965 coordinate (not in UNESCO XML)
# The UNESCO XML records the Lumbini WHS centroid at 83.27611°E (used for
# \LumbiniDevKm etc above).  The Pagoda itself has a distinct Wikidata entry
# at 83.2743°E / 27.4833°N.  These are two different points ~0.9 km apart.
# Beru deviation is computed with the equatorial convention (no cos-lat),
# matching every other site in the study.
WPP_LON = 83.2743            # Wikidata Q6540965 — Nipponzan-Myohoji pagoda


def beru_dev(lon: float, anchor: float) -> tuple[float, float, str]:
    """Return (dev_beru, dev_km, tier_label) for a site relative to anchor."""
    arc  = abs(lon - anchor)
    bv   = arc / BERU
    near = round(bv * 10) / 10
    dev  = abs(bv - near)
    km   = dev * BERU * KM_PER_DEG   # equatorial, no cos-lat
    if dev <= TIER_APP:
        tier = "A++"
    elif dev <= TIER_APLUS:
        tier = "A+"
    else:
        tier = "A" if dev <= TIER_A_MAX else "B"
    return dev, km, tier


def find_site(sites: list[dict], fragment: str) -> dict | None:
    frag = fragment.lower()
    for s in sites:
        if frag in s["name"].lower():
            return s
    return None


def _html_to_latex(text: str) -> str:
    """Convert UNESCO XML HTML markup and special characters to LaTeX equivalents."""
    import re as _re
    import unicodedata as _ud
    # HTML italic/bold
    text = _re.sub(r'<[Ii]>(.*?)</[Ii]>', r'\\textit{\1}', text)
    text = _re.sub(r'<[Bb]>(.*?)</[Bb]>', r'\\textbf{\1}', text)
    # Remove any remaining HTML tags
    text = _re.sub(r'<[^>]+>', '', text)
    # Unicode → LaTeX for common diacritics (expand as needed)
    _UNICODE_TO_LATEX = {
        'š': r'\v{s}', 'Š': r'\v{S}',
        'č': r'\v{c}', 'Č': r'\v{C}',
        'ž': r'\v{z}', 'Ž': r'\v{Z}',
        'ř': r'\v{r}', 'Ř': r'\v{R}',
        'ě': r'\v{e}', 'Ě': r'\v{E}',
        'ñ': r'\~{n}', 'Ñ': r'\~{N}',
        'ó': r"\'{o}", 'Ó': r"\'{O}",
        'é': r"\'{e}", 'É': r"\'{E}",
        'á': r"\'{a}", 'Á': r"\'{A}",
        'í': r"\'{i}", 'Í': r"\'{I}",
        'ú': r"\'{u}", 'Ú': r"\'{U}",
        'ü': r'\"u',   'Ü': r'\"U',
        'ö': r'\"o',   'Ö': r'\"O',
        'ä': r'\"a',   'Ä': r'\"A',
        'ō': r'\={o}', 'Ō': r'\={O}',
        'ū': r'\={u}', 'Ū': r'\={U}',
        'ā': r'\={a}', 'Ā': r'\={A}',
    }
    for char, latex in _UNICODE_TO_LATEX.items():
        text = text.replace(char, latex)
    return text


def _fmt_lon(lon: float, dp: int = 4) -> str:
    """Format a longitude to dp decimal places, stripping trailing zeros."""
    return f"{lon:.{dp}f}".rstrip('0').rstrip('.')


def main():
    # Load UNESCO corpus
    corpus = load_corpus()
    cultural = cultural_sites_with_coords(corpus)
    sites = [{"name": s.site, "lon": s.longitude, "lat": s.latitude}
             for s in cultural]
    lons = np.array([s["lon"] for s in sites])

    # ── Count A++ from Gerizim ────────────────────────────────────────────
    arcs = np.abs(lons - GERIZIM)
    arcs = np.minimum(arcs, 360 - arcs)
    bvs  = arcs / BERU
    near = np.round(bvs * 10) / 10
    devs = np.abs(bvs - near)
    n_app = int(np.sum(devs <= TIER_APP))

    # ── All A++ sites sorted by deviation (for \topHit* and \appSite* macros) ──
    _ORDINAL_LABELS = [
        "One", "Two", "Three", "Four", "Five", "Six", "Seven",
        "Eight", "Nine", "Ten",
    ]
    site_devs = [(sites[i], devs[i]) for i in range(len(sites))]
    app_sites_sorted = sorted(
        [(s, d) for s, d in site_devs if d <= TIER_APP],
        key=lambda x: x[1]
    )
    top3 = app_sites_sorted[:3]

    # ── Jerusalem A+ count (self-excluded — corpus minus Jerusalem itself) ─
    JERUSALEM = _CFG["anchors"]["jerusalem"]["longitude"]
    jer_arcs = np.abs(lons - JERUSALEM)
    jer_arcs = np.minimum(jer_arcs, 360 - jer_arcs)
    jer_bvs  = jer_arcs / BERU
    jer_near = np.round(jer_bvs * 10) / 10
    jer_devs = np.abs(jer_bvs - jer_near)
    # exclude Jerusalem itself (the site with lon closest to JERUSALEM)
    jer_self_idx = int(np.argmin(np.abs(lons - JERUSALEM)))
    jer_devs_excl = np.delete(jer_devs, jer_self_idx)
    jer_ap_count = int(np.sum(jer_devs_excl <= TIER_APLUS))

    # ── Mecca A+ count ────────────────────────────────────────────────────
    # config.json uses display-name keys like "Mecca (39.826°E)"
    _mecca_key = next((k for k in _CFG.get("notable_anchors", {}) if "mecca" in k.lower()), None)
    if _mecca_key is None:
        MECCA_LON = 39.826   # fallback: Mecca longitude
    else:
        MECCA_LON = _CFG["notable_anchors"][_mecca_key]
    mecca_arcs = np.abs(lons - MECCA_LON)
    mecca_arcs = np.minimum(mecca_arcs, 360 - mecca_arcs)
    mecca_bvs  = mecca_arcs / BERU
    mecca_near = np.round(mecca_bvs * 10) / 10
    mecca_devs = np.abs(mecca_bvs - mecca_near)
    mecca_ap   = int(np.sum(mecca_devs <= TIER_APLUS))
    from scipy.stats import binomtest as _binom
    mecca_p = _binom(mecca_ap, len(sites), P_NULL_AP, alternative="greater").pvalue

    # ── Mt. Meru A+ count (Tanzania; 36.750°E) ────────────────────────────
    # Physical Mt. Meru, Tanzania. The cosmological/Hindu "Mt. Meru" has no
    # fixed longitude; the most cited geographic candidate is Mt. Kailash
    # (81.312°E, Tibet), tested separately below.
    MERU_LON    = 36.750
    KAILASH_LON = 81.312
    def _anchor_ap(anchor_lon):
        arcs = np.abs(lons - anchor_lon)
        arcs = np.minimum(arcs, 360 - arcs)
        bvs  = arcs / BERU
        near = np.round(bvs * 10) / 10
        devs = np.abs(bvs - near)
        ap   = int(np.sum(devs <= TIER_APLUS))
        p    = _binom(ap, len(sites), P_NULL_AP, alternative="greater").pvalue
        return ap, p
    meru_ap,    meru_p    = _anchor_ap(MERU_LON)
    kailash_ap, kailash_p = _anchor_ap(KAILASH_LON)

    print("=" * 72)
    print("  GROUP 16 — Anchor-Comparison Macros")
    print("  Script: analysis/global/anchor_site_comparison.py")
    print("=" * 72)

    # ── Lumbini ───────────────────────────────────────────────────────────
    lumbini = find_site(sites, "lumbini")
    if lumbini is None:
        sys.exit("ERROR: Lumbini not found in UNESCO XML")
    dev_g, km_g, tier_g = beru_dev(lumbini["lon"], GERIZIM)
    dev_j, km_j, tier_j = beru_dev(lumbini["lon"], JERUSALEM)
    m_g = km_g * 1000

    print(f"\n  Lumbini ({lumbini['name']})")
    print(f"    lon = {lumbini['lon']}°E")
    print(f"    vs Gerizim:   δ = {dev_g:.6f} beru = {km_g:.2f} km = {m_g:.0f} m  [{tier_g}]")
    print(f"    vs Jerusalem: δ = {dev_j:.6f} beru = {km_j:.2f} km  [{tier_j}]")
    print(f"\n  % LaTeX macros:")
    print(f"  \\newcommand{{\\LumbiniDevKm}}{{{km_g:.2f}}}          % Lumbini deviation from Gerizim (km)")
    print(f"  \\newcommand{{\\LumbiniDevM}}{{{int(round(m_g))}}}           % Lumbini deviation from Gerizim (m)")
    print(f"  \\newcommand{{\\LumbiniDevBeru}}{{{dev_g:.6f}}}      % Lumbini deviation from Gerizim (beru)")
    print(f"  \\newcommand{{\\LumbiniTier}}{{{tier_g}}}            % Lumbini tier relative to Gerizim")
    print(f"  \\newcommand{{\\LumbiniJerDevKm}}{{{km_j:.2f}}}        % Lumbini deviation from Jerusalem (km)")
    print(f"  \\newcommand{{\\LumbiniJerDevBeru}}{{{dev_j:.6f}}}    % Lumbini deviation from Jerusalem (beru)")
    print(f"  \\newcommand{{\\LumbiniJerTier}}{{{tier_j}}}          % Lumbini tier relative to Jerusalem")

    # ── Takht-e Soleyman ─────────────────────────────────────────────────
    takht = find_site(sites, "takht")
    if takht is None:
        sys.exit("ERROR: Takht-e Soleyman not found in UNESCO XML")
    dev_g, km_g, tier_g = beru_dev(takht["lon"], GERIZIM)
    dev_j, km_j, tier_j = beru_dev(takht["lon"], JERUSALEM)
    m_j = km_j * 1000

    print(f"\n  Takht-e Soleyman ({takht['name']})")
    print(f"    lon = {takht['lon']}°E")
    print(f"    vs Gerizim:   δ = {dev_g:.6f} beru = {km_g:.2f} km  [{tier_g}]")
    print(f"    vs Jerusalem: δ = {dev_j:.6f} beru = {km_j:.2f} km = {m_j:.0f} m  [{tier_j}]")
    print(f"\n  % LaTeX macros:")
    print(f"  \\newcommand{{\\TakhtGerizimDevKm}}{{{km_g:.2f}}}       % Takht-e Soleyman deviation from Gerizim (km)")
    print(f"  \\newcommand{{\\TakhtGerizimDevBeru}}{{{dev_g:.6f}}}   % Takht-e Soleyman deviation from Gerizim (beru)")
    print(f"  \\newcommand{{\\TakhtGerizimTier}}{{{tier_g}}}         % Takht-e Soleyman tier relative to Gerizim")
    print(f"  \\newcommand{{\\TakhtJerDevKm}}{{{km_j:.2f}}}          % Takht-e Soleyman deviation from Jerusalem (km)")
    print(f"  \\newcommand{{\\TakhtJerDevM}}{{{int(round(m_j))}}}         % Takht-e Soleyman deviation from Jerusalem (m)")
    print(f"  \\newcommand{{\\TakhtJerDevBeru}}{{{dev_j:.6f}}}      % Takht-e Soleyman deviation from Jerusalem (beru)")
    print(f"  \\newcommand{{\\TakhtJerTier}}{{{tier_j}}}            % Takht-e Soleyman tier relative to Jerusalem")

    # ── Khoja Ahmed Yasawi ───────────────────────────────────────────────
    khoja = find_site(sites, "khoja")
    if khoja is None:
        khoja = find_site(sites, "yasawi")
    if khoja is None:
        sys.exit("ERROR: Khoja Ahmed Yasawi not found in UNESCO XML")
    dev_g, km_g, tier_g = beru_dev(khoja["lon"], GERIZIM)
    dev_j, km_j, tier_j = beru_dev(khoja["lon"], JERUSALEM)

    print(f"\n  Khoja Ahmed Yasawi ({khoja['name']})")
    print(f"    lon = {khoja['lon']}°E")
    print(f"    vs Gerizim:   δ = {dev_g:.6f} beru = {km_g:.2f} km  [{tier_g}]")
    print(f"    vs Jerusalem: δ = {dev_j:.6f} beru = {km_j:.2f} km  [{tier_j}]")
    print(f"\n  % LaTeX macros:")
    print(f"  \\newcommand{{\\KhojaGerizimDevKm}}{{{km_g:.2f}}}       % Khoja Ahmed Yasawi deviation from Gerizim (km)")
    print(f"  \\newcommand{{\\KhojaGerizimDevBeru}}{{{dev_g:.6f}}}   % Khoja Ahmed Yasawi deviation from Gerizim (beru)")
    print(f"  \\newcommand{{\\KhojaGerizimTier}}{{{tier_g}}}         % Khoja Ahmed Yasawi tier relative to Gerizim")
    print(f"  \\newcommand{{\\KhojaJerDevKm}}{{{km_j:.2f}}}          % Khoja Ahmed Yasawi deviation from Jerusalem (km)")
    print(f"  \\newcommand{{\\KhojaJerDevBeru}}{{{dev_j:.6f}}}      % Khoja Ahmed Yasawi deviation from Jerusalem (beru)")
    print(f"  \\newcommand{{\\KhojaJerTier}}{{{tier_j}}}            % Khoja Ahmed Yasawi tier relative to Jerusalem")

    # ── A++ count ─────────────────────────────────────────────────────────
    print(f"\n  A++ sites from Gerizim: {n_app}")
    print(f"  \\newcommand{{\\NappTotal}}{{{n_app}}}             % total A++ sites from Gerizim anchor")
    print(f"  \\newcommand{{\\GerizimAppCount}}{{{n_app}}}        % Gerizim A++ count (alias)")

    # ═════════════════════════════════════════════════════════════════════
    print("\n" + "=" * 72)
    print("  GROUP 22 — World Peace Pagoda (Lumbini)")
    print("  Script: analysis/global/anchor_site_comparison.py")
    print("=" * 72)

    dev_wpp, km_wpp, tier_wpp = beru_dev(WPP_LON, GERIZIM)
    print(f"\n  World Peace Pagoda")
    print(f"    lon = {WPP_LON}°E  (Wikidata / OSM coordinate)")
    print(f"    vs Gerizim: δ = {dev_wpp:.6f} beru = {km_wpp:.2f} km  [{tier_wpp}]")
    print(f"\n  % LaTeX macros:")
    print(f"  \\newcommand{{\\WPPDeltaKm}}{{{km_wpp:.2f}}}          % World Peace Pagoda deviation from Gerizim (km)")
    print(f"  \\newcommand{{\\WPPDeltaBeru}}{{{dev_wpp:.6f}}}      % World Peace Pagoda deviation from Gerizim (beru)")

    # ── All A++ sites — top-hit macros (\topHit*) and appendix macros (\appSite*) ──
    print("\n" + "=" * 72)
    print("  GROUP 23 — Top-Hit and AppSite Macros (all A++ sites, sorted by deviation)")
    print("  Script: analysis/global/anchor_site_comparison.py")
    print("=" * 72)
    import re as _re
    for (s, d), ord_label in zip(app_sites_sorted, _ORDINAL_LABELS):
        km_t  = d * BERU * KM_PER_DEG
        m_t   = km_t * 1000
        name  = _html_to_latex(s["name"])
        lon_s = _fmt_lon(s["lon"])
        short_raw = _re.split(r'[,(]', s["name"])[0].strip()
        # Strip HTML italic tags for the short name (used in body text)
        short = _re.sub(r'<[Ii]>(.*?)</[Ii]>', r'\1', short_raw).strip()
        # Strip common lead-in words so body text reads naturally.
        # Only strip if the name is clearly "Prefix + core-name".
        for _prefix in ("Mausoleum of ", "Church of ", "Ruins of ",
                        "Archaeological Site of ", "Monastery of "):
            if short.startswith(_prefix):
                short = short[len(_prefix):]
                break
        # For "X of Y and Z" patterns where X is a known abbreviated form,
        # keep only X.  Example: "Sacri Monti of Piedmont and Lombardy" → "Sacri Monti"
        # We detect this by checking if the name contains " of " and the part
        # before " of " is not itself a generic descriptor (handled above).
        _m = _re.match(r'^([\w\s]+) of (.+)$', short)
        if _m:
            _before = _m.group(1).strip()
            # Only shorten if the 'before' part is ≥2 words (i.e., proper name)
            if len(_before.split()) >= 2:
                short = _before
        print(f"\n  % Rank {ord_label}: {s['name'][:60]}")
        # \topHit* macros (top 3 only — manuscript body text references these)
        if ord_label in ("One", "Two", "Three"):
            print(f"  \\newcommand{{\\topHitName{ord_label}}}{{{name}}}  % rank-{ord_label.lower()} A++ site name")
            print(f"  \\newcommand{{\\topHitLon{ord_label}}}{{{lon_s}}}  % rank-{ord_label.lower()} A++ site longitude (°E)")
            print(f"  \\newcommand{{\\topHitDevBeru{ord_label}}}{{{d:.6f}}}  % rank-{ord_label.lower()} beru deviation")
            print(f"  \\newcommand{{\\topHitDevKm{ord_label}}}{{{km_t:.2f}}}  % rank-{ord_label.lower()} deviation (km)")
            print(f"  \\newcommand{{\\topHitDevM{ord_label}}}{{{int(round(m_t))}}}  % rank-{ord_label.lower()} deviation (m)")
            print(f"  \\newcommand{{\\topHit{ord_label}}}{{{short}}}  % rank-{ord_label.lower()} short name")
        # \appSite* macros — full A++ appendix table (all ranks)
        print(f"  \\newcommand{{\\appSiteName{ord_label}}}{{{name}}}  % A++ table rank {ord_label}: name")
        print(f"  \\newcommand{{\\appSiteLon{ord_label}}}{{{lon_s}}}  % A++ table rank {ord_label}: longitude")
        print(f"  \\newcommand{{\\appSiteDevBeru{ord_label}}}{{{d:.6f}}}  % A++ table rank {ord_label}: beru deviation")
        print(f"  \\newcommand{{\\appSiteDevKm{ord_label}}}{{{km_t:.2f}}}  % A++ table rank {ord_label}: km")
        print(f"  \\newcommand{{\\appSiteDevM{ord_label}}}{{{int(round(m_t))}}}  % A++ table rank {ord_label}: m")

    # ── Jerusalem and Mecca anchor counts ─────────────────────────────────
    print("\n" + "=" * 72)
    print("  GROUP 24 — Jerusalem and Mecca anchor A+ counts")
    print("=" * 72)
    # Jerusalem A+ p-value (self-excluded corpus, N = full corpus - 1)
    from scipy.stats import binomtest as _binom2
    jer_p = _binom2(jer_ap_count, len(sites) - 1, P_NULL_AP, alternative="greater").pvalue
    print(f"  \\newcommand{{\\JerusalemAp}}{{{jer_ap_count}}}  % Jerusalem A+ count (self-excluded corpus)")
    print(f"  \\newcommand{{\\NjerSelfExclCorpus}}{{{len(sites) - 1}}}  % corpus size, Jerusalem self-excluded (Sweep B)")
    print(f"  \\newcommand{{\\NgerizimSweepCorpus}}{{{len(sites) - 1}}}  % corpus size, Gerizim sweep (Jerusalem removed, Sweep A)")
    def _pfmt(p):
        """Format p-value for LaTeX: scientific notation if it would round to zero."""
        if p < 0.0001:
            from decimal import Decimal
            exp = int(f"{p:.2e}".split("e")[1])
            mantissa = p / (10 ** exp)
            return f"{mantissa:.2f} \\times 10^{{{exp}}}"
        return f"{p:.4f}"

    print(f"  \\newcommand{{\\pJerusalemAp}}{{{_pfmt(jer_p)}}}  % p-value, Jerusalem A+ binomial (self-excluded)")
    print(f"  \\newcommand{{\\MeccaAp}}{{{mecca_ap}}}  % Mecca A+ count (full corpus)")
    print(f"  \\newcommand{{\\pMeccaAnchor}}{{{_pfmt(mecca_p)}}}  % p-value, Mecca binomial (A+, one-sided)")
    print(f"  \\newcommand{{\\MeruAp}}{{{meru_ap}}}  % Mt. Meru (Tanzania, 36.750E) A+ count")
    print(f"  \\newcommand{{\\pMeruAnchor}}{{{_pfmt(meru_p)}}}  % p-value, Mt. Meru binomial (A+, one-sided)")
    print(f"  \\newcommand{{\\KailashAp}}{{{kailash_ap}}}  % Mt. Kailash (81.312E) A+ count")
    print(f"  \\newcommand{{\\pKailashAnchor}}{{{_pfmt(kailash_p)}}}  % p-value, Kailash binomial (A+, one-sided)")

    ResultsStore = __import__("lib.results_store", fromlist=["ResultsStore"]).ResultsStore
    ResultsStore().write_many({
        "JerusalemAp":          jer_ap_count,
        "NjerSelfExclCorpus":   len(sites) - 1,
        "NgerizimSweepCorpus":  len(sites) - 1,
        "pJerusalemAp":         jer_p,
        "MeccaAp":              mecca_ap,
        "pMeccaAnchor":         mecca_p,
        "MeruAp":               meru_ap,
        "pMeruAnchor":          meru_p,
        "KailashAp":            kailash_ap,
        "pKailashAnchor":       kailash_p,
    })

    print("\n" + "=" * 72)
    print("  DONE — GROUP 16, 22, 23, 24 macros printed above.")
    print("=" * 72)


if __name__ == "__main__":
    main()
