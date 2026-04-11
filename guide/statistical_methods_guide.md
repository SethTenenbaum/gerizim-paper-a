# Statistical Methods Guide

## A Plain Language Companion to Paper A: The Global 3 Degree Periodic Structure in the UNESCO World Heritage Longitude Distribution

*Tenenbaum (2026). This guide explains every concept, formula, and statistical test used in the manuscript in clear American English prose.*

---

## Table of Contents

1. [What This Study Is About](#what-this-study-is-about)
2. [The Beru Unit and the Harmonic Grid](#the-beru-unit-and-the-harmonic-grid)
3. [Beru Deviation: How Closeness Is Measured](#beru-deviation)
4. [Tier Classification](#tier-classification)
5. [The Geometric Null: The 4 Percent Baseline](#the-geometric-null)
6. [The Circularity Problem and Its Mathematical Consequences](#the-circularity-problem)
7. [Test 1: Full Corpus Enrichment (Descriptive)](#test-1-full-corpus)
8. [Test 2: Dome and Spherical Monument Enrichment (Primary)](#test-2-dome-enrichment)
9. [Test 2b: Hemispherical Mound Evolution (Sensitivity Variant)](#test-2b-mound-evolution)
10. [Test 3: Cluster Asymmetry (Primary)](#test-3-cluster-asymmetry)
11. [Test 4: Temporal Gradient (Descriptive)](#test-4-temporal-gradient)
12. [Test B: Rayleigh Phase Lock (Structural)](#test-b-rayleigh)
13. [Test E: Founding Site Classifier (Post Hoc)](#test-e-founding-site)
14. [The Binomial Enrichment Test](#the-binomial-enrichment-test)
15. [Fisher's Exact Test](#fishers-exact-test)
16. [The Cochran Armitage Trend Test](#cochran-armitage-trend-test)
17. [The Mantel Haenszel Stratified Test](#mantel-haenszel-stratified-test)
18. [Permutation Tests and Z Scores](#permutation-tests)
19. [The Rayleigh Test for Circular Uniformity](#the-rayleigh-test)
20. [Clopper Pearson Exact Confidence Intervals](#clopper-pearson-intervals)
21. [The Mann Whitney U Test](#mann-whitney-u-test)
22. [Simulation Null Models](#simulation-null-models)
23. [Geographic Concentration Null Models](#geographic-concentration-nulls)
24. [Spatial Autocorrelation and Effective N Correction](#spatial-autocorrelation)
25. [Multiple Comparisons: Bonferroni and FDR](#multiple-comparisons)
26. [The Unit Sensitivity Sweep](#unit-sensitivity-sweep)
27. [The Global Anchor Sweep](#global-anchor-sweep)
28. [The x.18 Degree Periodic Structure](#the-x18-periodic-structure)
29. [Leave One Out and Leave K Out Sensitivity](#leave-one-out)
30. [Enrichment Ratio and Effect Size](#enrichment-ratio)
31. [How the Tests Combine: The Convergence Argument](#convergence)
32. [Summary of All Tests at a Glance](#summary-table)
33. [Glossary](#glossary)

---

## 1. What This Study Is About {#what-this-study-is-about}

This study asks a single central question: does the list of UNESCO World Heritage Cultural and Mixed sites show a regular repeating pattern in their longitude positions, and if so, does that pattern match a known ancient unit of angular measurement?

The ancient unit in question is the Babylonian beru, equal to 30 degrees of arc, which is documented in the cuneiform astronomical compendium known as MUL.APIN and dates to roughly 1000 BCE. The study tests whether UNESCO sites tend to cluster at positions that are multiples of 0.1 beru (that is, 3 degrees) from a reference longitude in the Levant.

The study is entirely exploratory. No pre-registration was filed. All results are treated as preliminary pattern documentation, not as confirmed findings. The manuscript makes no claim about the mechanism that might have produced the pattern. It tests the statistics and reports what it finds.

The corpus consists of 1,011 UNESCO Cultural and Mixed World Heritage sites with valid longitude coordinates. The reference anchor is Mount Gerizim at 35.272 degrees east, chosen on historical grounds because of its documented role as the cosmological center point (the "navel of the earth") in the Samaritan tradition.

Three classes of explanation are evaluated throughout the paper. The first is geographic concentration: the idea that civilizations producing UNESCO caliber monuments simply developed in longitude bands that happen to align with the beru grid. The second is metrological transmission: the idea that a 3 degree angular reference system anchored near the Levant was deliberately used in siting major foundations. The third is unmodeled statistical artifact: the idea that some systematic bias in UNESCO curation, geocoding, or the discrete longitude distribution produces the apparent periodicity. Adjudicating between the first two explanations requires independent corpora and is beyond the scope of this study.

---

## 2. The Beru Unit and the Harmonic Grid {#the-beru-unit-and-the-harmonic-grid}

### What is a beru?

A beru is a Babylonian unit of arc equal to 30 degrees. It divides the full circle of 360 degrees into 12 equal segments. The beru is attested in the MUL.APIN compendium, where it serves as both a celestial and terrestrial unit of measurement. The terrestrial use is attested explicitly in MUL.APIN II ii 43 through iii 15, which converts beru of daylight into ground distances. At the equator, 30 degrees of longitude corresponds to roughly 3,330 kilometers.

The beru subdivides into 30 smaller units called ush (written as u-sh in the manuscript), where 1 ush equals 1 degree. The beru also functioned as a unit of territorial description. Royal inscriptions record distances in beru between major cities. On the Babylonian Imago Mundi, distances between the "islands" beyond the encircling Salt Sea are marked as multiples of beru. Boundary stones called kudurru record territorial extents in cuneiform, and many were deposited in temples as legally authoritative records.

### What is a harmonic?

This study tests the 0.1 beru subdivision, which equals 3 degrees or 3 ush. Starting from the anchor longitude of Mount Gerizim (35.272 degrees east), harmonic meridians repeat at every 3 degree interval in both directions. So the harmonic meridians include 32.272 degrees, 35.272 degrees, 38.272 degrees, 41.272 degrees, and so on, wrapping around the globe. There are 120 such harmonics spanning the full 360 degrees.

### Why 0.1 beru?

The 0.1 beru interval (3 degrees) is the smallest integer subdivision of the beru that is attested in cuneiform metrology. It corresponds to exactly 3 ush. The study tested this specific spacing because it is the finest documented subdivision that produces a tractable statistical test with a meaningful null rate of 4 percent. Coarser spacings like 0.5 beru or 1 beru are too broad to produce informative tests. Finer spacings such as 0.02 beru do not produce a signal. The 0.1 beru spacing was not selected after the fact to maximize results.

The 0.05 beru spacing (6 degrees) also reaches significance, but this is structurally expected because every 0.10 beru harmonic is also a 0.05 beru harmonic, making the 0.05 beru result a coarser sub-sampling of the same periodicity rather than a separate signal.

---

## 3. Beru Deviation: How Closeness Is Measured {#beru-deviation}

For every UNESCO site with a known longitude, the analysis computes how close that site falls to the nearest harmonic meridian. This measurement is called the beru deviation, written as the Greek letter delta.

### The calculation

Take a site at longitude L (in degrees). The steps are:

1. Compute the angular distance from the anchor: arc = |L minus 35.272|
2. Convert to beru: beru_value = arc / 30.0
3. Find the nearest 0.1 beru harmonic: nearest = round(beru_value / 0.1) times 0.1
4. Compute the deviation: delta = |beru_value minus nearest|

The result is always a number between 0 and 0.05 beru, because the farthest any point can be from a harmonic is exactly halfway between two harmonics (0.05 beru).

### Converting to kilometers

To translate a beru deviation into an approximate physical distance at the equator, multiply as follows: distance in kilometers = delta times 30 degrees times 111 kilometers per degree, which simplifies to delta times 3,330 km. A deviation of 0.002 beru works out to roughly 6.7 kilometers at the equator. At higher latitudes, the physical distance shrinks because longitude lines converge toward the poles, which means the angular threshold is actually more conservative (harder to meet) at higher latitudes. A geodesic sensitivity check confirms that applying a strict physical threshold would increase the A+ count, strengthening all primary results.

### The haversine formula

For precise distances on the earth's surface, the haversine formula is used:

```
a = sin^2(delta_phi / 2) + cos(phi_1) * cos(phi_2) * sin^2(delta_lambda / 2)
c = 2 * arctan2(sqrt(a), sqrt(1 - a))
d = R * c
```

where phi is latitude, lambda is longitude, and R = 6,371 km is the earth's mean radius. For longitude only comparisons (since beru deviation is a longitude concept), the simplified form is:

```
d approximately equals R * |delta_lambda_rad| * cos(phi_mean)
```

which is what the manuscript uses for approximate kilometer conversions.

### Worked example: Lumbini

Lumbini, the birthplace of the Buddha, has a longitude of 83.276 degrees east.

Arc from Gerizim: 83.276 minus 35.272 = 48.004 degrees.
Beru value: 48.004 / 30.0 = 1.6001 beru.
Nearest 0.1 beru harmonic: 1.6 beru.
Deviation: |1.6001 minus 1.6| = 0.0001 beru (roughly 0.33 kilometers).
Classification: Tier A++ (the closest possible tier).

---

## 4. Tier Classification {#tier-classification}

Sites are sorted into tiers based on how close their beru deviation is to zero, meaning how close they fall to a harmonic. All tier boundaries were defined before any beru computations were run.

| Tier | Deviation threshold | Approximate distance | Geometric null probability |
|------|---------------------|----------------------|----------------------------|
| A++  | 0.0002 beru or less | 0.67 km or less      | 0.4%  |
| A+   | 0.002 beru or less  | 6.7 km or less       | 4.0%  |
| A    | 0.010 beru or less  | 33 km or less        | 20.0% |
| B    | 0.050 beru or less  | 167 km or less       | 100%  |

The A+ tier is the primary analytical threshold used throughout the paper. It is the tier at which enrichment is most robust and the 4 percent null provides a meaningful baseline for testing.

The A++ tier is reported descriptively but is not tested statistically because the expected count under the null is too small for reliable inference.

Kilometer labels are equatorial approximations. At higher latitudes, the physical distances shrink because longitude lines converge, making the angular criterion stricter in physical terms.

---

## 5. The Geometric Null: The 4 Percent Baseline {#the-geometric-null}

### What is a null hypothesis?

Every statistical test needs a baseline assumption about what "no signal" looks like. This is called the null hypothesis. In this study, the null hypothesis is that UNESCO site longitudes have no special relationship to the beru grid. Under this assumption, a site's longitude is effectively random with respect to the harmonic positions.

### Deriving the 4 percent

The A+ tier captures any site whose beru deviation is 0.002 or less. Each harmonic repeats every 0.1 beru. The A+ window extends 0.002 beru on each side of a harmonic, for a total width of 0.004 beru. The fraction of the 0.1 beru interval that falls inside this window is:

```
(2 times 0.002) / 0.1 = 0.004 / 0.1 = 0.04 = 4%
```

An equivalent way to see it is that the maximum possible deviation from a harmonic is half a step, or 0.05 beru. The A+ window captures deviations from 0 to 0.002, which is 0.002 / 0.05 = 4 percent of the half step interval [0, 0.05]. Both derivations give the same result: p0 = 4 percent.

This means that if you placed a site at a completely random longitude, it would land in the A+ window 4 percent of the time purely by chance. This 4 percent figure is the denominator for every enrichment ratio in the paper and the baseline probability for every binomial test.

### Why it matters

If the observed A+ rate in some population significantly exceeds 4 percent, that is evidence of non-random clustering near harmonic positions. The geometric null is derived purely from the geometry of the grid. It does not assume anything about archaeology or ancient knowledge. It simply measures the fraction of "target space" relative to "total space" in the periodic grid.

### Simulation confirmation

The paper does not rely solely on the geometric null. A Gaussian kernel density estimate based on the actual, non-uniform distribution of UNESCO longitudes produces an expected A+ count of roughly 40.4 plus or minus 6.2 (approximately 4.0 percent), confirming that the 4 percent geometric baseline is not materially inflated by the fact that UNESCO sites cluster in certain longitude bands. In addition, the primary tests use simulation based null models (longitude permutation, bootstrap, KDE smoothed null) that preserve the empirical longitude distribution and do not depend on the uniformity assumption at all.

---

## 6. The Circularity Problem and Its Mathematical Consequences {#the-circularity-problem}

### What the circularity is

The study began with the hypothesis that the beru unit might be detectable in the UNESCO longitude data, using Mount Gerizim as the reference anchor. When the global anchor sweep was run, it turned out that Gerizim's longitude falls within the optimal phase band (the x.18 degree phase). If the sweep had revealed a different optimal phase, the paper would not have been written in its current form.

This creates a circularity: the anchor was chosen partly on historical grounds, but the decision to report the results was conditioned on the observation that the anchor happened to align with the optimal phase. There is no way to verify after the fact whether the phase alignment was discovered before or after Gerizim was selected as the anchor. Even a timestamped archive cannot resolve this, because the internal sequence of analytical steps is unverifiable by any outside party. The manuscript states this problem explicitly: "This identification problem cannot be closed by any statistical method applied to this corpus."

### What the circularity contaminates

The circularity affects the following results:

The Levantine corridor statistics (Tests 1, 2, and 3 as computed from the Gerizim anchor), because reporting these results was conditioned on observing the phase alignment.

The x.18 degree framing itself, which was identified in the sweep rather than predicted independently.

The Rayleigh Test B result, whose near perfect R value is a mathematical consequence of the A+ classification criterion, not an independent observation.

All narrative emphasis on the Levantine corridor as historically motivated.

### What the circularity does not contaminate

Certain results are structurally immune to the circularity. The existence of a 3 degree periodic structure in the UNESCO longitude file would be recovered by any researcher running an unanchored 120 phase scan without any reference to Gerizim. The dome enrichment (Test 2) is immune because all 120 x.18 degree anchors produce identical A+ counts for the dome population, so the result does not depend on the choice of anchor. The cluster asymmetry (Test 3) is a full corpus result that requires no keyword filter and does not depend on anchor selection. And the unit sensitivity sweep is immune because the sharp collapse of significance at plus or minus 0.3 percent of the canonical spacing is a property of the corpus structure, not of the anchor choice.

### Mathematical formalization of the circularity

The circularity can be understood formally as a form of outcome dependent reporting, sometimes called post-hoc selection or optional stopping in the statistical literature. The key idea is that the probability of the paper being written is not independent of the observed data, and this dependence biases the reported p values downward.

To make this precise, consider the following setup. Let theta denote the phase of the optimal anchor in the global sweep, and let theta_G denote the phase of Mount Gerizim's longitude modulo 3.0 degrees. Before the sweep was run, the investigator's intention was to use Gerizim as the anchor. The sweep then revealed that theta is approximately equal to theta_G (Gerizim's phase of approximately 2.269 degrees modulo 3.0 is within 0.09 degrees, or about 10 km, of the optimal phase band near 2.18 degrees).

Define two conditional probabilities. Let P(report | theta approximately equals theta_G) be the probability that this paper is written, given that the optimal phase aligned with Gerizim. Let P(report | theta far from theta_G) be the probability of writing this paper in its current form, given that the optimal phase fell elsewhere. These two probabilities are not equal. The first is close to 1 (the alignment confirms the starting hypothesis and motivates a full study). The second is close to 0 (a different optimal phase would have undermined the Gerizim motivated framing).

This inequality means that the act of writing and reporting the paper is itself conditioned on a favorable data outcome. The raw p values computed from the Gerizim anchored analysis do not account for this conditioning. In Bayesian terms, the prior probability of the hypothesis has been inflated by the selection event. In frequentist terms, the effective number of "looks" at the data is larger than 1, because the investigator implicitly tested whether the anchor matched the optimal phase before deciding to proceed.

The magnitude of this bias can be bounded. If we treat the anchor sweep as a single additional test (did Gerizim fall in the optimal band or not?), then the worst case bias is a factor of 2 (equivalent to one additional comparison). If we treat it as selecting one of 120 possible phase positions, the worst case would be a factor of 120, although this overestimates the bias because only anchors near the Levant would have been historically motivated. The Bonferroni correction at k = 2 (the primary family) provides a partial accounting, and the conservative k = 5 bound provides additional margin. However, no finite correction factor can fully remove the circularity because the reporting decision itself is the problem, not the number of tests.

### The protocol for handling the circularity

The paper's protocol addresses this problem through four measures. First, all results are treated as exploratory without exception. Second, the Bonferroni correction is applied as a descriptive bound on effect magnitude rather than a confirmatory threshold. Third, the Levantine corridor statistics are reported as the motivated entry point rather than as independent evidence. Fourth, all downstream sections reference this circularity discussion rather than restating the issue.

### Why this design avoids the worst pitfalls of archaeoastronomy

Prior metrological studies in archaeoastronomy were criticized for three reasons: investigators assembled their own site lists, angular intervals were tuned after the fact to maximize results, and tier boundaries were adjusted to produce significant outcomes. This study avoids all three. The UNESCO site list was compiled by an external body without any reference to this hypothesis. The angular interval (0.1 beru) was selected on metrological grounds from the documented sexagesimal subdivision of the beru. And tier boundaries and keyword sets were specified before any beru distances were computed.

---

## 7. Test 1: Full Corpus Enrichment (Descriptive) {#test-1-full-corpus}

### What it asks

Does the entire UNESCO Cultural and Mixed corpus show an A+ rate above the 4 percent geometric null?

### How it works

The full corpus has 1,011 sites. Under the null hypothesis, the expected number of A+ sites is 1,011 times 0.04 = 40.4. The observed count is 57 A+ sites, giving an A+ rate of 5.6 percent and an enrichment ratio of about 1.38 times the null expectation.

A one sided binomial test asks how likely it is to observe 57 or more A+ sites out of 1,011 if each site had only a 4 percent chance of being A+. The resulting p value is approximately 0.010. The 95 percent Clopper Pearson confidence interval for the true A+ rate is [4.2 percent, 7.1 percent], which excludes the 4 percent null at the lower end.

### Why it is classified as descriptive

Test 1 provides the corpus level backdrop for the stronger results in Tests 2 and 3. It shares A+ sites with Test 3 (every A+ site in Test 1 contributes to the cluster asymmetry calculation), so it is not independent. It is reported as context, not as primary evidence.

### Region stratified robustness

A Mantel Haenszel analysis by UNESCO region produces a common odds ratio and p value that reaches trend level significance. Tested against region specific empirical null rates, the enrichment remains significant. The pre-2000 corpus also retains significance.

### World religion sites

Sites associated with any world religion show an elevated A+ rate, and Christianity associated sites reach individual significance. Buddhism associated sites show directional elevation but do not reach significance. These populations were defined by religion name matching only.

---

## 8. Test 2: Dome and Spherical Monument Enrichment (Primary) {#test-2-dome-enrichment}

### What it asks

Among UNESCO sites whose descriptions contain morphological keywords related to domed or spherical architecture (stupa, stupas, tholos, dome, domed, domes, spherical), is the A+ rate significantly higher than 4 percent?

### How it works

The keyword filter identifies 90 sites in the raw sweep. A context validated version (Test 2x) removes non-architectural matches like Uluru's "rock domes" and Hiroshima's "Genbaku Dome," reducing the population to 83 sites. Of these 83 validated sites, 11 are Tier A+, giving an A+ rate of 13.3 percent and an enrichment of 3.31 times the null expectation.

Three statistical approaches are used. The first is the exact binomial test, which yields p = 0.0005. The second is longitude permutation with 100,000 trials, which yields Z = 2.81 and p = 0.0038. The third is bootstrap from the full corpus with 100,000 trials, which yields p = 0.013. The Bonferroni adjusted p value at k = 2 (the primary family size) is 0.0018, still highly significant.

### The A+ dome sites

The 11 Tier A+ dome and stupa sites are Khoja Ahmed Yasawi (0.2 km from harmonic), Boyana Church (0.3 km), Lumbini (0.8 km), Humayun's Tomb (2.0 km), Pyu Ancient Cities (2.3 km), Silk Roads: Chang'an to Tianshan (2.5 km), Trulli of Alberobello (3.6 km), Old City of Jerusalem (4.1 km), Borobudur (4.3 km), Kizhi Pogost (5.0 km), and Echmiatsin/Zvartnots (6.2 km).

### Why domes specifically?

The dome and stupa morphology was selected because hemispherical and spherical sacred architecture has a documented association with cosmological symbolism across multiple civilizations, including Buddhist stupas, Roman domes, Islamic domes, and tholos tombs. The keyword list was fixed before any beru computations were run. No site was excluded or included based on its beru deviation.

### Context validation (Test 2x)

Because some keyword matches are false positives (for example, "rock domes" in a geological description), a context validation step applies sentence level architectural confirmation using a two gate system. Gate 1 applies 28 negative patterns to reject non-architectural uses. Gate 2 applies 86 positive patterns to require architectural confirmation. This is a post hoc robustness check: the validation rules were developed after the initial sweep, so Test 2x is excluded from the primary Bonferroni family. The fact that removing false positive matches strengthens the signal (p drops from 0.0009 to 0.0005) is consistent with the rejected sites being non-architectural noise that dilutes the true enrichment.

### Geographic concentration null models for domes

Domed monuments are concentrated in the Eurasian corridor, with approximately 85 percent falling between 20 degrees east and 120 degrees east, compared to roughly 69 percent of the full corpus. Two dedicated null models test whether this alone explains the dome enrichment.

Null A (within dome bootstrap) draws 83 longitudes from the dome sites' own longitude list with replacement. It asks whether the dome sites' positions reproduce the observed A+ count when resampled. The result is not significant, meaning the dome sites' own longitude distribution is concentrated near harmonic positions and reproduces the enrichment.

Null C (restricted geographic draw) pools all full corpus sites within plus or minus 5 degrees of any dome site (about 623 sites). It draws 83 longitudes from this restricted pool and asks whether any monuments from that same geographic footprint show comparable enrichment. The pool shows 5.8 percent A+, while dome sites show 13.3 percent A+. The difference is significant (p = 0.003).

Together, these establish that proximity to harmonic longitudes is necessary to explain the dome enrichment (Null A) but that proximity alone is not sufficient (Null C). Dome sites significantly exceed the enrichment of other monuments at the same longitudes. The remaining question is whether this morphological specificity reflects metrological intent or some unmodeled feature of the civilizations that built hemispherical monuments.

---

## 9. Test 2b: Hemispherical Mound Evolution (Sensitivity Variant) {#test-2b-mound-evolution}

### What it asks

If the dome and stupa population is extended to include earthen mound stage structures (tumulus, tumuli, barrow, barrows, kofun, mound), which share hemispherical geometry, does the enrichment persist?

### How it works

The extended population has 104 sites. The combined A+ rate exceeds 4 percent and the result is significant (p less than 0.001). However, the three stages show different profiles. The mound stage alone does not reach significance (6.2 percent A+, p = 0.26, not significant). The stupa stage shows the strongest enrichment (22.2 percent A+, 5.56 times, p = 0.005). The dome stage falls between the two. The signal is carried by stupas and domes, not by mounds.

A context validated version (Test 2bx) reduces the population to 126 sites and yields stronger results (p = 0.0005, 2.78 times enrichment).

### Stupa geographic concentration null

Stupa sites are concentrated in South and Southeast Asia. The same Null A and Null C framework is applied. The within stupa bootstrap (Null A) is not significant, meaning stupa positions themselves reproduce the count. The restricted geographic draw (Null C) is significant, meaning other monuments at the same longitudes are less enriched. The structure parallels the dome result.

### Why it is classified as a sensitivity variant

Test 2b shares the entire dome and stupa core with Test 2 and adds only a subset of mound sites. It is not an independent hypothesis. It is not counted separately in the Bonferroni family. It tests whether the morphological boundary of the dome enrichment extends to related forms, and the answer is that it does for the combined population but not for mounds alone.

---

## 10. Test 3: Cluster Asymmetry (Primary) {#test-3-cluster-asymmetry}

### What it asks

Do A+ sites tend to be located in denser clusters of UNESCO sites than non-A+ sites? And do the harmonics that contain at least one A+ site host more total sites than harmonics with no A+ sites?

### How it works

For every UNESCO site, the analysis counts how many other sites fall within plus or minus 0.3 degrees of longitude. This produces a "cluster size" for each site. A+ sites have a higher mean cluster size than non-A+ sites. A Mann Whitney U test confirms this difference (one tailed), and a permutation test with 10,000 shuffles provides a non-parametric confirmation (Z and p value reported).

At the harmonic level, the result is even stronger. The analysis counts the total number of UNESCO sites assigned to each 0.1 beru harmonic. Harmonics that contain at least one A+ site host roughly 3 times as many total sites as harmonics without any A+ sites (p less than 0.0001).

A further breakdown compares sites at harmonic node positions to those at harmonic midpoints (anti-nodes), finding that harmonic nodes have significantly more neighbors.

### Subtest: does clustering produce the unit sensitivity peak?

A joint permutation test asks whether extreme clustering mechanically produces the canonical unit peak. The answer is no: among permutations that reproduce the clustering pattern, the fraction that also produces the sensitivity peak is nearly identical to the unconditional rate. The unit sensitivity peak and the clustering signal are statistically independent.

### Why it matters

This test uses the full corpus with no keyword filter. It asks a structurally different question from Test 2: not "are certain types of sites near harmonics?" but "are harmonic positions themselves located in areas where UNESCO sites are dense?" The fact that both tests are significant provides converging evidence from two different angles. The two primary tests are pre-specified as a family of k = 2.

---

## 11. Test 4: Temporal Gradient (Descriptive) {#test-4-temporal-gradient}

### What it asks

Does the A+ rate decline as UNESCO inscription year increases? The reasoning is that the earliest inscribed sites (1978 to 1984) represent the most universally recognized ancient monuments, and if these monuments were placed with reference to a harmonic grid, they should show higher A+ rates than sites inscribed more recently.

### How it works

Sites are grouped into three cohorts: 1978 to 1984 (founding canon), 1985 to 1999, and 2000 to 2025. The A+ rate in the founding canon is approximately 8 percent, in the middle cohort approximately 6.8 percent, and in the modern era approximately 4.5 percent (indistinguishable from the null). A Cochran Armitage trend test confirms a statistically significant monotonic decline (Z approximately negative 1.72, p = 0.043, one sided). A five cohort version also reaches significance.

A Fisher exact test comparing the founding canon (1978 to 1984) directly with the modern era (2000 to 2025) yields a trend level result.

### Inscription year stratification

The five cohort breakdown reveals that only the founding canon (1978 to 1984) reaches individual significance. The 1985 to 1991 cohort breaks the otherwise monotonic decline, likely reflecting a disproportionate share of Mixed category sites in that period. Post-2000 cohorts are indistinguishable from the 4 percent null.

### Sequential significance

When the corpus is reconstructed in the order UNESCO inscribed it, the cumulative A+ rate reaches significance by 1982 and remains significant through the full corpus. Peak significance occurs around 2003, and then gradually declines as more modern (null rate) sites dilute the signal.

### Why it is classified as descriptive

A within region diagnostic reveals that the temporal gradient is largely explained by regional composition. The founding canon is concentrated in Europe and the Mediterranean, where the beru harmonic signal is strongest. Later cohorts add proportionally more sites from Africa, the Americas, and Oceania, where the null hypothesis predicts lower A+ rates regardless of any metrological effect. After region adjustment using the Mantel Haenszel method (common OR and chi squared statistic reported), the gradient weakens. No individual region replicates the pooled gradient at conventional significance.

The temporal gradient is therefore treated as descriptive context, not as independent evidence. The substantive observation is that UNESCO's own ordering of importance correlates with the beru harmonic signal, but this correlation is substantially explained by the geographic composition changing over time.

---

## 12. Test B: Rayleigh Phase Lock (Structural) {#test-b-rayleigh}

### What it asks

Do the A+ sites share a common phase (angular position) relative to the 3 degree grid?

### Why it is excluded from all family counts

The near perfect Rayleigh R value of the A+ sub-sample is a mathematical consequence of the A+ classification criterion. Because A+ sites are defined as those within 0.002 beru of a harmonic, they are all, by definition, clustered at the same phase position (near zero deviation from the harmonic). This means the Rayleigh test's result is guaranteed by the classification rule and does not provide independent evidence of periodicity. It is reported as a structural description, not as a test of the hypothesis. Accordingly, it is excluded from all Bonferroni family counts and the FDR audit.

A permutation check confirms that a random draw of the same number of longitudes from the full corpus does not reproduce the observed phase concentration, and the full corpus Rayleigh R is small with a nonsignificant permutation p value, confirming no spurious global periodicity in the raw longitude file.

---

## 13. Test E: Founding Site Classifier (Post Hoc) {#test-e-founding-site}

### What it asks

Are sites classified as founding capitals, sacred origins, or founding monuments by a keyword based classifier over-represented among A+ sites?

### How it works

The founding site classifier assigns sites to five categories: Founding Capital, Sacred Origin, Founding Monument, Founding Axis, and Ancient Continuous Landscape, using unambiguous and context gated keywords. Applied to the full corpus, the classifier finds an elevated A+ rate among keyword founding sites compared to non-founding sites (Fisher exact test OR greater than 1, p approximately 0.011).

### Why it is classified as post hoc

The founding site classifier in the code file lib/founding_filter.py was iteratively tuned after examining misclassified sites. Because the classifier was refined after data inspection, the enrichment result is directly affected by this tuning and must be treated as fully post hoc. It carries no weight in the primary evidential chain. The two primary results (Tests 2 and 3) and all keyword independent sub-tests are unaffected by this classifier. The complete version history of the classifier is preserved in the project's git log.

---

## 14. The Binomial Enrichment Test {#the-binomial-enrichment-test}

### What it does

The binomial test is the most basic statistical tool in this study. It asks: given a population of n sites, is the observed number of A+ sites (k) significantly higher than what you would expect if each site had only a 4 percent chance of being A+?

### The mathematical model

Under the null hypothesis, each site is independently either A+ or not A+, with probability p0 = 0.04 of being A+. The number of A+ sites in a group of n follows a binomial distribution: Binomial(n, p0).

The one sided p value is the probability of observing k or more A+ sites purely by chance:

```
P(X >= k | n, p0) = sum from j = k to n of C(n, j) * p0^j * (1 - p0)^(n - j)
```

where C(n, j) is the binomial coefficient "n choose j," meaning the number of ways to select j items from n.

In Python, this is computed as: `scipy.stats.binomtest(k, n, p0, alternative='greater').pvalue`

### Why one sided?

All binomial tests in this manuscript are one sided (testing for "greater than") because the hypothesis is directional. The beru hypothesis predicts that certain sub-populations should be over-represented in the A+ tier, not under-represented. A two sided test would waste statistical power testing for depletion, which is not what the hypothesis predicts.

### Example: dome population

For 83 dome/stupa sites with 11 observed A+ sites and a null probability of 0.04:

Expected count: 83 times 0.04 = 3.32.
Observed count: 11.
Enrichment ratio: 11 / 3.32 = 3.31 times.
p value: P(X >= 11 | 83, 0.04) = 0.0005.

The probability of getting 11 or more A+ hits out of 83 sites by chance alone is 0.05 percent, which is very strong evidence against the null hypothesis for this sub-population.

### Building intuition

Under the null, the expected count is 3.32 and the standard deviation is the square root of (83 times 0.04 times 0.96), which is approximately 1.79. Observing 11 is (11 minus 3.32) divided by 1.79 = 4.3 standard deviations above the mean. That is deep in the tail of the distribution. If you dealt 83 cards from a deck where 4 percent are marked, getting 11 marked cards would be like being dealt a full house three times in a row: theoretically possible but very improbable.

---

## 15. Fisher's Exact Test {#fishers-exact-test}

### What it does

Fisher's exact test compares A+ rates between two groups without relying on large sample approximations. It is used in this study when comparing founding class sites versus non-founding sites, dome sites versus non-dome sites within a geographic band, and different inscription era cohorts against each other.

### The setup

A two by two contingency table:

```
                  A+      Not A+     Total
Group 1            a        b         n1
Group 2            c        d         n2
Total              k       N - k      N
```

Fisher's test computes the exact probability of observing this table (and all more extreme tables) under the hypergeometric distribution, assuming the row and column totals are fixed. The probability of a specific table is:

```
P(table | marginals) = C(n1, a) * C(n2, c) / C(N, k)
```

The one sided p value sums over all configurations with a greater than or equal to the observed value, holding marginals fixed.

In Python: `scipy.stats.fisher_exact([[a, b], [c, d]], alternative='greater')`

### The odds ratio

The test also produces an odds ratio: OR = (a times d) / (b times c). An odds ratio greater than 1 means Group 1 has a higher A+ rate than Group 2. An odds ratio of 2.0, for example, means the first group is about twice as likely to contain an A+ site.

### Why Fisher instead of chi squared?

Fisher's exact test is preferred here because some cell counts in the two by two tables are small (less than 5), where the chi squared approximation becomes unreliable. Fisher's test gives exact probabilities based on the hypergeometric distribution, with no large sample assumption required.

### When to use Fisher versus binomial

The binomial test asks whether a group's A+ rate exceeds the 4 percent geometric null (a fixed theoretical baseline). Fisher's test compares two groups within the data against each other, without assuming the 4 percent baseline. These are complementary: the binomial measures absolute enrichment above the geometric expectation, while Fisher measures relative enrichment between two observed groups.

---

## 16. The Cochran Armitage Trend Test {#cochran-armitage-trend-test}

### What it does

The Cochran Armitage test checks whether the A+ rate changes monotonically across ordered groups. In this study, it tests whether the A+ rate decreases as UNESCO inscription year increases.

### Why not just compare pairs?

Testing each pair of cohorts separately would require multiple comparison corrections and would throw away the ordering information. The Cochran Armitage test uses the ordering directly, making it more powerful for detecting monotonic trends with a single test rather than multiple pairwise comparisons.

### The mechanics

Sites are divided into K ordered cohorts (for example, three inscription year groups scored 1, 2, 3). Let n_i be the number of sites in cohort i, a_i the number of A+ sites, and s_i the ordinal score. Define:

```
N = sum of all n_i (total sites)
R = sum of all a_i (total A+ sites)
p_bar = R / N (overall A+ rate)
s_bar = sum of (s_i * n_i) / N (weighted mean score)
```

The numerator of the test statistic captures the covariation between the ordinal score and the A+ count:

```
T = sum of (s_i * a_i) minus R * s_bar
```

The variance of T under the null hypothesis (no trend) is:

```
Var(T) = p_bar * (1 minus p_bar) * [sum of (s_i^2 * n_i) minus N * s_bar^2]
```

The standardized Z statistic is:

```
Z = T / sqrt(Var(T))
```

Under the null, Z follows a standard normal distribution. For a decreasing trend (the A+ rate declines as the score increases), the one sided p value is the left tail probability:

```
p = Phi(Z)
```

where Phi is the standard normal cumulative distribution function. Because a decreasing trend produces Z less than 0, and Phi of a negative number is less than 0.5, the p value is small when the decreasing trend is strong.

### In this study

Three cohort result: Z = negative 1.72, p = 0.043 (one sided). The A+ rate drops from roughly 8 percent in the founding canon to roughly 4.5 percent in the modern era. The five cohort version also reaches significance. However, after region stratification, the gradient is largely explained by the changing geographic composition of the corpus over time.

---

## 17. The Mantel Haenszel Stratified Test {#mantel-haenszel-stratified-test}

### What it does

The Mantel Haenszel test asks whether an enrichment effect persists after stratifying by a potential confounding variable, in this case geographic region. If the A+ enrichment were purely an artifact of European sites being concentrated at certain longitudes, it should disappear within individual regions. If it persists across regions, it cannot be explained by geography alone.

### How it works

For each of K geographic strata (for example, UNESCO administrative regions crossed with longitude bins), a separate two by two table of A+ versus not A+ is constructed. The Mantel Haenszel estimator computes a common odds ratio that pools across strata:

```
OR_MH = sum of (a_i * d_i / N_i) divided by sum of (b_i * c_i / N_i)
```

This weights each stratum by its contribution to the overall information, down-weighting strata with small samples.

Under the null hypothesis (common OR = 1 in every stratum), the expected value of a_i in each stratum is:

```
E(a_i) = n1_i * k_i / N_i
```

The test statistic is:

```
chi_squared_MH = [sum of (a_i minus E(a_i))]^2 / [sum of Var(a_i)]
```

which under the null follows a chi squared distribution with one degree of freedom.

### In this study

The Mantel Haenszel test is used in two main contexts. First, it is used to test whether the temporal gradient persists after stratifying by region. When the temporal gradient is tested with region as a stratifying variable, the adjusted OR weakens toward 1, indicating that regional composition is the dominant driver. This is why the temporal gradient is treated as descriptive.

Second, it is used in the full corpus analysis to assess whether the overall A+ enrichment persists across regions. The region adjusted common OR is reported alongside a chi squared statistic and p value.

### Key assumptions

The Mantel Haenszel test assumes stratum homogeneity (the true odds ratio is the same in all strata), independence of observations within a stratum, and sufficiently large strata for the chi squared approximation.

---

## 18. Permutation Tests and Z Scores {#permutation-tests}

### What they do

Permutation tests provide a non-parametric null distribution that makes no assumptions about how site longitudes are distributed. They are the primary evidence against the chance explanation for Tests 2 and 3.

### The procedure

1. Compute the observed test statistic (for example, the count of A+ dome sites).
2. For each of B permutation trials (100,000 for the simulation nulls, 10,000 for other permutation checks), shuffle all 1,011 corpus longitudes and assign them to the dome indexed positions. Count how many A+ sites result from the shuffled data.
3. The permutation p value is the fraction of trials where the shuffled statistic equals or exceeds the observed statistic.
4. The permutation Z score standardizes the observed statistic against the null distribution: Z = (observed minus mean of null) / (standard deviation of null).

### Why permutation instead of the geometric null?

The geometric null (4 percent) assumes site longitudes are uniformly distributed with respect to the beru grid. In reality, UNESCO sites cluster heavily in Europe and Asia, where beru harmonics intersect certain longitude bands more densely. The permutation test uses the actual UNESCO longitude distribution as its null, making it robust to this geographic confound. This is why the permutation p value is typically more conservative (larger) than the binomial p value.

### Key result

For the dome population: observed = 11 A+ sites, permutation null mean = approximately 3.3, Z = 2.81, p = 0.0038. The permutation p value is more conservative than the binomial p value (0.0038 versus 0.0005) because the permutation accounts for the non-uniform longitude distribution. Even so, the dome enrichment survives comfortably.

### Other permutation tests in the paper

Permutation tests are also used for the cluster asymmetry (Test 3), the Rayleigh phase lock verification (Test B), the circular shift anchor null (Test C in the x.18 degree formal tests), the sensitivity slope sharpness test (Test I), and the canonical unit specificity test (Test II). In each case, the same logic applies: shuffle the data to build a null, then compare the observed statistic to the shuffled distribution.

---

## 19. The Rayleigh Test for Circular Uniformity {#the-rayleigh-test}

### What it does

The Rayleigh test checks whether a set of angular values cluster in a preferred direction rather than being uniformly distributed around a circle. In this study, it is applied to the phases of site longitudes modulo 3 degrees (the harmonic period).

### The statistic

For each site, the longitude is converted to a phase angle theta_j within the 3 degree period. The mean resultant length R measures the concentration of the angles:

```
R = (1 / n) * |sum of e^(i * theta_j)|
```

where i is the imaginary unit and the sum runs over all n sites. R = 0 means perfectly uniform (no preferred phase). R = 1 means all angles are identical (perfect phase lock). The Rayleigh test statistic is Z = n times R squared, which under the null (uniform distribution on the circle) has an approximately exponential distribution. The p value is:

```
p approximately equals exp(negative Z) * [1 + correction terms involving Z and n]
```

### In this study

For the A+ sub-sample, R is near 1. However, this is a mathematical consequence of the A+ definition: all A+ sites have deviations less than 0.002 beru, so they all cluster at the same phase position (near zero deviation from the harmonic). The Rayleigh test's result is guaranteed by the classification rule and does not provide independent evidence. It is excluded from all family counts.

For the full corpus, R is small and the permutation p value is not significant, confirming that there is no spurious global periodicity in the raw longitude file. This is an important negative result because it means the 3 degree signal is specific to the A+ tier and not an artifact of the overall longitude distribution.

---

## 20. Clopper Pearson Exact Confidence Intervals {#clopper-pearson-intervals}

### What they are

A confidence interval gives a range of plausible values for the true A+ rate in a population, based on the observed data. The Clopper Pearson method is an exact method (no large sample approximation) based on inverting the binomial cumulative distribution function. It is conservative, meaning it never under-covers the true rate.

### How it works

Given k successes in n trials, the 95 percent Clopper Pearson interval [L, U] is defined by:

The lower bound L is the value of p such that P(X >= k | n, p) = 0.025.
The upper bound U is the value of p such that P(X <= k | n, p) = 0.025.

In Python: `scipy.stats.binom.ppf` and related functions are used to compute the bounds.

### In this study

For the full corpus (57 A+ out of 1,011), the 95 percent interval is [4.2 percent, 7.1 percent]. This interval excludes 4 percent at the lower end, consistent with the significant binomial test. For the dome population (11 out of 83), the interval is wider, reflecting the smaller sample size, but still excludes 4 percent.

These intervals are reported throughout the manuscript alongside p values to provide a measure of precision and to show how far above the null the true rate plausibly lies.

---

## 21. The Mann Whitney U Test {#mann-whitney-u-test}

### What it does

The Mann Whitney U test compares two distributions without assuming normality. In this study, it compares the cluster sizes (number of neighboring sites within plus or minus 0.3 degrees longitude) of A+ sites versus non-A+ sites.

### How it works

The test ranks all cluster sizes from both groups together, then sums the ranks for one group. If one group tends to have larger cluster sizes, its rank sum will be larger than expected by chance. The test statistic U captures this difference, and the p value is computed from the distribution of U under the null hypothesis that the two groups are drawn from the same population.

### Why it is used

Cluster sizes are not normally distributed (they are right skewed because some longitude bands are very dense and others are nearly empty). The Mann Whitney test is robust to non-normality because it works with ranks rather than raw values. It is the appropriate non-parametric alternative to the two sample t-test.

### In this study

A+ sites have significantly higher mean cluster sizes than non-A+ sites (one tailed test). This is confirmed by both the Mann Whitney U test and by an independent permutation test, providing converging non-parametric evidence that A+ sites tend to sit in the densest parts of the longitude distribution.

---

## 22. Simulation Null Models {#simulation-null-models}

### What they are

The study uses three simulation approaches, each run for 100,000 trials with a fixed random seed (seed = 42 for reproducibility), to build null distributions that are more realistic than the geometric 4 percent baseline. These simulation nulls preserve the actual longitude distribution of UNESCO sites rather than assuming uniform placement.

### Longitude permutation

All 1,011 corpus longitudes are shuffled randomly. For each trial, the dome indexed positions (the positions in the shuffled list that correspond to the dome sites in the original list) are read off, and the A+ count is computed. This preserves the total set of longitudes but randomizes which sites occupy which longitude positions. It tests whether dome sites are more enriched than expected given the full corpus longitude distribution, including its Eurasian concentration.

For the dome population: p_perm = 0.0038 (significant).
For the founding canon: p = 0.14 (not significant).
For the pre-2000 corpus: p = 0.074 (trend level).

### Bootstrap from full corpus

For each trial, 83 longitudes (matching the dome population size) are drawn randomly with replacement from the full 1,011 site longitude list. This tests whether a random sub-sample of the same size would show comparable enrichment. For the dome population: p_boot = 0.013 (significant).

### KDE smoothed null

A Gaussian kernel density estimate of the full corpus longitude distribution is used to generate synthetic longitudes. This produces a smoother null that removes any fine scale discreteness in the observed distribution. The expected count is 40.4 plus or minus 6.2 (approximately 4.0 percent), confirming the geometric null. The observed 57 A+ sites yields Z = 2.50 and p = 0.011.

### Key result

The dome enrichment survives all three simulation null models. The canon and pre-2000 signals are more sensitive to the non-uniform longitude distribution and are treated as less robust. The fact that the dome signal persists even when the empirical longitude distribution is used as the null is the central statistical finding of the paper.

---

## 23. Geographic Concentration Null Models {#geographic-concentration-nulls}

### The problem

Domed monuments are concentrated in the Eurasian corridor. About 85 percent fall between 20 degrees east and 120 degrees east, compared to roughly 69 percent of the full corpus. Could this geographic concentration alone explain their high A+ rate?

### Null A: within dome bootstrap

Draw 83 longitudes from the dome sites' own longitude list (with replacement). This tests whether the dome sites' positions reproduce the observed A+ count when resampled. Result: the p value is not significant. The bootstrap mean matches the observed count. This tells us that the dome sites' own geographic distribution is concentrated near harmonic positions, and simply resampling from those positions reproduces the enrichment. Geographic concentration is necessary to explain the result.

### Null C: restricted geographic draw

Pool all full corpus sites within plus or minus 5 degrees of any dome site (about 623 sites). Then draw 83 longitudes from this restricted pool. This tests whether any monument from the dome longitude footprint shows comparable enrichment. Result: the pool shows 5.8 percent A+, while dome sites show 13.3 percent A+. The difference is significant (p = 0.003, Z score reported). Sensitivity checks at plus or minus 2 degrees and plus or minus 10 degrees yield consistent results. Geographic concentration is not sufficient.

### Reconciliation

Together, Null A and Null C establish a clear picture. Null A shows why dome sites are near harmonics: their own longitude positions are concentrated there. Null C shows that other monuments at those same longitude positions are significantly less enriched. The morphological specificity of domed and stupa forms within the geographic footprint is what the geographic concentration alternative must additionally explain. Answering this requires independent corpora.

### Region conditioned permutation

A within group permutation test shuffles sites within UNESCO administrative region crossed with 20 degree longitude bins. In every trial, the observed A+ count is reproduced exactly. This is mathematically inevitable: a bin of 5 degrees or wider spans at least one full 3 degree harmonic period, so shuffling within the bin leaves each site's harmonic residue unchanged. Within-bin permutation therefore cannot discriminate a real alignment from geographic concentration at this resolution. The geographic concentration question is properly addressed by Null A and Null C, not by within-bin permutation.

---

## 24. Spatial Autocorrelation and Effective N Correction {#spatial-autocorrelation}

### The problem

UNESCO sites are not spatially independent. A cluster of European sites at similar longitudes will all tend to have similar beru deviations. This violates the independence assumption of the binomial test and can inflate significance.

### Design effect correction (DEFF)

The design effect measures how much spatial clustering inflates the variance compared to independent sampling. It is computed using the intraclass correlation coefficient (ICC) within longitude bands of plus or minus 0.5 degrees. The effective sample size is N divided by DEFF, which is smaller than the raw N. The corrected p value using the effective N still excludes the null.

Formally, the design effect is:

```
DEFF = 1 + rho * (m_bar minus 1)
```

where rho is the ICC and m_bar is the average cluster size. The effective sample size is:

```
N_eff = N / DEFF
```

### Block bootstrap

As an independent check, the corpus is divided into non-overlapping 3 degree blocks (matching the harmonic period). Within each block, sites are resampled with replacement 100,000 times. The 95 percent confidence interval from this block bootstrap excludes the 4 percent null, confirming that spatial autocorrelation does not explain the enrichment.

### Why this matters

If spatial autocorrelation were severe enough, the apparent enrichment could be an artifact of non-independent sampling. The DEFF correction and block bootstrap provide two independent ways to check this. Both confirm that the signal survives correction for spatial clustering.

---

## 25. Multiple Comparisons: Bonferroni and FDR {#multiple-comparisons}

### The problem

When you run multiple statistical tests, some will appear significant by chance. If you run 20 tests at a 5 percent significance level, you expect about 1 false positive even if there is no real signal. The study addresses this with two approaches.

### Bonferroni correction

The Bonferroni correction divides the significance threshold by the number of tests. The primary pre-specified family consists of k = 2 tests (Test 2 for dome enrichment and Test 3 for cluster asymmetry). The corrected threshold is 0.05 / 2 = 0.025. Both primary tests pass this threshold.

The adjusted p value for each test is: p_adj = raw p times k. If p_adj is less than 0.05, the result survives the correction. Dome enrichment: p_adj = 0.0018 (passes comfortably). Cluster asymmetry: p_adj = 0.0068 (passes). Both also survive the more conservative k = 5 bound (dome p_adj = 0.005, cluster p_adj = 0.017).

### Benjamini Hochberg false discovery rate (FDR)

The BH-FDR procedure controls the expected proportion of false discoveries among all significant results, rather than controlling the probability of any false positive. It is less conservative than Bonferroni and is the standard approach when many tests are reported. The procedure works by sorting all p values in ascending order, then comparing each to a threshold that increases linearly with rank:

```
p_(i) <= (i / m) * q
```

where m is the total number of tests and q is the target false discovery rate (typically 0.05). Applied across all reported tests, the FDR audit retains the primary results at q less than 0.05.

The full FDR audit covers all reported tests (excluding Test B, which is structural). The number of tests surviving FDR at q less than 0.05 and the number surviving Bonferroni at k equal to the total number of tests are both reported in the manuscript.

### Why k = 2 rather than a larger number?

Tests 2 and 3 ask structurally different questions on different populations (dome sub-population versus full corpus clustering). They share no A+ sites in common (Test 2 examines dome A+ sites; Test 3 examines the spatial distribution of all sites around all harmonics). The remaining tests are either descriptive (Test 1, Test 4), sensitivity variants (Test 2b), post hoc (Test E), or structurally guaranteed (Test B). Including these in the primary family would be over-conservative and would not reflect the actual structure of the hypothesis space. The k = 5 bound is also reported for transparency.

---

## 26. The Unit Sensitivity Sweep {#unit-sensitivity-sweep}

### What it does

The harmonic spacing is varied from 0.05 beru to 0.20 beru while the physical hit window is held fixed at 6.7 kilometers. For each spacing, the A+ rate is computed for both the dome and full corpus populations. This tests whether the 0.1 beru spacing is special or whether any nearby spacing would produce comparable results.

### Key finding

The 0.1 beru spacing produces the strongest joint significance across both populations. Adjacent non-harmonic spacings (0.06, 0.07, 0.09, 0.11, 0.12 beru) are not significant in either population. The 0.05 and 0.20 beru spacings also reach significance, but this is structurally expected: every 0.10 beru harmonic is also a 0.05 beru harmonic and a 0.20 beru harmonic, so these are sub-samplings of the same periodicity rather than separate signals. The 0.08 beru spacing reaches nominal significance in the dome population but not in the full corpus, and inspection shows nearly all of its dome hits are a subset of the 0.10 beru A+ dome sites.

### Fine grained sensitivity

Within plus or minus 1 percent of 0.10 beru, a minus 0.5 percent shift (0.0995 beru) collapses the dome population to non-significance, and a plus or minus 0.3 percent shift collapses both populations. Only the canonical value and its immediate neighbor (0.1001 beru) retain dome significance. This is an extremely sharp peak.

### Permutation context

Two permutation tests evaluate whether the sharp peak at the canonical unit is unusual.

Test I (sharpness permutation) repeats the fine sweep on thousands of longitude permuted datasets. The sharp peak pattern occurs with moderate frequency in permuted data, making the test nonsignificant. Geographic concentration is sufficient to reproduce the collapse pattern.

Test II (canonical unit specificity) asks whether 0.1000 beru specifically achieves the best joint significance among all fine sweep spacings. The fraction of permutations where 0.1000 is the unique best joint significance spacing is not significant. Geographic concentration alone produces the canonical spacing as the best spacing in a non-trivial fraction of permutations.

Both permutation tests are non-significant. The unit sensitivity observation is consistent with the metrological hypothesis but does not independently confirm it. The 0.1 beru harmonic interval derives from the documented sexagesimal subdivision of the Babylonian beru, not from the data.

---

## 27. The Global Anchor Sweep {#global-anchor-sweep}

### What it does

Every longitude from 0 to 360 degrees (at 0.01 degree resolution, totaling 36,000 trial anchors) is tested as a candidate reference point. For each trial anchor, the A+ count is computed for the full corpus. Two parallel sweeps are run: one with Jerusalem removed (Sweep A), and one with Jerusalem retained but self-excluded when it serves as the anchor (Sweep B).

### Key finding

Gerizim ranks at approximately the 97th percentile (top 3 percent). Jerusalem produces the same 56 A+ count. The single site difference between sweeps (55 versus 56 A+) is precisely the Old City of Jerusalem, which sits at 0.0012 beru from Gerizim's harmonic. The enrichment belongs to the roughly 35 degree east corridor, not to the specific choice of anchor within it. No other historical meridian reaches significance except Mecca, which shows trend level results but does not reach conventional significance.

### Landmark ranking

When each of the 1,011 UNESCO sites is used as an anchor, Gerizim (an external anchor not in the inscribed corpus) ranks in the top 4 percent. Excluding sites at the x.18 degree phase, Gerizim ranks number 24 of 965 (top 2.5 percent).

### Global corridor comparison

The 35 degree east corridor leads globally, with the most sites scoring at or above 55 A+ when used as beru anchors. The closest rival corridor (23.2 degrees east) spans four countries with no shared cosmological tradition.

---

## 28. The x.18 Degree Periodic Structure {#the-x18-periodic-structure}

### What it is

The most striking quantitative result of the global anchor sweep is that every longitude of the form x.18 degrees east (all 120 such anchors) produces identically the maximum A+ count. The spacing between consecutive x.18 degree anchors is exactly 3 degrees = 0.1 beru. The label "x.18 degrees east" is a rounded shorthand for a phase value of approximately 2.18 degrees modulo 3.0 degrees. Mount Gerizim's longitude modulo 3.0 degrees falls within combined coordinate and bandwidth uncertainty of this optimal phase.

### Two interpretations

Interpretation 1 (beru alignment): a global grid with 3 degree spacing anchored to the approximately 35 degree east corridor would produce exactly this pattern.

Interpretation 2 (geographic concentration): Eurasian heritage sites cluster in longitude bands that happen to align with the x.18 degree phase. This cannot be refuted on the UNESCO corpus alone.

### Formal tests

Three permutation tests (10,000 shuffles) are applied. Test A checks the full corpus Rayleigh R modulo 3 degrees and finds no spurious global periodicity (not significant). This is important because it rules out the possibility that the 3 degree signal is an artifact of the overall longitude distribution. Test B checks the A+ only Rayleigh R and finds it significant, but this is structurally expected because of the A+ classification criterion. Test C uses a circular shift anchor null and finds that a randomly placed anchor recovers only about 40 A+ sites on average versus 56 at Gerizim (significant).

### Dome specific periodicity audit

A separate dome specific audit confirms that Gerizim scores at the top percentile among dome sites. The dome site phases are approximately uniform (no evidence they cluster at a single phase more than the full corpus), and all dome A+ sites have phases within plus or minus 0.15 degrees of the optimal phase window. Within the x.18 degree window, dome sites and non-dome sites both show elevated A+ rates, consistent with both the geographic concentration and beru alignment readings.

### Non-independence

The 3 degree periodicity and the unit sensitivity sweep are two views of a single pattern and do not constitute independent evidence. The periodicity tells you there is a 3 degree structure. The unit sweep tells you it is sharpest at 0.1 beru. These are the same observation expressed in different ways.

---

## 29. Leave One Out and Leave K Out Sensitivity {#leave-one-out}

### What they do

These checks test whether a single influential site drives the entire dome result. In a leave one out analysis, each of the 11 A+ dome sites is removed in turn, and the test is rerun. If the result survives every single omission, no single site is responsible for the finding.

### Results

Leave one out: removing any single A+ site yields 10 A+ out of 82, with p approximately 0.002. The result survives every omission. Leave two out: the worst case across all 55 pair combinations still yields a significant result. Leave three out: the worst case across all 165 triple combinations also survives.

These checks confirm that the dome enrichment is distributed across multiple independent sites and is not driven by any small subset. This is important because it rules out the possibility that a single anomalous site (say, one with an incorrectly geocoded longitude) could be responsible for the entire signal.

---

## 30. Enrichment Ratio and Effect Size {#enrichment-ratio}

### What the enrichment ratio is

The enrichment ratio is the observed A+ rate divided by the null expectation (4 percent). It provides a simple, intuitive measure of effect size. An enrichment ratio of 1.0 means the observed rate exactly matches the null. A ratio of 3.31 (the dome result) means 3.31 times more A+ sites than expected by chance.

### Enrichment ratios in this study

The full corpus enrichment ratio is 1.38, which is modest. It means there are approximately 15 additional A+ sites above chance expectation in 1,011, consistent with a scenario in which only a fraction of sites were placed with reference to a grid. The dome population (3.31 times) and stupa population (5.56 times) show much stronger effect sizes, with the caveat that smaller populations have wider confidence intervals.

### Interpreting small versus large enrichment

A small enrichment ratio (like 1.38 for the full corpus) can still be statistically significant if the sample is large enough, because the test becomes more sensitive with more data. Conversely, a large enrichment ratio (like 5.56 for stupas) may not be highly significant if the sample is small (only 9 stupa sites). Both the enrichment ratio and the p value should be considered together.

---

## 31. How the Tests Combine: The Convergence Argument {#convergence}

### The logic

No single test in this study is decisive. The strength of the analysis comes from the convergence of independent lines of evidence:

First, the full corpus shows an A+ rate modestly above the null (Test 1, descriptive).

Second, dome and spherical monuments, identified by a mechanical keyword filter, show an A+ rate more than 3 times the null expectation, surviving all three simulation nulls and leave-k-out checks (Test 2, primary).

Third, A+ bearing harmonics host roughly 3 times as many total sites as non-A+ harmonics, using the full corpus with no keyword filter (Test 3, primary).

Fourth, the earliest inscribed sites show the highest A+ rates, though this is largely explained by regional composition (Test 4, descriptive).

Fifth, the 0.1 beru spacing uniquely produces joint significance in both populations, while all neighboring non-harmonic spacings collapse (unit sensitivity sweep), though permutation tests show geographic concentration can reproduce this pattern.

Sixth, the geographic concentration null models show that proximity to harmonic longitudes is necessary but not sufficient to explain the dome enrichment (Null A and Null C).

Seventh, the Americas sub-corpus shows directional depletion (2.1 percent A+, not significant), consistent with a Eurasian anchored grid.

### What the convergence does not prove

The convergence does not prove that ancient societies used a beru grid. It establishes a reproducible quantitative regularity whose spacing matches a known metrological unit. The geographic concentration alternative remains viable because the distribution of Eurasian civilizations and the beru harmonics from the Levant are not independent. This study offers the result as a quantitative basis for further investigation.

The two findings most resistant to a purely geographic account are the morphological specificity (dome enrichment exceeding the restricted geographic draw, Null C) and the cluster asymmetry (Test 3). Together with the global 3 degree periodicity, these constitute a coherent exploratory pattern consistent with a 0.1 beru harmonic grid, awaiting independent historical confirmation.

### Criteria for confirmation

Independent replication using a non-UNESCO corpus would provide the most direct test. The hypothesis would gain substantially from cuneiform or other texts explicitly assigning beru longitude coordinates to distant sites, from archaeological survey instruments calibrated in beru subdivisions, or from textual evidence from Hellenistic, Sanskrit, or Islamic sources referencing a Levantine centered angular reference frame.

---

## 32. Summary of All Tests at a Glance {#summary-table}

| Test | Description | Population | Status | Key statistic |
|------|-------------|------------|--------|---------------|
| Test 1 | Full corpus A+ rate | 1,011 sites | Descriptive | Binomial p = 0.010 |
| Test 2 | Dome/spherical A+ rate | 83 to 90 sites | Primary | Permutation p = 0.0038, enrichment 3.31x |
| Test 2b | Hemispherical mound evolution | 104 sites | Sensitivity variant | p < 0.001 |
| Test 2x | Context validated dome | 83 sites | Post hoc robustness | p = 0.0005, enrichment 3.31x |
| Test 3 | Cluster asymmetry | 1,011 sites | Primary | Permutation p = 0.0034, harmonic ratio 3x |
| Test 4 | Temporal gradient | 1,011 sites | Descriptive | Cochran Armitage p = 0.043 |
| Test B | Rayleigh phase lock | 57 A+ sites | Structural (excluded) | R near 1 by construction |
| Test E | Founding site classifier | 1,011 sites | Post hoc | Fisher p = 0.011 |
| Test I | Sensitivity sharpness permutation | Permutation null | Robustness | Not significant |
| Test II | Canonical unit specificity | Permutation null | Robustness | Not significant |
| Null A | Within dome bootstrap | 83 dome sites | Geographic control | Not significant |
| Null C | Restricted geographic draw | ~623 sites | Geographic control | p = 0.003 |
| Stupa Null A | Within stupa bootstrap | Stupa sites | Geographic control | Not significant |
| Stupa Null C | Restricted geographic draw for stupas | Restricted pool | Geographic control | Significant |

Primary Bonferroni family: Tests 2 and 3 (k = 2). Both survive at p_adj < 0.025. Both also survive the conservative k = 5 bound.

Full FDR audit: across all reported tests, the primary results survive BH-FDR at q < 0.05, and the strongest results survive Bonferroni at k = total number of tests.

---

## 33. Glossary {#glossary}

**A+ tier**: The primary analytical threshold. A site is A+ if its beru deviation is 0.002 or less (roughly within 6.7 km of a harmonic at the equator).

**Anchor**: The reference longitude from which beru distances are measured. Mount Gerizim at 35.272 degrees east is the primary anchor.

**Beru**: A Babylonian unit of arc equal to 30 degrees, attested in the MUL.APIN cuneiform compendium.

**Beru deviation (delta)**: The angular distance between a site's beru value and the nearest 0.1 beru harmonic. Always between 0 and 0.05 beru.

**Binomial test**: Tests whether the count of A+ sites in a group exceeds the expectation from the geometric null (4 percent).

**Block bootstrap**: A resampling method that respects spatial structure by dividing the data into non-overlapping 3 degree blocks matching the harmonic period.

**Bonferroni correction**: Divides the significance threshold by the number of tests to control the family-wise error rate.

**Bootstrap**: Drawing random samples with replacement from an observed dataset to build a null distribution.

**Circularity**: The problem that the decision to report results was conditioned on the observation that the anchor longitude fell in the optimal phase band.

**Clopper Pearson interval**: An exact confidence interval for a binomial proportion, computed by inverting the binomial cumulative distribution function.

**Cochran Armitage test**: Tests for a monotonic trend in proportions across ordered groups.

**Context validation (Test 2x)**: A post hoc check that removes non-architectural keyword matches from the dome population.

**DEFF (design effect)**: The factor by which spatial clustering inflates variance compared to independent sampling.

**Enrichment ratio**: The observed A+ rate divided by the null expectation (4 percent). A ratio of 3x means three times more A+ sites than expected by chance.

**FDR (false discovery rate)**: The expected proportion of false positives among all significant results. Controlled by the Benjamini Hochberg procedure.

**Fisher's exact test**: Compares A+ rates between two groups using the hypergeometric distribution, without large sample approximations.

**Geometric null**: The 4 percent probability that a randomly placed site falls within the A+ window, derived purely from the grid geometry.

**Harmonic**: A longitude that falls exactly on a multiple of 0.1 beru (3 degrees) from the anchor.

**Haversine formula**: The standard formula for computing great circle distances between two points on a sphere.

**ICC (intraclass correlation coefficient)**: Measures the similarity of beru deviations among sites within the same longitude band. Used to compute the design effect.

**KDE (kernel density estimate)**: A smoothing method used to estimate the continuous probability distribution of UNESCO site longitudes.

**Leave one out (LOO)**: Removing one site at a time to check that no single site drives the result.

**Mann Whitney U test**: A non-parametric test comparing two distributions by their ranks.

**Mantel Haenszel test**: Tests for an association in stratified two by two tables, pooling information across strata.

**MUL.APIN**: A Babylonian cuneiform compendium of astronomical observations, dated to roughly 1000 BCE.

**Null hypothesis (H0)**: The assumption that there is no beru signal, meaning site longitudes are random with respect to the harmonic grid.

**Null A (within dome bootstrap)**: A geographic concentration null model that resamples from the dome sites' own longitudes. Tests whether geographic concentration is necessary.

**Null C (restricted geographic draw)**: A geographic concentration null model that resamples from all sites near dome site longitudes. Tests whether geographic concentration is sufficient.

**Odds ratio (OR)**: The ratio of odds of A+ status in one group versus another. OR greater than 1 means the first group is more likely to be A+.

**p value**: The probability of observing a result as extreme as the actual result, assuming the null hypothesis is true. Smaller p values provide stronger evidence against the null.

**Permutation test**: A non-parametric method that shuffles the data to build a null distribution, making no assumptions about the underlying distribution.

**Phase**: The position of a longitude within the 3 degree harmonic period, computed as longitude modulo 3.0 degrees.

**Rayleigh test**: Tests whether circular data (angles) cluster in a preferred direction rather than being uniformly distributed.

**Ush**: A Babylonian sub-unit of the beru. 1 ush = 1 degree. 30 ush = 1 beru.

**x.18 degrees east**: A shorthand for the optimal phase band. Every longitude whose value modulo 3.0 is approximately 2.18 degrees produces the maximum A+ count. The spacing between consecutive x.18 degree anchors is exactly 3 degrees, equaling 0.1 beru.
