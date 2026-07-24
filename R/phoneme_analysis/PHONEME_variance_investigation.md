# Phoneme variance investigation — `[1.5]` / `[1.6]` driver tables

Investigation of the unexpected variance introduced by individual phonemes in the
`[1.5]_PHONEME_drivers.R` output tables, plus the `[1.6]` unrelated-skew side
question. All figures reconstruct the IDF-weighted cosine pipeline from
`RUHLENdf_PH.csv` + `phoneme_freq_ruhlen_austronesian.csv` and cross-check against
the stored `data/PHONEME_cossim.csv` values.

**Baseline means (from `PHONEME_cossim.csv`, 58 PH languages):**

| baseline | mean cosine to PH |
|---|---|
| Japanese | **0.159** |
| Spanish | 0.117 |
| English | 0.095 |
| Unrelated null | 0.087 |

Japanese sits above every other baseline *and* the unrelated null — that gap is the
central anomaly.

---

## Q0 — All phonemes present in the Philippine language set

Of the full 728-phoneme Ruhlen inventory, only **49 are present in at least one of the
58 PH languages** (`ph_prevalence > 0`); the other 679 (93%) never occur in any PH
language and are irrelevant to the PH-side of every cosine. All 49, descending by
`ph_prevalence`:

| ipa | class | ph_prevalence | IDF | austronesian_prevalence | unrelated_prevalence | span_has | jap_has | eng_has |
|---|---|---|---|---|---|---|---|---|
| ɡ | consonant | 1.000 | 0.649 | 0.521 | 0.607 | ✓ | ✓ | ✓ |
| d | consonant | 1.000 | 0.571 | 0.563 | 0.668 | ✓ | ✓ | ✓ |
| b | consonant | 1.000 | 0.547 | 0.577 | 0.682 | ✓ | ✓ | ✓ |
| j | consonant | 1.000 | 0.411 | 0.662 | 0.953 | ✓ | ✓ | ✓ |
| w | consonant | 1.000 | 0.236 | 0.789 | 0.763 | ✓ | ✓ | ✓ |
| ŋ | consonant | 1.000 | 0.111 | 0.894 | 0.488 | | | ✓ |
| n | consonant | 1.000 | 0.018 | 0.982 | 0.986 | ✓ | ✓ | ✓ |
| t | consonant | 1.000 | 0.014 | 0.986 | 0.981 | ✓ | ✓ | ✓ |
| k | consonant | 1.000 | 0.014 | 0.986 | 0.976 | ✓ | ✓ | ✓ |
| m | consonant | 1.000 | 0.007 | 0.993 | 0.972 | ✓ | ✓ | ✓ |
| ʔ | consonant | 0.983 | 0.616 | 0.539 | 0.417 | | | |
| l | consonant | 0.983 | 0.155 | 0.856 | 0.825 | ✓ | | |
| s | consonant | 0.966 | 0.139 | 0.870 | 0.872 | ✓ | ✓ | ✓ |
| a | vowel | 0.966 | 0.054 | 0.947 | 0.877 | ✓ | ✓ | ✓ |
| i | vowel | 0.931 | 0.050 | 0.951 | 0.863 | ✓ | ✓ | |
| p | consonant | 0.879 | 0.115 | 0.891 | 0.896 | ✓ | ✓ | ✓ |
| u | vowel | 0.776 | 0.092 | 0.912 | 0.834 | ✓ | | |
| h | consonant | 0.603 | 0.499 | 0.606 | 0.635 | | ✓ | ✓ |
| o | vowel | 0.534 | 0.268 | 0.764 | 0.678 | ✓ | ✓ | |
| **ɨ** (barred i) | vowel | 0.345 | 1.74 | 0.173 | 0.142 | | | |
| **vː** (long V) | mod_vowel | 0.345 | 0.857 | 0.423 | 0.569 | | ✓ | |
| **ɾ** (tap) | consonant | 0.293 | 1.63 | 0.194 | 0.199 | ✓ | ✓ | |
| r (trill) | consonant | 0.293 | 0.577 | 0.560 | 0.545 | ✓ | | |
| e | vowel | 0.293 | 0.385 | 0.680 | 0.640 | ✓ | ✓ | |
| **ə** (schwa) | vowel | 0.276 | 1.34 | 0.261 | 0.256 | | | ✓ |
| **cː** (geminate C) | mod_consonant | 0.190 | 2.25 | 0.102 | 0.246 | | ✓ | |
| ɸ | consonant | 0.121 | 2.94 | 0.049 | 0.052 | | | |
| ɔ | vowel | 0.121 | 1.43 | 0.236 | 0.322 | | | ✓ |
| ɛ | vowel | 0.121 | 1.40 | 0.243 | 0.351 | | | ✓ |
| ɣ | consonant | 0.103 | 1.53 | 0.215 | 0.128 | | | |
| β | consonant | 0.086 | 1.99 | 0.134 | 0.043 | | | |
| ɡ̌ | consonant | 0.086 | 1.91 | 0.144 | 0.417 | | ✓ | ✓ |
| I | vowel | 0.069 | 2.71 | 0.063 | 0.147 | | | ✓ |
| f | consonant | 0.069 | 1.45 | 0.232 | 0.374 | ✓ | | ✓ |
| ž (ʒ) | consonant | 0.052 | 3.86 | 0.018 | 0.218 | | | ✓ |
| c | consonant | 0.052 | 2.56 | 0.074 | 0.100 | | | |
| ñ (n˜) | consonant | 0.052 | 1.66 | 0.187 | 0.393 | ✓ | | |
| ɟ | consonant | 0.034 | 3.86 | 0.018 | 0.057 | | | |
| ʊ | vowel | 0.034 | 2.94 | 0.049 | 0.137 | | | ✓ |
| x | consonant | 0.034 | 2.61 | 0.070 | 0.332 | ✓ | | |
| č (tʃ) | consonant | 0.034 | 1.74 | 0.173 | 0.668 | ✓ | ✓ | ✓ |
| v | consonant | 0.034 | 1.11 | 0.327 | 0.346 | | | ✓ |
| ẓ | consonant | 0.017 | 4.96 | 0.004 | 0.014 | | | |
| v̆ | mod_vowel | 0.017 | 4.96 | 0.004 | 0.000 | | | |
| βʷ | consonant | 0.017 | 4.55 | 0.007 | 0.000 | | | |
| ḷ | consonant | 0.017 | 4.27 | 0.011 | 0.062 | | | |
| ɐ | vowel | 0.017 | 3.35 | 0.032 | 0.019 | | | ✓ |
| ø | vowel | 0.017 | 3.17 | 0.039 | 0.085 | | | |
| cʷ | mod_consonant | 0.017 | 1.28 | 0.275 | 0.052 | | | |

**Reading the table.** The top 19 rows (ph_prevalence ≥ 0.53) are all low-IDF (<0.65) —
these are the "backbone" phonemes examined in the weighting experiment. Below that, the
table drops sharply into a long tail of 30 phonemes each present in fewer than 35% of PH
languages, where IDF climbs quickly (up to 4.96) — this is where nearly every phoneme
flagged elsewhere in this report lives: the Q1 Japanese drivers (**ɨ, vː, ɾ, ə, cː,
ɡ̌**), the Q2 Spanish loan phonemes (**x, f, ñ, č**), and the Q3 unrelated-tail drivers
(**ž, ḷ, ẓ, ɟ, ø**) are all packed into this bottom 30. Only 13 of the 49 are shared by
all three baselines simultaneously (m, t, k, n, a, p, s, w, j, b, d, ɡ, č); the
marked/rare half of the table is where baseline-to-baseline differences — and the
anomalies investigated below — actually originate.

---

## Q1 — Japanese: what drives the high similarity scores

Reconstructed mean cos(JP, PH) = **0.1586** (matches stored). The score is carried by
a handful of phonemes — **top 5 = 71%, top 10 = 95%** of the total:

| rank | ipa | class | IDF | ph_prevalence | austronesian_prevalence | contrib | share_of_cos |
|---|---|---|---|---|---|---|---|
| 1 | **cː** (geminate C) | mod_consonant | 2.25 | 0.190 | 0.102 | 0.0344 | **21.7%** |
| 2 | **ɾ** (tap) | consonant | 1.63 | 0.293 | 0.194 | 0.0302 | **19.1%** |
| 3 | ɡ | consonant | 0.649 | 1.000 | 0.521 | 0.0193 | 12.2% |
| 4 | d | consonant | 0.571 | 1.000 | 0.563 | 0.0150 | 9.5% |
| 5 | b | consonant | 0.547 | 1.000 | 0.577 | 0.0137 | 8.7% |
| 6 | **vː** (long V) | mod_vowel | 0.857 | 0.345 | 0.423 | 0.0117 | 7.4% |
| 7 | ɡ̌ | consonant | 1.91 | 0.086 | 0.144 | 0.0091 | 5.7% |
| 8 | j | consonant | 0.411 | 1.000 | 0.662 | 0.0078 | 4.9% |
| 9 | h | consonant | 0.499 | 0.603 | 0.606 | 0.0067 | 4.2% |
| 10 | w | consonant | 0.236 | 1.000 | 0.789 | 0.0026 | 1.6% |

**Finding.** ~**29% of Japanese's similarity comes from two collapsed length /
gemination binaries** (`cː` geminate consonant + `vː` long vowel) and another
**19% from the tap `ɾ`**. All three are high-IDF, so a single shared yes/no flag is
worth roughly a fifth of the whole score each. The match is typologically
coincidental — Japanese moraic length/gemination and a Japanese tap-/r/ align against
the 19–35% of PH languages that happen to trip the same flags. The voiced stops
ɡ/d/b contribute only through sheer ubiquity (ph_prevalence = 1.0, low IDF). The
`cː`/`vː` drivers are the same collapsed-modification binaries examined in Q4.

### The GMM fit for Japanese

The same `[2]` bimodality test run on Japanese picks up a much larger influenced
group — **20 of 58 PH languages (34.5%)**, nearly double Spanish's 10:

*Isnag, Yogad, North Kalinga, South Kalinga, Balangaw, Central Bontok, Inibaloi,
Itbayaten, Ivatanen, Kapampangan, Iraya, Ata, Maranao, Hanunoo, Buhid, Batak, Naga,
Kalagan, Hiligaynon, Kuyunon*

| group | n | mean cossim_jap | sd cossim_jap | mean cossim_unr |
|---|---|---|---|---|
| Non-influenced | 38 | 0.100 | 0.051 | 0.080 |
| **Influenced** | 20 | **0.271** | 0.055 | 0.099 |

A 2.7× gap between the two groups on the real cosine scale — comparable in ratio to
Spanish's 2.5×, but spread across twice as many languages. Unlike Spanish, the two
groups' null means are somewhat more separated (0.080 vs 0.099), so part of the
Japanese-influenced group's edge is also riding a slightly richer unrelated-null
baseline, not purely Japanese-specific similarity. **9 of the 10 Spanish-influenced
languages** — all except *Tiruray* — are also flagged Japanese-influenced
(*Itbayaten, Kapampangan, Iraya, Maranao, Hanunoo, Buhid, Batak, Hiligaynon, Kuyunon*),
i.e. almost the entire Spanish-influenced set doubles as Japanese-influenced, consistent
with the [1.6] finding that similarity leverage keeps tracing back to the *same* handful
of marked phonemes (`cː`, `ɾ`, `vː`) regardless of which baseline is being tested against.

---

## Q2 — Spanish: do the expected loan phonemes have PH prevalence?

Of Spanish's 25 phonemes, **23 appear in ≥1 PH language; only 2 are wholly absent**
(ʎ, ʮ). To probe the *loan* signal directly, `ph_inf_prevalence` gives the prevalence
among the **10 Spanish-influenced PH languages** (the high-similarity component from
the Gaussian-mixture classification in `[2]_PHONEME_cosine_distribution_analysis.R` —
not `[1.5]`, correcting an earlier mislabel — `span_influenced == TRUE`):
*Itbayaten, Kapampangan, Iraya, Maranao, Hanunoo, Buhid, Batak, Hiligaynon, Kuyunon,
Tiruray*. `inf_lift = ph_inf_prevalence − ph_prevalence` (positive = enriched in the
influenced subset).

### The GMM fit behind the classification

`[2]` fits a 1–3 component Gaussian mixture (BIC-selected) per baseline — it fits on
`cossim_<baseline> − cossim_unr` internally, but the group means below are reported
on the **actual `cossim_span` scale** (real cosine similarity, not the difference
transform):

| group | n | mean cossim_span | sd cossim_span | mean cossim_unr |
|---|---|---|---|---|
| Non-influenced | 48 | **0.092** | 0.039 | 0.086 |
| **Influenced** | 10 | **0.235** | 0.047 | 0.091 |

The two components are cleanly separated on the real similarity scale: **non-influenced
languages average 0.092**, **influenced languages average 0.235** — a 2.5× gap. The two
groups' own unrelated-null means are nearly identical (0.086 vs 0.091), confirming the
split is driven by *Spanish* similarity itself, not by those 10 languages having an
unusually low null to begin with.

`share_of_cos` below decomposes mean cos(Spanish, PH) = **0.1169** (reconstruction
matches the stored mean) the same way Q1 decomposed Japanese's score —
`contrib_p = mean over ph( U[Spanish,p]·U[ph,p] )`, `share_of_cos_p = contrib_p / Σcontrib`.

| ipa | class | IDF | ph_prevalence | **ph_inf_prevalence** | inf_lift | austronesian_prev | unrelated_prev | **share_of_cos** |
|---|---|---|---|---|---|---|---|---|
| ʎ (ll) | consonant | 5.65 | 0.000 | 0.0 | 0.000 | 0.000 | 0.076 | 0.0% |
| ʮ | consonant | 2.88 | 0.000 | 0.0 | 0.000 | 0.053 | 0.076 | 0.0% |
| **x** | consonant | 2.61 | 0.034 | **0.2** | **+0.166** | 0.070 | 0.332 | 5.2% |
| č (tʃ) | consonant | 1.74 | 0.034 | 0.0 | −0.034 | 0.173 | 0.668 | 2.2% |
| ñ (n˜) | consonant | 1.66 | 0.052 | 0.1 | +0.048 | 0.187 | 0.393 | 3.3% |
| **ɾ** (tap) | consonant | 1.63 | 0.293 | **0.9** | **+0.607** | 0.194 | 0.199 | **25.5%** |
| **f** | consonant | 1.45 | 0.069 | 0.1 | +0.031 | 0.232 | 0.374 | 5.3% |
| ɡ | consonant | 0.649 | 1.000 | 1.0 | 0.000 | 0.521 | 0.607 | **16.3%** |
| r (trill) | consonant | 0.577 | 0.293 | 0.1 | **−0.193** | 0.560 | 0.545 | 3.7% |
| d | consonant | 0.571 | 1.000 | 1.0 | 0.000 | 0.563 | 0.668 | 12.6% |
| b | consonant | 0.547 | 1.000 | 1.0 | 0.000 | 0.577 | 0.682 | 11.6% |
| j | consonant | 0.411 | 1.000 | 1.0 | 0.000 | 0.662 | 0.953 | 6.5% |
| e | vowel | 0.385 | 0.293 | 0.2 | −0.093 | 0.680 | 0.640 | 1.7% |
| o | vowel | 0.268 | 0.534 | 0.4 | −0.134 | 0.764 | 0.678 | 1.5% |
| w | consonant | 0.236 | 1.000 | 1.0 | 0.000 | 0.789 | 0.763 | 2.2% |
| l | consonant | 0.155 | 0.983 | 1.0 | +0.017 | 0.856 | 0.825 | 0.9% |
| s | consonant | 0.139 | 0.966 | 1.0 | +0.034 | 0.870 | 0.872 | 0.7% |
| p | consonant | 0.115 | 0.879 | 0.9 | +0.021 | 0.891 | 0.896 | 0.5% |
| u | vowel | 0.092 | 0.776 | 0.8 | +0.024 | 0.912 | 0.834 | 0.3% |
| a | vowel | 0.054 | 0.966 | 1.0 | +0.034 | 0.947 | 0.877 | 0.1% |
| i | vowel | 0.050 | 0.931 | 1.0 | +0.069 | 0.951 | 0.863 | 0.1% |
| n | consonant | 0.018 | 1.000 | 1.0 | 0.000 | 0.982 | 0.986 | 0.0% |
| t | consonant | 0.014 | 1.000 | 1.0 | 0.000 | 0.986 | 0.981 | 0.0% |
| k | consonant | 0.014 | 1.000 | 1.0 | 0.000 | 0.986 | 0.976 | 0.0% |
| m | consonant | 0.007 | 1.000 | 1.0 | 0.000 | 0.993 | 0.972 | 0.0% |

**Finding.** Unlike Japanese (where two collapsed length binaries plus a single tap
account for ~50% of the score), Spanish's similarity is dominated by the tap **`/ɾ/`
alone at 25.5% of the cosine** — the same phoneme the `ph_inf_prevalence` column
already flagged as the real loan marker (0.29 → 0.90 in the influenced subset) — plus
the near-universal `ɡ/d/b` trio contributing another 40.5% through sheer ubiquity
rather than any Spanish-specific signal. The other loan phonemes (`x`, `f`, `ñ`, `č`)
each contribute only 2–5%, consistent with their low, contact-driven prevalence.

**Findings.**

- The classic loan phonemes are present but at the **low minority prevalence contact
  borrowing predicts** — `/f/` 7%, `/tʃ/` 3%, `/x/` 3%, `/ɲ/` 5%, `/ʎ/` 0%. They are
  *not* spuriously widespread across PH, so Spanish's similarity is built from shared
  common phonemes, not loan-phoneme contamination.
- The **influenced subset localises the loan signal cleanly**: the tap **`/ɾ/` leaps
  from 0.29 → 0.90 (inf_lift +0.61)** and **`/x/` from 0.03 → 0.20 (+0.17)** among the
  10 influenced languages. These two are the genuine Spanish loan markers.
- Counter-intuitively, the **trill `/r/` drops** in the influenced subset (0.29 → 0.10,
  −0.19), as do mid vowels `/e/` (−0.09) and `/o/` (−0.13). So the "Spanish influence"
  the GMM picks up is the **tap and velar fricative**, not the trill or the 5-vowel
  system — worth noting, since the trill is the more stereotyped Spanish import.

---

## Q3 — Unrelated: high-IDF phonemes with >50% PH prevalence

**The requested set is empty.** There are **0 phonemes with IDF > 1 and
ph_prevalence > 0.5**; the highest-IDF phoneme that reaches a PH majority is `/ɡ/`
at IDF 0.649. Every majority-in-PH phoneme is low-IDF, so the weighting is *not*
producing a "common phoneme carrying huge weight" pathology.

Two structural facts underlie this:

- The **IDF distribution is degenerate**: **546 of 728 columns (75%) are pinned at the
  maximum IDF 5.652** — phonemes absent from all 284 Austronesian reference languages —
  and **only 49 columns have any PH prevalence at all**.
- Unrelated-null inflation therefore comes from *moderate*-IDF, *moderate*-prevalence
  phonemes, which the table's own `null_skew_leverage` score is built to catch.

### `null_skew_leverage` — what it measures and which phonemes it flags

`null_skew_leverage` estimates a phoneme's (proportional) contribution to the **third
central moment (right-skew)** of a PH language's cosine-to-unrelated null. A shared
phoneme adds `IDF²` to the cosine numerator; over the unrelated set that term is
`Bernoulli(q)` with `q = unrelated_prevalence`, giving a skew contribution of

```
null_skew_leverage = IDF⁶ · q(1−q)(1−2q)     (gated to ph_prevalence > 0)
```

It is **positive for q < 0.5** (peaks near q ≈ 0.21, a moderate minority), **zero at
q = 0** (a phoneme in ~no unrelated language can't form a tail) and **at q = 0.5**.
The `drives_unrelated_tail` flag fires above the **90th percentile of candidate
phonemes** (leverage cutoff = **59.76**).

**Exactly 5 phonemes are flagged** (`drives_unrelated_tail == TRUE`):

| ipa | class | IDF | ph_prevalence | unrelated_prevalence | **null_skew_leverage** | flag |
|---|---|---|---|---|---|---|
| **ž** (ʒ) | consonant | 3.86 | 0.052 | 0.218 | **318.2** | both (baseline + unrelated-tail) |
| **ḷ** (retroflex l) | consonant | 4.27 | 0.017 | 0.062 | **306.4** | unrelated-tail driver |
| **ẓ** (retroflex z) | consonant | 4.96 | 0.017 | 0.014 | **203.4** | unrelated-tail driver |
| **ɟ** (voiced palatal stop) | consonant | 3.86 | 0.034 | 0.057 | **157.2** | unrelated-tail driver |
| **ø** (front rounded V) | vowel | 3.17 | 0.017 | 0.085 | **65.4** | unrelated-tail driver |

Full flag tally across all 728 phonemes:

| flag | n |
|---|---|
| neither | 677 |
| baseline driver | 46 |
| unrelated-tail driver | 4 |
| both (baseline + unrelated-tail) | 1 |

**Finding.** The tail drivers are all **rare, high-IDF marked consonants** (ž, ḷ, ẓ, ɟ)
plus one front-rounded vowel (ø): present in only 1–5% of PH languages but shared with
a moderate minority of the unrelated controls, so their `IDF⁶` weight pumps the null's
right tail. **`ž` is the single "both" phoneme** — it inflates a baseline's observed
similarity *and* the null that similarity is tested against, making it the one to
scrutinise most.

---

## Q4 — Section 4 table: has the Ruhlen coding collapse introduced variance?

Yes, measurably. `PHONEME_driver_phoible.csv` compares each phoneme's IDF under Ruhlen
vs PHOIBLE. **93 distinct PHOIBLE modified segments collapse into 19 Ruhlen
modification binaries** (mean 4.9 segments/flag, max 16). Largest collapses:

| Ruhlen binary | class | # PHOIBLE segments folded in | PHOIBLE IDF spread erased | median idf_shift |
|---|---|---|---|---|
| dental | mod_consonant | 16 | 0.396 | — |
| long | mod_vowel | 13 | 0.564 | −0.167 |
| creaky | mod_consonant | 10 | 0.000 | — |
| long | mod_consonant | 10 | 0.183 | −0.201 |
| nasalized | mod_vowel | 9 | 0.107 | — |
| syllabic | mod_consonant | 6 | 0.000 | — |
| palatalized | mod_consonant | 5 | 0.107 | −0.172 |
| labialized | mod_consonant | 3 | 0.183 | −0.591 |
| breathy | mod_vowel | 3 | 0.000 | −0.065 |

Two distinct effects:

1. **Variance destruction.** Up to **0.56 of normalized IDF spread** among the folded
   segments (long vowels) is flattened to one value. A rare pharyngealized-dental and
   a common plain-dental become the same binary — the modification's internal
   informativeness is erased.
2. **Systematic under-weighting.** `idf_shift = IDF_ruhlen − IDF_phoible` is **negative
   far more often than positive** — matched consonants 55 negative vs 19 positive;
   mod_consonant 18 neg vs 2 pos; mod_vowel 10 neg vs 6 pos. Because a Ruhlen flag
   fires if *any* qualifying segment is present, the collapsed flag is more prevalent
   than any single PHOIBLE segment → **lower IDF** → the modification is treated as more
   common and weighted down relative to per-segment coding.

Direction of `variance_shift_dir` across all matched rows:

| class | + (Ruhlen higher) | − (Ruhlen lower) | 0 |
|---|---|---|---|
| consonant | 19 | 55 | 15 |
| mod_consonant | 2 | 18 | 1 |
| mod_vowel | 6 | 10 | 0 |
| vowel | 6 | 16 | 8 |

**Finding.** The collapse resolves the Japanese paradox from Q1: `cː`/`vː` remain
high-IDF *within Austronesian* (genuinely rare there) so they dominate the cosine —
while being coarse "does this language have *any* geminate / *any* long vowel"
binaries that lump typologically distinct length phenomena into one heavily-weighted
yes/no. That coarseness is the mechanism inflating Japanese's score.

---

## Side question — `[1.6]` `PHONEME_unrelated_skew_candidates.csv`

Top skew-driving unrelated languages (most often exceeding a PH language's own
`median + 2·sd` cutoff): **Kanuri (34), Fur (31), Bodo (25), Turkmen (21), Kodagu (21),
Kanakuru (18), Micmac (18), Dongolawi (16), Coptic (16), Yaaku (16), Alabama (16)** —
African, Dravidian, Turkic, and Algonquian languages.

Aggregating each PH↔(top-skew-unrelated) cosine phoneme by phoneme, the drivers of the
right skew are:

| ipa | class | IDF | ph_prevalence | unrelated_prevalence | mean_contrib |
|---|---|---|---|---|---|
| **cː** (geminate C) | mod_consonant | 2.25 | 0.190 | 0.246 | 0.0218 |
| **ə** (schwa) | vowel | 1.34 | 0.276 | 0.256 | 0.0178 |
| ɡ | consonant | 0.649 | 1.000 | 0.607 | 0.0171 |
| **ɨ** (barred i) | vowel | 1.74 | 0.345 | 0.142 | 0.0146 |
| b | consonant | 0.547 | 1.000 | 0.682 | 0.0135 |
| d | consonant | 0.571 | 1.000 | 0.668 | 0.0133 |
| **vː** (long V) | mod_vowel | 0.857 | 0.345 | 0.569 | 0.0131 |
| j | consonant | 0.411 | 1.000 | 0.953 | 0.0098 |
| ɾ (tap) | consonant | 1.63 | 0.293 | 0.199 | 0.0069 |
| ɡ̌ | consonant | 1.91 | 0.086 | 0.417 | 0.0067 |
| ɸ | consonant | 2.94 | 0.121 | 0.052 | 0.0056 |
| ž (ʒ) | consonant | 3.86 | 0.052 | 0.218 | 0.0033 |
| ɔ | vowel | 1.43 | 0.121 | 0.322 | 0.0023 |

**Finding.** The skew is driven by **central / high vowels the vowel-rich unrelated
languages share with a PH minority — schwa `ə`, barred-i `ɨ`, `ɔ` — plus the same
length/gemination flags (`cː`, `vː`), the tap `ɾ`, and marked fricatives (`ɸ`, `ž`).**
The central vowels `ə` and `ɨ` are the distinctive additions relative to the Japanese
profile.

---

## Bottom line

The **same small set of features** drives every anomaly — Japanese's lead, the
unrelated null's right tail, and the `[1.6]` skew candidates: the collapsed
length/gemination binaries (`cː`, `vː`), the tap `ɾ`, central vowels (`ə`, `ɨ`), and
marked fricatives (`ɸ`, `ž`, `ɟ`, `ḷ`, `ẓ`). Two structural causes:

1. **75% of the 728 columns are max-IDF singletons** (absent from all Austronesian
   reference languages) contributing nothing, so all signal concentrates in a handful
   of mid-IDF marked features.
2. **Modifications collapsed into single binaries** carry disproportionate, coarsely
   defined weight and are systematically under-weighted vs. per-segment PHOIBLE coding.

Spanish, by contrast, checks out: its loan phonemes appear only at the low minority
prevalence contact borrowing predicts, with the tap `/ɾ/` and velar `/x/` cleanly
localised to the 10 GMM-flagged influenced languages.
