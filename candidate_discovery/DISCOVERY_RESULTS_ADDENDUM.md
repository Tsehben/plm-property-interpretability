# Addendum to DISCOVERY_RESULTS.md — interpretation corrections and a redundancy analysis

`DISCOVERY_RESULTS.md` and everything in `v1/` are preserved unchanged. The hashed candidate
list is unchanged and was re-verified at the end of this analysis: selected-list SHA-256
`070d0f43a9aea7a4de9102e377d7eb70b3ca7e44578ff709fddd3beecddbb917`. This addendum
characterises the existing 20 candidates; it does not replace them.

New artefacts, all from saved data — no new ESM extraction, no screening, no interventions,
no P3/P4, no network:

| artefact | path |
|---|---|
| analysis script | `candidate_discovery/redundancy_analysis.py` |
| pairwise behaviour, decoder similarity, overlap (190 pairs) | `v2_redundancy/candidate_pairwise.tsv` |
| stratified signed AUROCs, P1 and P2 separately | `v2_redundancy/stratified_auroc.tsv` |
| four-cell annotation distributions | `v2_redundancy/signal_annotation_cells.tsv` |
| interval-exception bookkeeping | `v2_redundancy/interval_exceptions.tsv` |
| summary, definitions, hashes | `v2_redundancy/redundancy_summary.json` |

---

## 1. Corrections to the original interpretation

**1.1 Truncation can affect the first-30 representations.** `DISCOVERY_RESULTS.md` §1 says
"Truncation does not touch the primary summary… `mean_nterm` is unaffected for every
protein." That is wrong as a statement about representations. It is true only of *positions*:
the first-1,022-residue policy removes C-terminal residues, so the 30 positions we pool over
are the same residues either way. But ESM-2 is bidirectional, so the layer-8 representation
*at* each of those positions is computed with attention over the whole input. Removing 500
C-terminal residues changes the context and therefore changes the first-30 activations.
553 discovery proteins are truncated (50 Extracellular, 503 Cytoplasm), and the truncation
rate is strongly class-imbalanced, so this is a differential, not a uniform, perturbation.
Corrected statement: *the primary summary pools a fixed set of positions that truncation
never removes, but the representations at those positions are context-dependent and are
affected by truncation.* The length-restricted stratification in §3 below is the direct
check: AUROCs on the 4,637 non-truncated proteins differ from the full-partition values by at
most 0.0051 for any candidate in either partition.

**1.2 Experimental-evidence agreement does not remove shared bias.** The original report
says the experimental subset "agrees, which weakens but does not remove the shared-bias
concern" — correct in direction, but the 240-protein subset is not an independent reference.
It shares the same literature-curation pipeline, the same organisms, the same Swiss-Prot
sampling and the same biology as the computationally evidenced majority. Proteins with
experimentally characterised signal peptides are preferentially well-studied, abundant,
secreted eukaryotic proteins. Agreement across evidence tiers rules out *one* specific
artefact — that the pattern exists only where an unnamed predictor made a call — and rules
out nothing else.

**1.3 Length correlation is not evidence of mediation.** The original §4 limitation says
feature 3900's "class association is partly length-mediated". That claim was not tested and
is not supported by a correlation. A correlation between activation and length, alongside a
correlation between length and class, is consistent with mediation, with confounding, with
reverse ordering, and with coincidence. The stratification in §3 is the relevant evidence,
and it does not support mediation: feature 3900's signed AUROC is 0.035 / 0.048 / 0.017 /
0.065 across the four P1 length quartiles and 0.011 / 0.030 / 0.023 / 0.031 across the P2
quartiles, against 0.037 and 0.015 unstratified. The association is essentially unchanged
within strata. Corrected statement: *3900 is length-correlated (ρ = +0.599), and length is
itself class-associated (AUROC 0.239), but conditioning on length does not attenuate the
feature's association, so length does not account for it. Nothing here establishes a causal
path in either direction.*

**1.4 Homology dependence and effective sample size are unmeasured.** No clustering exists
for this cohort. Every AUROC, every correlation and every count in the original report and
in this addendum treats 5,190 proteins as 5,190 observations. They are not: homologous
proteins share both sequence and label. The effective independent sample size is unknown and
smaller, by an unknown factor, than the nominal n. This affects precision, not point
estimates, and it is why no interval or significance statement is made anywhere.

**1.5 Signal-peptide association does not establish biological identity.** Saying a feature
"is a signal-peptide feature" is a claim about what it encodes. What we have is that its
activation co-occurs with annotated signal-peptide intervals in this cohort. The feature may
respond to hydrophobicity, to N-terminal charge, to a positional prior, to amino-acid
composition shared with transmembrane helices and transit peptides, or to a correlate of all
of these. Co-occurrence with an annotation is compatible with every one of those and
distinguishes none of them.

---

## 2. Redundancy: similar behaviour, different directions

Three quantities, kept conceptually separate because they answer different questions.

| measure | definition |
|---|---|
| **Spearman** | correlation of per-protein `mean_nterm` across the stated population — *behaviour* |
| **Decoder cosine** | cosine between SAE decoder rows `w_dec[f]` — *intervention direction in residual space*. Rows are not unit-norm in this checkpoint; norms 0.469–0.625 are recorded per feature |
| **Jaccard** | \|{proteins where A is active in the first 30} ∩ {same for B}\| / \|union\|, over the 5,190 discovery proteins — *co-firing* |

190 pairs:

| measure | median | min | max | notable |
|---|---|---|---|---|
| Spearman, all proteins | **0.697** | −0.725 | 0.833 | 91 pairs > 0.7; **0 pairs > 0.9** |
| Spearman, within Extracellular | **0.139** | −0.562 | 0.667 | **0 pairs > 0.7** |
| Spearman, within Cytoplasm | **0.179** | −0.293 | 0.464 | **0 pairs > 0.7** |
| Decoder cosine | **0.008** | −0.098 | **0.305** | 0 pairs > 0.5; 1 pair > 0.3 |
| Jaccard (active-protein sets) | **0.536** | 0.233 | 0.714 | — |

**Answer to the question asked: similar behaviour, different directions — not both, not
neither.** The candidates co-fire heavily (median Jaccard 0.54) and correlate strongly across
the whole sample (median 0.70), but that correlation is almost entirely produced by the class
split: conditioning on class collapses it to 0.14 and 0.18, with no pair above 0.7 in either
class. So they agree about *which proteins* are secreted-looking and disagree about
everything finer. Their decoder directions are nearly orthogonal (median cosine 0.008,
maximum 0.305 over 190 pairs).

Consequences, stated plainly:

- **These are not 20 independent mechanisms.** As detectors they are largely interchangeable,
  and their discovery statistics (AUROC 0.951–0.975 for the 19 Extracellular-directed ones)
  are near-duplicates of one another. Treating 20 candidates as 20 independent pieces of
  evidence would be double-counting.
- **They are, however, 20 distinct intervention directions.** Near-orthogonal decoders mean
  that editing feature 642 and editing feature 1446 inject nearly unrelated vectors into the
  residual stream, even though the two fire on nearly the same proteins. An intervention
  result on one therefore does *not* transfer to another, and "we tested a signal-peptide
  feature" is not a claim that generalises across the set.
- Both statements are discovery-side observations about dictionary geometry and co-firing.
  Neither predicts an intervention outcome.

---

## 3. Associations under simple controls

Signed AUROCs, P1 and P2 reported separately and never pooled. Length-quartile boundaries
come from the 5,190 discovery proteins only: **Q1 ≤ 223, Q2 223–407, Q3 407–672, Q4 > 672**
residues. Minimum to score a cell: 10 per class. **Every cell met it** — no cell was sparse
or single-class, so none is marked unsuitable; exact class counts are in the table.

Cell populations (identical for every feature):

| partition | stratum | n | Extracellular | Cytoplasm |
|---|---|---|---|---|
| P1 | all | 2484 | 519 | 1965 |
| P1 | original length ≤ 1022 | 2207 | 494 | 1713 |
| P1 | Metazoa / Fungi / Viridiplantae / Other | 1515 / 620 / 268 / 81 | 317 / 112 / 35 / 55 | 1198 / 508 / 233 / 26 |
| P1 | Q1 / Q2 / Q3 / Q4 | 513 / 652 / 696 / 623 | 205 / 125 / 125 / 64 | 308 / 527 / 571 / 559 |
| P2 | all | 2706 | 817 | 1889 |
| P2 | original length ≤ 1022 | 2430 | 792 | 1638 |
| P2 | Metazoa / Fungi / Viridiplantae / Other | 1792 / 546 / 288 / 80 | 635 / 92 / 38 / 52 | 1157 / 454 / 250 / 28 |
| P2 | Q1 / Q2 / Q3 / Q4 | 791 / 641 / 599 / 675 | 509 / 152 / 106 / 50 | 282 / 489 / 493 / 625 |

Representative rows (full matrix: `v2_redundancy/stratified_auroc.tsv`):

| feature | P | all | ≤1022 | Metazoa | Fungi | Viridipl. | Other | Q1 | Q2 | Q3 | Q4 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 1446 | P1 | 0.965 | 0.965 | 0.958 | 0.984 | 0.948 | 0.980 | 0.962 | 0.963 | 0.982 | 0.962 |
| 1446 | P2 | 0.979 | 0.978 | 0.979 | 0.978 | 0.964 | 0.973 | 0.987 | 0.966 | 0.983 | 0.975 |
| 642 | P1 | 0.963 | 0.966 | 0.958 | 0.976 | 0.951 | 0.982 | 0.967 | 0.963 | 0.985 | 0.905 |
| 642 | P2 | 0.981 | 0.981 | 0.982 | 0.985 | 0.923 | 0.991 | 0.994 | 0.960 | 0.967 | 0.977 |
| 2564 | P1 | 0.945 | 0.948 | 0.932 | 0.978 | 0.929 | 0.950 | 0.946 | 0.934 | 0.974 | 0.899 |
| 2564 | P2 | 0.964 | 0.964 | 0.964 | 0.986 | 0.904 | 0.950 | 0.979 | 0.935 | 0.964 | 0.969 |
| **3900** | P1 | 0.037 | 0.035 | 0.043 | 0.025 | 0.041 | 0.009 | 0.035 | 0.048 | 0.017 | 0.065 |
| **3900** | P2 | 0.015 | 0.015 | 0.013 | 0.012 | 0.064 | 0.023 | 0.011 | 0.030 | 0.023 | 0.031 |

Every candidate's association survives every stratification, in both partitions. The
weakest cells are Viridiplantae (0.885–0.975, n = 35 and 38 Extracellular) and P1 Q4
(0.878–0.962, n = 64 Extracellular) — still far from chance.

**This is not complete confounding removal.** Each stratification conditions on one measured
variable at a time; within-stratum variation in that variable remains, kingdom and length are
not conditioned on jointly, unmeasured confounders are untouched, and homology dependence
(§1.4) is present inside every cell. **No subgroup AUROCs were pooled**: there is no weighted
average in the table, because any such figure would need a stated weighting and target
population, and neither has been chosen.

---

## 4. Signal-annotation counterexamples — the four cells

Verified annotation-conditioned discovery sample, n = 5,009. `not_annotated` means a record
was retrieved carrying no signal-peptide annotation; **it is not confirmed absence** and is
never pooled with `unverifiable`.

| cell | n |
|---|---|
| Extracellular, annotated_present | 1,252 |
| Cytoplasm, annotated_present | **59** |
| Extracellular, not_annotated | **53** |
| Cytoplasm, not_annotated | 3,645 |
| Extracellular, annotated_present, experimental evidence | 236 |
| Cytoplasm, annotated_present, experimental evidence | **4** — too few to characterise; reported, not interpreted |

Median activation / fraction of proteins active:

| feature | E, annotated | **C, annotated** | **E, not_annot** | C, not_annot | E, annot, experimental |
|---|---|---|---|---|---|
| 1446 | 4.37 / 1.00 | **4.42 / 0.98** | **0.03 / 0.51** | 0.00 / 0.21 | 4.25 / 1.00 |
| 243 | 4.34 / 0.99 | **4.46 / 0.98** | 0.04 / 0.58 | 0.00 / 0.34 | 4.63 / 1.00 |
| 642 | 2.67 / 0.99 | 2.28 / 0.98 | 0.05 / 0.58 | 0.00 / 0.34 | 2.76 / 1.00 |
| 900 | 2.27 / 0.99 | **3.05 / 0.98** | 0.05 / 0.57 | 0.00 / 0.31 | 2.38 / 1.00 |
| 2564 | 1.74 / 0.98 | 1.45 / 0.88 | 0.00 / 0.47 | 0.00 / 0.26 | 1.81 / 0.98 |
| 3900 | 0.27 / 0.80 | 0.84 / 0.97 | 3.52 / 1.00 | 4.09 / 1.00 | 0.22 / 0.76 |

**What these comparisons show.** For all 19 Extracellular-directed candidates the activation
splits on **annotation status, not on localisation class**. Cytoplasm-labelled proteins that
carry a signal annotation activate at Extracellular levels (1446: median 4.42 against 4.37;
243: 4.46 against 4.34; 900 is *higher* in the Cytoplasm cell). Extracellular proteins
without an annotation sit at the not-annotated baseline (1446: median 0.03, active in 51%).
Feature 3900 shows the mirror image. This is exactly the pattern a signal-peptide-associated
feature would produce, and it is **incompatible** with a feature that tracks the localisation
label itself.

**What they cannot distinguish.** (a) Whether the feature responds to the signal peptide as a
functional element or to correlated N-terminal biochemistry — hydrophobicity, charge,
composition — shared with transmembrane segments, transit peptides and untargeted hydrophobic
N-termini. (b) Whether the 53 Extracellular `not_annotated` proteins are genuinely
signal-peptide-independent secretion or simply unannotated; `not_annotated` cannot separate
those, and the 51% that *are* active in that cell are consistent with either. (c) Anything
causal: these are co-occurrence distributions on discovery data. (d) The experimental-evidence
split is informative on the Extracellular side (n = 236) and uninformative on the Cytoplasm
side (n = 4).

---

## 5. The interval exceptions

**First, a correction.** `DISCOVERY_RESULTS.md` §3.2 says "Eighteen of twenty candidates are
strongly enriched inside the annotated interval". The count is wrong: **19 of 20 are enriched
on both criteria** — positive mean per-protein difference *and* a majority of proteins with
inside > outside. There is **one** exception, not two.

| criterion | features failing it |
|---|---|
| mean per-protein difference > 0 | **2564** only |
| inside > outside in a majority of the 1,306 proteins | **2564** only |

The nearest thing to a second exception is **3900**, the weakest enrichment in the set:

| feature | mean inside | mean outside | mean difference | inside > outside | fraction | zero inside |
|---|---|---|---|---|---|---|
| **2564** | 0.190 | 1.270 | **−1.079** | 166 / 1,306 | **0.127** | 547 |
| **3900** | 0.495 | 0.046 | +0.449 | 903 / 1,306 | 0.691 | 382 |
| 1052 | 0.516 | 0.207 | +0.308 | 1,041 / 1,306 | 0.797 | 104 |
| (for contrast) 1446 | 5.587 | 0.123 | +5.464 | 1,305 / 1,306 | 0.999 | 1 |

**The comparison uses the full visible sequence, not matched-position windows.** "Outside" is
every visible residue not in [start, end]. The two regions therefore differ in length (a
~20–25 residue interval against a median of roughly 180–470 remaining residues), in position
(N-terminal against everything downstream) and in sequence context. Several ordinary
explanations for outside-enrichment are left open by this design and are **not** excluded:

- A feature could fire on the mature N-terminus immediately after the cleavage site, which
  falls in "outside" but is still signal-peptide-adjacent.
- A longer "outside" region gives more opportunities for sparse firing, which inflates its
  mean for a feature with a diffuse profile.
- A feature could be position-sensitive rather than content-sensitive.

So **outside-enrichment is not evidence of a distinct mechanism.** Settling it needs a
matched-position window comparison — equal-length windows at matched offsets inside and
outside the interval — which the saved artefacts cannot supply, because per-residue
activations were streamed and deliberately not retained. That is a bounded recompute over
1,306 proteins (~3 minutes), not performed here.

---

## 6. Results table

| feature | rank | AUROC (P1 / P2) | direction | prevalence | interval enrichment | behaviour redundancy | decoder direction | notes |
|---|---|---|---|---|---|---|---|---|
| 1446 | 3 | 0.965 / 0.979 | Extracellular | 0.420 | **+5.46**, 1305/1306 | high (co-fires with the other 18) | near-orthogonal to all | cleanest interval-localised candidate |
| 243 | 4 | 0.959 / 0.979 | Extracellular | 0.509 | +6.15, 1301/1306 | high | near-orthogonal | largest inside mean |
| 642 | 2 | 0.963 / 0.981 | Extracellular | 0.512 | +1.85, 1289/1306 | high | near-orthogonal | highest P2 AUROC |
| 2564 | 14 | 0.945 / 0.964 | Extracellular | 0.444 | **−1.08**, 166/1306 | high as a detector | near-orthogonal | the single interval exception |
| 3900 | 1 | 0.037 / 0.015 | **Cytoplasm** | 0.951 | +0.45, 903/1306 | inverse-correlated with the rest | near-orthogonal | ρ(length) = +0.599; association survives length strata |
| other 15 | 5–20 | 0.937–0.956 / 0.954–0.980 | Extracellular | 0.28–0.62 | +0.31 to +3.81, ≥1041/1306 | high | near-orthogonal | near-duplicate discovery statistics |

---

## 7. Recommendation for the first intervention pilot — at most three existing candidates

Discovery data only. No strengths chosen, no controls drawn, nothing run.

**Recommended: 1446, 2564, 3900.**

1. **1446 — the representative signal-associated candidate.** AUROC 0.974 overall (P1 0.965,
   P2 0.979), top-ranked on P1 and 5th on P2, the cleanest interval localisation in the set
   (inside > outside in 1,305 of 1,306 proteins, only 1 protein with zero activation inside),
   and the four-cell table shows it tracks annotation rather than label. It is chosen as a
   *representative* of a group of near-duplicates, not as the unique best: 243, 642, 963 and
   2040 would serve nearly identically, and nothing in the discovery data distinguishes them
   meaningfully. 1446 wins on interval cleanliness and on having the largest decoder norm
   (0.625) in the set, which keeps the edit magnitude interpretable.

2. **2564 — the candidate with meaningfully different observed behaviour.** It reaches
   AUROC 0.957 while being the only candidate *not* enriched inside the annotated interval
   (inside 0.190 against outside 1.270, outside greater in 1,137 of 1,306 proteins). Whatever
   it responds to is positioned differently from the rest. **Support for "different" is
   qualified**: §5 shows the inside/outside comparison is unmatched for region length and
   position, so the difference is established as a *difference in where activation sits
   relative to the annotated interval*, not as a different mechanism. Including it tests
   whether the intervention effects we find are specific to interval-localised features or
   general to any high-AUROC N-terminal feature — which is the more informative question than
   running two interchangeable interval features.

3. **3900 — the strongest cytoplasmic-directed candidate, with its length association
   explicitly considered.** Highest |AUROC − 0.5| in the entire eligible set (0.477; AUROC
   0.023), and it is the only direction-reversed candidate in the top 20. Its length
   association is real — ρ = +0.599 with original length, and length alone reaches AUROC
   0.239 against this label — but conditioning on length does not attenuate it: within-quartile
   AUROCs are 0.035 / 0.048 / 0.017 / 0.065 on P1 and 0.011 / 0.030 / 0.023 / 0.031 on P2,
   against 0.037 and 0.015 unstratified. So length does not account for the association,
   though nothing establishes what does. Two cautions carry into the pilot design: it is
   near-ubiquitous (active in 95.1% of discovery proteins, mean 23.7 of 30 positions), so
   ablating it is a much broader edit than ablating 1446; and 382 of 1,306 annotated proteins
   show zero activation inside the interval, so it is a diffuse feature, not a focal one.

**Redundancy verdict, stated rather than worked around.** These three are not three
established mechanisms. 1446 and 2564 co-fire heavily (Jaccard in the 0.5 range, like most
pairs) and are both high-AUROC N-terminal features; 3900 is their inverse. What the three
*do* give are three genuinely different intervention directions — pairwise decoder cosines
near zero — plus one contrast that discovery data can already distinguish (interval-localised
against not). If the pilot has budget for one candidate only, use **1446**; if for two, add
**2564**, because the interval contrast is the only behavioural difference the discovery data
actually supports. 3900 is third because its breadth makes it the hardest to interpret, not
because it is the weakest association.
