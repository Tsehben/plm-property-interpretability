# Bounded P3 sequence-variant pilot — results

Run once at the frozen sample and replicate counts. P4 was not read; P0 remains
reserved and was never read. No sequence left this machine: WoLF PSORT ran in
the local container with `--network none`. Nothing committed, nothing pushed,
no earlier output altered.

| | |
|---|---|
| protocol hash | `2cffcda67aa857dd68dab44e248f335d3f2189ea924aba3257f3be2d6a6b16f2` |
| sample hash | `c121ffa56d2d507101b861cc2e6e87c877676368681454eb3ee3a2a24d3a1ce1` |
| frozen at | 2026-11-22T08:56:05−0500, **before any variant was generated** |
| sample | 40 Metazoa P3 proteins, **20 Extracellular / 20 Cytoplasm**, from a pool of 1,343 (362 / 981 available; quotas asserted and met) |
| generation | 8,000 forward passes, 1,645 s (27.4 min) |
| scoring | 840 sequences, animal mode, 267 s (4.4 min) |
| candidate | feature 1446, selected on P1+P2 during discovery, before any P4 output existed |

Scripts: [generate_variants.py](generate_variants.py), [score_variants.py](score_variants.py),
[analyse_pilot.py](analyse_pilot.py). Outputs and raw predictor output in [v1/](v1/).

---

## 1. A correction found before interpreting anything: omitted ≠ zero

The brief asked to verify from the implementation that a class omitted from the
output truly means zero. **It does not, and the frozen readout's extraction rule
had to be corrected because of it.**

`runWolfPsortSummary` prints only a *selected subset* of the classes that have a
non-zero score. Measured over all 840 records in this run:

| | |
|---|---|
| records with ≥ 1 non-zero class omitted from the printed line | **840 of 840** |
| mean non-zero classes omitted per record | **7.60** |
| most omitted in one record | **15** |
| records where **`extr`** was non-zero but omitted | **119** |
| records where **`cyto`** was non-zero but omitted | **89** |

Reading the printed line and treating an absent class as zero would therefore
have mis-measured the margin on 208 of 840 records.

The fix keeps the readout's *definition* unchanged and corrects only the
*extraction*: scores for every class are recomputed from the k = 32 nearest
neighbours that `--print-neighbors` reports, through the shipped similarity
matrix `data/animal.wolfu`, using the formula the package documents:

```
score(c) = sum over the k neighbours n of  U[c, class(n)]
```

This is verified against the implementation, not substituted for it: **every one
of the 2,650 scores the program itself printed is asserted to equal the
recomputation, and all 2,650 match.** Neighbour counts equal k in all 840
records.

These remain **predictor scores** — similarity-weighted neighbour vote counts.
Not probabilities, not biological measurements.

---

## 2. Integrity — all checks pass

| check | result |
|---|---|
| forward passes | **8,000 of 8,000 expected**; 8,000 step records; accounting complete |
| arm A is an exact no-op | **max ‖Δh‖ = 0 exactly**, 0 latent coordinates changed |
| unintended top-K changes, any arm | **0** |
| latent coordinates changed in the random-direction arms D/E | **0** (raw residual injection; the sparse code is never touched) |
| Δh confined to the edited token | **8,000 / 8,000** |
| non-finite values | **0** |
| sequence length constant; mutations only at allowed positions | **yes, all 800 variants** |
| valid replicate cells | **800 of 800** (40 × 4 × 5) |
| sequences scored | **840 of 840**; 0 missing, 0 unparsable, preflight file empty, stderr empty |
| control norm matching, max relative error | **3.17e-07** |
| random direction vs candidate Δh direction | n = 1,495, mean cos **0.0016**, max \|cos\| **0.070** — near-orthogonal, as an isotropic draw in ℝ¹²⁸⁰ should be |

No variant was regenerated. No setting was changed at any point.

---

## 3. Reference overlap — 5 of 40, retained

Exact frozen-sequence SHA-256 match against the 17,242 distinct sequences shipped
with WoLF PSORT:

| accession | class | matched reference | reference label | model |
|---|---|---|---|---|
| P61366 | Extracellular | OSTN_HUMAN | extr | animal |
| P09535 | Extracellular | IGF2_MOUSE | extr | animal |
| Q9NPH9 | Extracellular | IL26_HUMAN | extr | animal |
| P11672 | Extracellular | NGAL_MOUSE | extr | animal |
| Q08257 | Cytoplasm | QOR_HUMAN | cyto | animal |

**5 of 40 proteins (12.5%) are members of the predictor's own reference set**, all
in the animal model, all with reference labels agreeing with our frozen class.
Each is partly its own nearest neighbour.

The frozen sample is retained; no match was removed. Overlap-stratified contrasts
are reported in §5 **descriptively only** — with n = 5 against n = 35 they cannot
separate anything.

**Exact non-matching is not evidence of homology independence.** The other 35
proteins may still have close homologues in the reference set; nothing here
measures that.

---

## 4. Each arm's paired original-to-variant D

Protein-level values: replicates averaged within protein first, then summarised
across the 40 proteins. `D = margin(variant) − margin(original)`.

| arm | mean D | median D | sd | +/−/tied |
|---|---|---|---|---|
| **A** unmodified sampling | **+0.934** | +0.563 | 4.80 | 21 / 16 / 3 |
| **B** 1446 ×2 | +0.972 | +0.625 | 4.89 | 22 / 14 / 4 |
| **C** 1446 ×0 | +0.972 | +0.438 | 4.94 | 22 / 15 / 3 |
| **D** random +u | +0.872 | +0.188 | 4.82 | 20 / 17 / 3 |
| **E** random −u | +0.866 | +0.500 | 4.97 | 21 / 15 / 4 |

**The dominant effect is masked resampling itself, not the intervention.**
Arm A — which applies no edit at all — moves the margin by +0.93 on average.
All five arms sit within 0.11 of each other.

For scale: the original sequences' margins have a median of **+23.0** for
Extracellular and **−7.25** for Cytoplasm, a between-class gap of about 30
predictor-score units. Everything measured in §5 is on the order of 0.1 of that.

At replicate level (n = 800): mean +0.923, **median exactly 0**, sd 6.07, range
−35 to +23, and **173 of 800 (21.6%) are exact ties**.

---

## 5. Primary paired comparisons

Contrasts are formed **within a replicate** — the arms share that replicate's
uniforms — then averaged over replicates within protein, then summarised across
proteins. n = 40 proteins, 160 valid replicate pairs, for every contrast.

| contrast | mean | **median of paired differences** | sd | + | − | **exact ties** |
|---|---|---|---|---|---|---|
| **B − A** (1446 ×2 vs unmodified) | +0.0375 | **0.0000** | 0.636 | 9 | 3 | **28 / 40** |
| **C − A** (1446 ×0 vs unmodified) | +0.0375 | **0.0000** | 0.749 | 8 | 6 | **26 / 40** |
| **B − D** (1446 ×2 vs random +u) | +0.1000 | **0.0000** | 0.740 | 9 | 4 | **27 / 40** |
| **C − E** (1446 ×0 vs random −u) | +0.1063 | **0.0000** | 0.667 | 5 | 5 | **30 / 40** |

**Every contrast has a median paired difference of exactly zero**, with 26–30 of
40 proteins exactly tied. The means are 0.04–0.11 against protein-level standard
deviations of 0.64–0.75 and against an arm-A sampling effect of +0.93.

### The "typical protein" correction, concretely

The median *of the paired differences* is the right statistic. The difference
*between separately computed medians* is not, and here the two disagree badly:

| contrast | median of paired differences | difference of separate medians (**do not use**) |
|---|---|---|
| B − A | **0.0000** | +0.0625 |
| C − A | **0.0000** | −0.1250 |
| **B − D** | **0.0000** | **+0.4375** |
| C − E | **0.0000** | −0.0625 |

Using the wrong statistic would have reported B − D as **+0.44**, more than four
times its own mean and in flat contradiction to the fact that 27 of 40 proteins
show exactly no difference. It would also have flipped the sign of C − A. Both
columns are written to [v1/contrasts.tsv](v1/contrasts.tsv) side by side so the
distinction stays checkable.

### By class and by reference overlap — descriptive

| contrast | Extracellular (n=20) | Cytoplasm (n=20) | ref-match (n=5) | no match (n=35) |
|---|---|---|---|---|
| B − A | +0.063, med 0, 9 tied | +0.013, med 0, **19 tied** | 0.000, med 0 | +0.043, med 0 |
| C − A | +0.050, med 0, 7 tied | +0.025, med 0, **19 tied** | +0.150, med 0 | +0.021, med 0 |
| B − D | **+0.200**, med 0, 7 tied | **0.000, med 0, all 20 tied** | −0.025, med 0 | +0.118, med 0 |
| C − E | +0.188, med 0, 11 tied | +0.025, med 0, 19 tied | +0.150, med 0 | +0.100, med 0 |

What little movement exists sits in the Extracellular proteins, which is the
direction 1446's discovery association points — but every class-specific median
is also exactly zero, and the Cytoplasm cell for B − D is 20 of 20 tied. The
overlap strata have n = 5 against n = 35 and separate nothing.

---

## 6. Realised sequences

| arm | mean mutations (of 10 steps) | median | mean identity | identical to original | mean Δ KD whole | mean Δ KD N-term 30 |
|---|---|---|---|---|---|---|
| A | 7.244 | 8.5 | 0.9678 | 4.38% | −0.00486 | −0.0391 |
| B | 7.231 | 8.5 | 0.9680 | 4.38% | −0.00421 | −0.0334 |
| C | 7.250 | 8.5 | 0.9679 | 4.38% | −0.00515 | −0.0416 |
| D | 7.250 | 8.5 | 0.9678 | 4.38% | −0.00560 | −0.0451 |
| E | 7.231 | 8.5 | 0.9680 | 4.38% | −0.00415 | −0.0348 |

Ten mutations were never forced; sampling the original residue was permitted and
happened at about 2.75 of the 10 steps on average, and 4.38% of variants came
back **identical to their original**.

**These mutation counts were not matched by design and are not claimed to be
matched.** The design fixes the edit *budget* — ten masked positions — not the
realised number of substitutions. That they landed within 0.02 of each other
across arms is an observation about this run, not a property the design
guarantees.

Hydrophobicity is a biochemical readout only. The Δ KD values are small and run
in the same direction in all five arms, including the unedited one.

## 7. Feature activity and perturbation norms

| arm | steps | 1446 active | % | mean active steps per 10-step trajectory | mean ‖Δh‖ over non-zero | median | max |
|---|---|---|---|---|---|---|---|
| A | 1,600 | 747 | 46.69 | 4.67 | — (no perturbation) | — | 0 |
| B | 1,600 | 748 | 46.75 | 4.68 | 8.138 | 7.370 | 23.569 |
| C | 1,600 | 750 | 46.88 | 4.69 | 8.150 | 7.373 | 23.569 |
| D | 1,600 | 748 | 46.75 | 4.68 | 8.161 | 7.382 | 23.569 |
| E | 1,600 | 747 | 46.69 | 4.67 | 8.166 | 7.390 | 23.569 |

1446 is active at roughly 47% of masked positions, so about half of every
trajectory's steps carry a perturbation and the rest are exact no-ops by the
frozen rule. D and E carry slightly larger mean norms than B and C because each
control arm's norm is computed on **its own** current input, and the trajectories
diverge. That is the matching rule working as specified, not a mismatch: the rule
matches the intervention on the current input, **not whole generation
trajectories**.

---

## 8. Limitations

> **Three interpretation statements in this section were corrected on 2026-11-22.**
> The results above are unchanged; items 2, 4 and 5 below overstate what the design
> shows. See [../LOCALISATION_EVIDENCE_SUMMARY.md](../LOCALISATION_EVIDENCE_SUMMARY.md)
> §2 for the corrected wording.


1. **Exploratory P3 pilot evidence.** No efficacy threshold, no confirmatory
   claim, no decision taken on any interval.
2. **The intervention effect is small against the sampling effect it sits in.**
   Arm A alone moves the margin by +0.93; the four contrasts move it by
   0.04–0.11; the between-class baseline gap is about 30.
3. **Zero medians with heavy tie structure.** Predictor scores move in steps of
   ½ and ⅓ vote, so 26–30 of 40 proteins show exactly no difference and 21.6% of
   replicate cells are exact ties. A mean of +0.1 on such a distribution is
   carried by a handful of proteins.
4. **The random-direction controls are weaker than the ones that were
   unavailable.** Feature 1446 had **0 of 5** eligible activation-matched SAE
   controls under the prespecified ±20% rule. A random isotropic direction
   matches no activation statistic and lies outside the sparse code; B − D and
   C − E answer "does perturbation magnitude alone explain it?", not "is 1446
   special among features that behave like it?" The second question is still
   unanswered.
5. **Common random numbers are a coupling, not exact pairing.** Arms share each
   (protein, replicate, step) uniform, but once their sequences diverge the
   conditional distributions being inverted differ. This is variance reduction
   only.
6. **12.5% of the sample is in the predictor's reference set**, and homology
   dependence among the other 35 is unmeasured. No effective sample size is
   claimed.
7. **Not claimed:** matched mutation counts, preservation of family or domain
   structure, foldability, actual relocation of any protein, or any independent
   biological judge. WoLF PSORT is a kNN predictor over UniProt-derived
   annotations; agreement or disagreement with it is agreement or disagreement
   between two sequence-based models.
8. **No compound-class sensitivity rule was defined in the frozen protocol
   beyond excluding compound scores from the extr and cyto terms**, so none has
   been introduced after seeing results. Compound outputs are reported
   descriptively: the top recomputed class is a compound in 25 of 800 variant
   records (3.1%).
9. **Licence, stated factually:** the package's `LICENSE` file says
   `LICENSE: GPLv3` (updated 2023-04-04) while also retaining a 2005
   "Copyright: All Rights Reserved" header and an academic/commercial clause that
   permits academic on-site use and academic redistribution but excludes
   providing a prediction service over the internet. That is what the file says.
   **No legal confirmation has been obtained and none is claimed.** Nothing was
   redistributed and nothing was published in this task.

---

## 9. Files

| path | contents |
|---|---|
| [v1/frozen_protocol.json](v1/frozen_protocol.json) | protocol, sample, direction seeds, input hashes — written before generation |
| [v1/pilot_sample.tsv](v1/pilot_sample.tsv) | the 40 proteins with class, length, sequence hash |
| [v1/reference_overlap.tsv](v1/reference_overlap.tsv) | per-protein exact-match result against the reference sequences |
| [v1/per_protein/](v1/per_protein/) | original and final sequences, every variant, and all 8,000 step records (sampled residue, uniform, probabilities, candidate activation, edit norms, seeds) |
| [v1/wolfpsort_queries.fasta](v1/wolfpsort_queries.fasta) | the 840 queries |
| [v1/wolfpsort_raw/](v1/wolfpsort_raw/) | raw `summary.txt` (with neighbours), `summary.err`, `preflight.txt`, provenance |
| [v1/wolfpsort_scores.tsv](v1/wolfpsort_scores.tsv) | per query: status, full recomputed scores, margin, printed line, compound flags |
| [v1/arm_level_D.tsv](v1/arm_level_D.tsv), [v1/contrasts.tsv](v1/contrasts.tsv) | §4 and §5, with per-class and overlap strata |
| [v1/per_protein_D.tsv](v1/per_protein_D.tsv), [v1/per_protein_contrasts.tsv](v1/per_protein_contrasts.tsv) | protein-level values |
| [v1/sequence_descriptives.tsv](v1/sequence_descriptives.tsv), [v1/activity_and_norms.tsv](v1/activity_and_norms.tsv) | §6 and §7 |
| [v1/analysis_summary.json](v1/analysis_summary.json) | integrity block and output hashes |
