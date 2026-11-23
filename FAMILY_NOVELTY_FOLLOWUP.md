# Follow-up: layer-26 probe on P4 proteins with no previously recorded Pfam family

> **This is a FOLLOW-UP analysis on previously examined P4 data.** P4 is development data that
> has already carried the 400-protein intervention evaluation and the layer-sweep evaluation.
> It is not a fresh test set and no figure here is a held-out result.
>
> **"Unseen" means:** the protein's recorded, applicable **Pfam family-type** signatures do not
> appear in the probe's fitting partitions (P1+P2) or in the layer-selection partition (P3),
> under this annotation snapshot. It is **not** a claim of evolutionary independence.

Saved outputs only — no refitting, no extraction, no layer re-selection. **P0 untouched.**
Results: [`subgroup_results.tsv`](family_novelty_followup/subgroup_results.tsv) ·
[`subgroup_results.json`](family_novelty_followup/subgroup_results.json) ·
membership [`subgroup_membership.tsv`](family_novelty_followup/subgroup_membership.tsv)
(SHA-256 `dacd8a55…`) · [`subgroup_frozen.json`](family_novelty_followup/subgroup_frozen.json).

---

## 1. Membership frozen before any score was read

[`freeze_subgroup.py`](family_novelty_followup/freeze_subgroup.py) opened **no prediction
file**. It selected P4 proteins carrying at least one *applicable* Pfam family-type signature,
none of which occurs in P1+P2+P3, and verified every row against the frozen manifest and the
actual sequence.

| Check | Result |
|---|---|
| Members | **443** |
| **Extracellular / Cytoplasmic** | **183 / 260** — matches expectation |
| Accession present in the core manifest | 443/443 |
| Frozen SHA-256 matches the manifest | 443/443 |
| Frozen SHA-256 matches the actual sequence | 443/443 |
| Target matches the manifest | 443/443 |
| Partition = P4, role = evaluation | 443/443 |
| Verification failures | **0** |

Targets and 30 %-cluster assignments were independently re-checked against the saved
prediction file after loading; both agreed for all 443 rows.

## 2. Results

Protein-weighted, on the 443-protein subgroup. Both models are the **frozen** fits from the
sweep — the layer-26 ESM-only probe and the matched 41-feature amino-acid composition baseline.

| Model | n | EX | CY | Clusters (30 %) | AUROC | AP |
|---|---|---|---|---|---|---|
| **Layer-26 ESM-only** | 443 | 183 | 260 | 234 | **0.986500** | **0.985300** |
| **Composition baseline** | 443 | 183 | 260 | 234 | 0.975000 | 0.963800 |

**Paired difference, ESM-only − composition**, by the sweep's paired 30 %-identity cluster
bootstrap (whole clusters drawn with replacement, same resample indices for both models,
mixed-class clusters never split, unstratified, B = 10,000, seed 42, percentile 95 %):

| Metric | Difference | 95 % CI | Excludes 0 |
|---|---|---|---|
| **AUROC** | **+0.011517** | **[−0.005361, +0.037176]** | **no** |
| AP | +0.021484 | [+0.000276, +0.108820] | yes, marginally |

| Resampling diagnostics | |
|---|---|
| Clusters resampled | **234** |
| **Largest cluster's share** | **83 / 443 = 18.74 %** (`Q2VBP4`) |
| Undefined resamples excluded | **0** of 10,000 |

**The largest cluster deserves separate mention.** `Q2VBP4` — the snake-venom three-finger-toxin
cluster — contributes 83 proteins, **all extracellular**. That is 18.7 % of the subgroup but
**45.4 % of its 183 extracellular proteins**. The subgroup also has 181 singletons and 59.1 % of
its members in clusters of ≥ 2.

## 3. Context: the full-P4 difference

| | n | EX / CY | Prevalence | Clusters | ΔAUROC | ΔAP |
|---|---|---|---|---|---|---|
| **Full P4** | 2,366 | 571 / 1,795 | 0.2413 | 1,547 | **+0.013299** | +0.045458 |
| **This subgroup** | 443 | 183 / 260 | **0.4131** | 234 | **+0.011517** | +0.021484 |

The full-P4 figure recomputes as +0.013299 against the +0.013298 reported in the sweep; the
difference is 1 × 10⁻⁶ and comes from recomputing on the float32 saved predictions rather than
the float64 in-memory arrays.

⚠️ **AP is not comparable between these two rows.** Average precision depends on the base rate,
and the subgroup's prevalence (0.4131) is nearly double full P4's (0.2413). The subgroup's
higher absolute AP for both models reflects that shift, not better performance. **AUROC is the
metric to compare across the two rows.**

## 4. What this does and does not show

**The layer-26 probe clearly distinguishes the two classes in this subset.** AUROC **0.9865**
on proteins whose recorded families were all absent from the fitting and layer-selection
partitions. The composition baseline also separates them well, at 0.9750. Neither is near
chance; the ability to tell extracellular from cytoplasmic does not collapse when the recorded
family identifiers are novel.

**Its advantage over composition is not resolved here.** The AUROC difference is +0.0115 with a
95 % interval of **[−0.0054, +0.0372]**, which includes zero. On this subgroup the data are
compatible with no advantage and with an advantage roughly three times the full-P4 estimate.
The AP interval excludes zero, but its lower bound is **+0.000276** — close enough to zero that
it should not be leaned on, and it was computed without any multiplicity adjustment alongside
the AUROC interval.

**The two point estimates are close.** +0.0115 here against +0.0133 on full P4 — a gap of
0.0018, far inside the subgroup interval's width of 0.0426. What changed between the two rows
is mostly **precision**, not the estimate: 443 proteins and 234 clusters instead of 2,366 and
1,547.

> **A difference between this subset and full P4 cannot by itself establish that family novelty
> caused a change in performance.** The subgroup differs from full P4 on several things at
> once, none of them isolated here: it is 5.3× smaller with 6.6× fewer clusters; its class
> prevalence is nearly double; it excludes every protein lacking an applicable family-type
> signature — **1,633 of 2,366 P4 proteins (69 %)**, which is a large and non-random slice; and
> it is far more concentrated, with 45 % of its extracellular class in one toxin cluster
> against 16 % in full P4. Any of these could move a subgroup estimate on its own.

**Two further limits on "unseen".** It is defined on *recorded* Pfam family-type signatures
under one annotation snapshot: Pfam coverage is incomplete by construction, and two distinct
family accessions can sit in the same clan. And unseen-by-Pfam does not imply unseen-by-sequence
— these proteins were selected on family identifiers, not on homology, and some share 30 %
clusters with the fitting partitions.

**Bottom line.** On P4 proteins whose recorded families are novel to the probe's fitting and
selection partitions, the layer-26 ESM probe still separates extracellular from cytoplasmic
proteins (AUROC 0.9865), and its margin over amino-acid composition is in the same direction
and of similar size as on full P4 — but on this subgroup that margin is **not statistically
resolved**, and the comparison cannot attribute any part of the change to family novelty.
