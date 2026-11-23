# Sequence-variant pilot — specification (P3 only, **not executed**)

Prepared, not run. No variant has been generated and no project sequence has
been given to any predictor. This document fixes the design before any outcome
exists; it is the thing to freeze and hash at the start of the run.

Depends on: [../wolfpsort_setup/WOLFPSORT_READINESS.md](../wolfpsort_setup/WOLFPSORT_READINESS.md)
(installed and verified) and the frozen readout defined in its §5.

---

## 1. Candidate

**Primary: SAE feature 1446.** Chosen because discovery on P1+P2 selected it as a
representative signal-associated feature — rank 3, AUROC 0.974, higher in
Extracellular, active in the first 30 residues in 42.0% of discovery proteins
and at 7.16 of those positions on average when active, mean active magnitude
2.948. All of that was established **before any P4 output existed**.

**2564 is not promoted.** Its P4 paired contrast excluded zero, but selecting
the next study's candidate on that basis would choose a candidate using the
outcome of the evaluation it is supposed to be independent of, and would inherit
that contrast's coverage restriction (61 of 400 proteins, 55 of them
Extracellular, 50 of 61 cells at residue 20 — see
[../evaluation_p4/EVALUATION_ADDENDUM.md](../evaluation_p4/EVALUATION_ADDENDUM.md)).
1446 remains primary. 2564 and 3900 are not part of this pilot.

---

## 2. Sample — P3 only

* Partition **P3** (tuning). **P4 is not touched; P0 remains reserved.**
* **40 proteins, 20 Extracellular / 20 Cytoplasm**, deliberately balanced; this
  is not natural prevalence and no prevalence claim follows from it.
* **Kingdom restricted to Metazoa**, so that one WoLF PSORT organism model
  (`animal`) applies to every protein. Scores are not comparable across organism
  models (different `k`, disjoint reference sets, different class vocabularies),
  so per-protein model switching is not an option. Eligible P3 Metazoa pool:
  362 Extracellular, 981 Cytoplasm.
* Eligibility: original length **60–1022**, standard amino acids only, from the
  prepared analysis view; **engineering proteins and all `p0_base_collision`
  records excluded**.
* Selection rule: **ascending frozen-sequence SHA-256**, taking the first 20 per
  class. Deterministic, and **not conditioned on 1446 activation, signal
  annotation, length, or any outcome.** Quotas are asserted; if a quota cannot be
  met the run stops and reports it rather than changing the rule.

---

## 3. Fixed masking and edit budget

* **Schedule: 10 positions, `{3, 6, 9, 12, 15, 18, 21, 24, 27, 30}`**, one-based
  residue coordinates, applied in that order. Identical for every protein and
  every arm. Position 1 is excluded so the initiator methionine is never
  resampled.
* Positions lie inside the first 30 residues — where 1446 was screened and where
  it is active.
* **Edit budget: at most 10 substituted residues per variant**, one per schedule
  step. Substitutions only, so **length is never changed**.
* Generation, per schedule step: mask the position, run one forward pass under
  the arm's intervention, sample a residue from the AA-20 renormalised
  distribution at temperature 1.0, write it into the sequence, continue. Ten
  forward passes per variant.
* Token offsets are checked against the residue at each position before editing,
  as in the P4 run.

Because substitutions preserve length (≥ 60) and draw only from the 20 standard
amino acids, **no variant can trigger WoLF PSORT's silent-drop rules** (length
< 30, length > 10,000, residue `O`). Input/output reconciliation is still run on
every batch; see §7.

---

## 4. Arms — five, all on the same inputs

| arm | intervention at the masked token |
|---|---|
| **A** | **unmodified sampling control** — no edit at all; same schedule, same masks, same random numbers |
| **B** | 1446 **× 2** (amplify) |
| **C** | 1446 **× 0** (ablate) |
| **D** | **norm-matched random-direction control**, matched to arm B's realised ‖Δh‖ |
| **E** | **norm-matched random-direction control**, matched to arm C's realised ‖Δh‖ |

Edits use the validated direct decoder-delta implementation. Arm A is the
reference for "what the model does anyway"; arms D and E isolate the effect of
*perturbation magnitude* from the effect of *this particular direction*.

**Random-direction control construction.** At the edited token, draw
`g ~ N(0, I₁₂₈₀)` from a seed derived from (sample hash, protein, draw, step,
arm), and set `Δh_rand = (g/‖g‖₂) · ‖Δh_candidate‖₂`, using the ‖Δh‖ actually
realised by the paired candidate arm on that same forward pass. Relative norm
tolerance 1e-5, as in the P4 run. If the candidate is inactive at that token
there is no norm to match and the condition is recorded unavailable — **never
imputed**.

### These are not the controls that were unavailable for 1446

Stated plainly, because the substitution is the weakest point of the design:

Under the prespecified activation-matching rule (±20% on both discovery
activation frequency and mean active magnitude, no AUROC restriction), feature
**1446 had 0 eligible matched SAE controls out of 5 requested** — only 3
non-candidate eligible features lie within ±20% of its mean active magnitude of
2.948, and none of those also matches its N-terminal prevalence of 0.420. That
is why 1446 was descriptive-only in the P3 pilot and in the P4 evaluation.

A random isotropic direction is a **different and weaker control**:

* it is not an SAE decoder direction and almost surely lies outside the top-K
  support, so it does not respect the model's sparse code at all;
* it matches **no** activation statistic — not frequency, not magnitude, not
  which proteins or positions the feature fires on;
* it carries no sign or directional semantics, so there is no meaningful
  "ablate" versus "amplify" for it — arms D and E differ only in the magnitude
  they are matched to;
* it controls for exactly one thing: **the size of the residual-space
  perturbation at the edited token**.

On direction alone the two are comparable — a random vector is near-orthogonal
to the decoder row (|cos| ≈ 0.028 typically), and the real SAE controls were
near-orthogonal too (median cosine 0.008 among the 20 candidates). On everything
else they are not. **A contrast of B against D answers "does magnitude alone
explain it?", not "is 1446 special among features that behave like it?"** The
second question remains unanswerable for 1446 with the current feature set, and
this pilot does not answer it.

---

## 5. Pairing the stochastic sampling

**Common random numbers.** For each (protein *p*, draw *d*, schedule step *j*),
one uniform deviate `u = PRNG(sample_sha256, p, d, j)` is derived and used by
**all five arms** via inverse-CDF sampling over the AA-20 distribution with
tokens in a fixed canonical order. The seed does not depend on the arm, so at a
given step every arm consumes the same `u`.

* At **step 1** all five arms see an identical input sequence and an identical
  `u`. Any difference in the sampled residue there is attributable to the shifted
  distribution alone, with zero sampling noise.
* From step 2 onward the arms' sequences may have diverged, so the inputs differ
  and the pairing is no longer exact. Shared `u` still removes a large part of
  the sampling variance, but **this is variance reduction, not exact pairing**,
  and the degradation is a declared limitation rather than something the design
  removes.
* **4 independent draws** per protein per arm, with the draw index entering the
  seed.

Analysis therefore reports step-1 divergence separately from whole-variant
outcomes, since step 1 is the only exactly paired comparison.

---

## 6. What is recorded

Per (protein, arm, draw):

* **Realised mutations** — for each of the 10 schedule positions: original
  residue, sampled residue, changed flag; total Hamming distance from the
  original; a flag for variants that came back **identical** to the original
  (a real and informative outcome, since the model often resamples the original
  residue — it is not a failure and is not discarded).
* **Hydrophobicity** — Kyte–Doolittle mean over the whole realised variant and
  over residues 1–30, with Δ against the original. A biochemical readout, not a
  localisation prediction.
* **Intervention bookkeeping**, as in P4 — latent before/after, number of latent
  coordinates changed, unintended top-K evictions/insertions, ‖Δh‖, realised norm
  match error, confinement of Δh to the edited token, finiteness.
* **Predictor outputs** — raw `summary.txt`, `summary.err`, `preflight.txt` and
  `run.json` preserved per batch; parsed per-record site→score map; `E`, `C`,
  `Δ`; the compound-class flags; and an explicit status for **every input
  identifier**, reconciled against the input FASTA.

Full float precision throughout; configuration, sample and seeds hashed before
generation begins.

---

## 7. Readout and analysis

The readout is the one frozen in
[../wolfpsort_setup/WOLFPSORT_READINESS.md §5](../wolfpsort_setup/WOLFPSORT_READINESS.md):
`Δ(s) = score(extr) − score(cyto)` over single-site classes, and
`D(o, v) = Δ(v) − Δ(o)`.

* Organism model **`animal`**, fixed for every protein.
* Aggregation: average the 4 draws within arm within protein, then contrast arms
  **within protein** — B−A, C−A, B−D, C−E — then summarise across proteins. The
  **protein is the sampling unit**; draws, positions and arms are not independent
  units.
* A protein contributes to a contrast only if both arms produced a prediction for
  that draw; otherwise the cell is dropped, counted and reported — never imputed
  as zero. An absent `extr` or `cyto` class is a genuine zero and is kept.
* Pre-registered sensitivity analysis: repeat excluding proteins whose original's
  top class is a compound.
* Reported alongside: the no-change rate, the Hamming-distance distribution, and
  the step-1 exactly-paired divergence.

**This is a P3 tuning pilot.** It is exploratory. It produces no confirmatory
claim, and no decision is taken on whether an interval contains zero. Its purpose
is to establish whether the generation loop, the controls and the readout are
technically usable at all — and whether `D` moves enough to be distinguishable
from its 0.5-vote granularity.

---

## 8. Estimated cost

Measured rates, not guesses: ESM-2 650M on this machine ran 13,500 forward passes
in 2,592.6 s during the P4 evaluation (**0.192 s/pass**); WoLF PSORT `animal`
under emulation costs **11.3 s fixed + 0.303 s/sequence** (benchmarked at 100 and
1,040 packaged reference sequences: 41.6 s and 326.5 s).

| item | count | time |
|---|---|---|
| forward passes | 40 proteins × 5 arms × 4 draws × 10 steps = **8,000** | ≈ 1,540 s (**26 min**) |
| variants generated | 40 × 5 × 4 = **800** | — |
| sequences scored | 800 variants + 40 originals = **840**, one `animal` batch | ≈ 266 s (**4.5 min**) |
| parsing, reconciliation, analysis | — | < 1 min |
| **total wall clock** | | **≈ 35 min** |

Storage ≈ 60 MB (per-protein JSON at P4 density, plus raw predictor output).
No network, no external service, no cost beyond local compute.

Doubling the draws to 8 or the sample to 80 roughly doubles it, still under 1.5
hours — but the sample size is frozen at 40 before the run and **is not to be
expanded after seeing results**.

---

## 9. Open before execution

* The reference-overlap check from
  [../wolfpsort_setup/WOLFPSORT_READINESS.md §4](../wolfpsort_setup/WOLFPSORT_READINESS.md)
  has not been run. Exact-sequence overlap between the P3 sample and the 17,242
  WoLF PSORT reference sequences should be measured and reported **before** the
  pilot is interpreted, since an overlapping protein is partly its own neighbour.
  It must not touch P0.
* Licence confirmation (§1 of the readiness document) before any output is
  redistributed or published.
* Unmeasured homology dependence among the 40 proteins, as everywhere else in
  this project. No effective sample size is claimed.
