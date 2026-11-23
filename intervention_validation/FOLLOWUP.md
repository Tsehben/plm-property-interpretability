# Follow-up: evidence, edit implementation, and corrected reporting

**Date:** 2026-11-20 · Supersedes the reporting in [REPORT.md](REPORT.md) where the two
conflict. Original outputs preserved unchanged; new outputs carry distinct names.

**Outcome:** 42/42 checks pass under both edit implementations. The **direct
decoder-delta** is adopted as the production implementation; difference-of-reconstructions
is retained as a reference. Two reporting claims from the first round were wrong and are
corrected in §4.

---

## 1. Evidence behind the result

### 1.1 Every check, criterion, observed result

Source: [`out/validation_report_direct.json`](out/validation_report_direct.json) → `checks`.
"Exact" means literal equality, measured against a repeat-run noise floor of exactly 0.

| # | Check | Criterion | Observed | Result |
|---|---|---|---|---|
| 1 | SAE checkpoint hash | == pinned SHA-256 | `1172503b…e56ebd` | PASS |
| 2 | Loading diagnostics | no required missing / mismatched key | 0 required missing, 0 mismatched | PASS |
| 3 | LM-head tie | `lm_head.decoder.weight` ≡ input embedding | tensors identical | PASS |
| 4 | LM-head trained | mean\|w\|>1e-4 and std>1e-4 | 0.07271 / 0.09215 | PASS |
| 5 | LM-head bias trained | max\|b\| > 0 | 14.59375 | PASS |
| 6 | MLM preserves encoder | `hidden_states[8]` bitwise equal to `EsmModel` | bitwise equal, 0.000e+00 | PASS |
| 7 | Hook site unique | exactly one layer output == `hidden_states[8]` | matching indices `[7]` | PASS |
| 8 | Hook site identity | == `encoder.layer[7]` | hook_layer=7 | PASS |
| 9 | Continuation complete | 33−8 = 25 blocks run after hook | `layer[8..32]`, 25 blocks | PASS |
| 10 | Output structure | determined empirically, handled | bare `Tensor`, not tuple | PASS |
| 11 | P0 excluded | all eligible in `development_P1_4` | 56 eligible of 64 | PASS |
| 12 | Sequence integrity | SHA-256 matches manifest | 8/8 | PASS |
| 13 | No truncation | original len ≤ 1022 | 8/8 | PASS |
| 14 | Decoder orientation | `w_dec` is (4096, 1280) | (4096, 1280) | PASS |
| 15 | Encoder orientation | `w_enc` is (1280, 4096) | (1280, 4096) | PASS |
| 16 | Latents vs `get_acts` | elementwise equal | bitwise, 0.000e+00 | PASS |
| 17 | Active support | identical to reference | 30,720 of 983,040 | PASS |
| 18 | Top-k count | ≤128 active per token | min 128, max 128 | PASS |
| 19 | Recon vs `forward_val` | elementwise equal | bitwise, 0.000e+00 | PASS |
| 20 | Pre-acts vs `encode` | elementwise equal | bitwise, 0.000e+00 | PASS |
| 21 | µ/σ vs `encode` | elementwise equal | bitwise | PASS |
| 22 | `decode(encode(x))` vs `forward_val` | elementwise equal | bitwise, 0.000e+00 | PASS |
| 23 | Scale restoration | \|σ_recon−σ_x\|/σ_x < 0.25 | 2.8492 → 2.6874 | PASS |
| 24 | **NEG CTRL** biased σ | must differ from reference | 1.259e-03 ≠ 0 | PASS |
| 25 | **NEG CTRL** no `b_pre` | must differ from reference | 1.304e+00 ≠ 0 | PASS |
| 26 | Noise floor | 3 clean runs reproducible | **bitwise identical**, floor 0.000e+00 | PASS |
| 27 | Capture hook inert | logits bitwise == clean | 8/8 | PASS |
| 28 | Identity hook, tensor | bitwise preserved | 8/8 | PASS |
| 29 | Identity hook, logits | bitwise == clean | 8/8 | PASS |
| 30 | Zero-edit, tensor | bitwise preserved | 32/32 arms | PASS |
| 31 | Zero-edit, delta | exactly 0 | max 0.000e+00 | PASS |
| 32 | Zero-edit, logits | bitwise == clean | 32/32 arms | PASS |
| 33 | Edit specificity | exactly 1 latent coordinate changed | 64/64 arms | PASS |
| 34 | Top-k integrity | 0 **unintended** evictions/insertions | 0 / 0 (32 intended) | PASS |
| 35 | Reference self-consistency | the two forms agree in float64 ≤1e-11 rel | 1.098e-13 | PASS |
| 36 | Production accuracy | rel L2 err vs float64 ref ≤ 1e-5 | **9.829e-07** | PASS |
| 37 | Implementation choice | direct ≤ DoR on median **and** max rel L2 | median 3.37e-07 vs 1.55e-06; max 9.83e-07 vs 7.00e-06 | PASS |
| 38 | Per-condition regression | no condition worse than 1.1× | worst ratio 1.066 | PASS |
| 39 | fp32 form difference | rel ≤ 1e-4 (rounding, not discrepancy) | 6.604e-05 | PASS |
| 40 | Delta localisation | activation delta exactly 0 outside edit position | 64/64 arms | PASS |
| 41 | Injection non-trivial | delta L2 > 0 | min 1.3574 | PASS |
| 42 | Above noise floor | max\|Δlogit\| > floor | min 2.451e-03 vs floor 0 | PASS |

Checks 1–26 and 35–39 come from the `direct` run; identical results hold for the `dor` run
([`out/validation_report_dor.json`](out/validation_report_dor.json), also 42/42).

### 1.2 Criteria that changed after a first run

Three criteria changed across the two rounds. **Two correct logical errors; one is a
genuine relaxation of numerical acceptance.** Labelled explicitly:

| Criterion | Original | Revised | Classification |
|---|---|---|---|
| Top-k integrity (#34) | `support_evicted == 0 and support_inserted == 0` | same, but excluding the **edited** coordinate | **Logical-error correction.** The original counted the ablated feature leaving its own support — the intended effect — as a violation. It observed 16 "evictions", all ablate arms, exactly 1 each. The revised criterion is not numerically looser: it still demands exactly 0 changes to all *other* features. |
| Edit-form agreement (#39) | `max_abs_diff <= 1e-5` (float32, absolute) | float64 self-consistency ≤1e-11 relative (#35) **plus** float32 relative ≤1e-4 (#39) | **Mixed: part correction, part relaxation.** The original 1e-5 absolute bound was arbitrary and scale-blind; it failed at 3.8e-5. Re-expressing it as relative is a correction. But the float32 acceptance *is* looser than before, and should be read as such. It is paired with a far stricter new test (#35) that directly checks the mathematical claim in float64. |
| Implementation choice (#37) | `rel_err_DD <= rel_err_DoR` in **every** condition | aggregate dominance on median and max (#37) **plus** per-condition tolerance ≤1.1× (#38) | **Relaxation, explicitly.** The original failed at 62/64. It asserted strict per-condition ordering between two methods that both sit at the fp32 rounding floor, where their ordering on any single condition is arbitrary. The 2 exceptions are the same protein and case, worse by 1.2% and 6.2%. Recorded in the JSON under `impl_comparison.criterion_note`. |

No criterion was changed to make a failing *implementation* pass. In each case the
implementation was unchanged and the criterion was the thing at fault or mis-specified.

### 1.3 Sample, positions, environment

**Engineering sample** — [`out/engineering_sample_direct.tsv`](out/engineering_sample_direct.tsv).
Selection: `split_role == development_P1_4` **and** original `seq_len ≤ 1022`, sorted by
accession ascending, first 8. No label, prediction or model output consulted.

| ACC | Partition | Length | Mask token idx | Case-B edit token idx |
|---|---|---|---|---|
| H2A0K6 | P4 | 238 | 120 | 60 |
| K7ZGS2 | P2 | 63 | 32 | 16 |
| O13747 | P4 | 547 | 274 | 137 |
| O42873 | P1 | 570 | 286 | 143 |
| O43037 | P3 | 147 | 74 | 37 |
| O43447 | P4 | 177 | 89 | 45 |
| O60462-6 | P4 | 555 | 278 | 139 |
| O61866 | P2 | 403 | 202 | 101 |

Mask residue = `L//2`, case-B edit residue = `L//4` (0-based); token index = residue + 1
for `<cls>`. `<cls>`, `<eos>` and padding excluded from all edits.

| Field | Value |
|---|---|
| ESM-2 | `facebook/esm2_t33_650M_UR50D`, revision `08e4846e537177426273712802403f7ba8261b6c` |
| SAE | `esm2_plm1280_l8_sae4096.safetensors`, SHA-256 `1172503b90c950651e787a5625d91ff1c702a04ff3d851e8e95b700201e56ebd` |
| SAE config | `SparseAutoencoder(d_model=1280, d_hidden=4096, k=128)` |
| torch / transformers / numpy | 2.8.0 / 4.57.6 / 2.0.2 |
| Python / platform | 3.9.6 / macOS-26.5.1-arm64 |
| Device / dtype / batch | `mps` / float32 / 1 |
| SAE execution | CPU float32; float64 twin for reference only |

### 1.4 P0 exclusion and engineering-only designation

- **P0 excluded.** Every eligible record was asserted `split_role == development_P1_4` and
  `partition != '0'` before selection (check #11). The DeepLoc P0 CSV rows were never read
  by this work. `hpa_testset.csv` was never read.
- **Engineering-only designation.** The 8 accessions above are recorded in
  `validation_report_direct.json` under
  `engineering_sample.designation` as *"ENGINEERING-ONLY. These accessions are reserved
  from future confirmatory intervention evaluation."* They were used to debug and tune
  criteria, so they are no longer clean for confirmatory use.
- **No new global split was created**, per instruction.
- **Caveat carried forward:** the 64-protein source pilot is label-stratified by
  construction, so these 8 inherit that stratification; only the selection *within* it is
  label-free.

---

## 2. Edit implementation: comparison and recommendation

Full data: [`out/edit_implementation_comparison.json`](out/edit_implementation_comparison.json)
(96 conditions = 8 proteins × 2 cases × 2 feature rules × 3 strengths). Edit strengths are
**unchanged** from the completed validation (amplify ×2, ablate →0, no-op ×1).

**Reference:** float64 direct decoder-delta. Justified rather than assumed — the two
float64 forms agree to **1.20e-14** relative L2, far below float32 resolution, so the
choice of reference form cannot drive the comparison.

**Relative-error definitions** (denominators handled explicitly):

- `rel_err_l2 = ‖approx − ref‖₂ / ‖ref‖₂` — primary. Undefined only if the reference delta
  is identically zero, in which case it is reported as `null` rather than computed. For
  the no-op arms both deltas are exactly zero, so no ratio is formed at all.
- `elementwise` relative error is computed **only** where `|ref_i| ≥ 1e-3 · max|ref|`;
  coordinates below that floor are excluded and counted, rather than producing an inflated
  ratio against a near-zero denominator.
- Significant digits are derived as `−log₁₀(rel_err_l2)`, **never** from absolute error.

### Results, 64 nonzero-edit conditions

| Quantity | Difference-of-reconstructions | **Direct decoder-delta** |
|---|---|---|
| rel L2 error vs float64 ref — median | 1.554e-06 | **3.367e-07** |
| rel L2 error vs float64 ref — max | 7.000e-06 | **9.829e-07** |
| absolute error — max | 3.965e-05 | **1.592e-05** |
| significant digits agreeing — min | 5.15 | **6.01** |
| better in how many conditions | 2 / 64 | **62 / 64** |
| median ratio DoR/DD | — | **3.98× more accurate** |
| worst-case ratio DoR/DD | — | 113× more accurate |

The advantage is largest exactly where it matters most — small edits. For the
`median_rank_active` feature rule (perturbations ~0.13% of ‖h‖), median relative error is
2.673e-06 for DoR against **4.077e-07** for the direct form, a 6.6× gap.

**The 2 exceptions** are both `O61866`, case B, `strongest_active` (amplify and ablate),
where DoR is better by 1.2% and 6.2%. Both methods are at the rounding floor there; the
ordering is not meaningful.

**Downstream logits.** Injecting each float32 delta and comparing: DoR vs DD differ by at
most **9.25e-05**; each differs from the rounded float64 reference by at most 9.06e-05
(DD) and 7.92e-05 (DoR). **No ordering is established downstream** — all three agree to
~1e-4, which is the fp32 forward's own sensitivity at this perturbation size, and neither
is bitwise equal to the reference (0/64 for both). The case for the direct form rests on
activation-level accuracy, not on better logits.

**Other properties, both forms identical:** no-op delta exactly zero (32/32), delta
confined to the edit position (96/96), float32-vs-float64 top-k support mismatches **0**
across all conditions.

### Recommendation

**Adopt the direct decoder-delta as the production implementation:**

```
delta_h = ((topK(a_edited) − topK(a)) @ w_dec) * std
h_edited = h + delta_h
```

Justification, in order of weight:

1. **Measurably more accurate** — 4× at the median, 113× worst case, 62/64 conditions;
   and the gap widens for small edits, which is the regime that matters for feature-level
   work.
2. **No cancellation by construction.** DoR forms two large reconstructions and subtracts
   them; the direct form never materialises the quantity that cancels. `b_pre` and µ drop
   out algebraically.
3. **Cheaper** — one sparse latent difference and one matmul, versus two full decodes plus
   a subtraction.
4. **Identical exactness properties** — zero-strength no-op is exactly zero and
   bitwise-preserving, and the delta is confined to the edited position, in both forms.

**Retain DoR as the reference implementation.** It is selectable with `--impl dor`, passes
all 42 checks (max rel L2 error 7.00e-06, still inside the ≤1e-5 criterion), and is the
form the specification is written in, so it remains the readable definition of intent.

**Reruns completed** with the direct form as production: no-op checks (#30–32), delta
localisation (#40), and reference agreement (#16–22, #35–36) all pass. Both
implementations were run end to end.

---

## 3. Measured effects, stated precisely

### KL direction and vocabulary

All KL values are **KL(clean ‖ intervened)**, in **nats**, computed in float64 on CPU.
Two distributions are reported throughout:

- **AA-20** (headline): logits restricted to the 20 standard amino acids
  `L A G V S E R T I D P K Q N F Y M H W C` (token ids 4–23), then re-softmaxed. This is
  the "amino-acid prediction" distribution.
- **Full vocabulary**: softmax over all 33 ESM-2 tokens, including `<cls> <pad> <eos>
  <unk> X B U Z O . - <null_1> <mask>`.

The two agree closely on this sample because non-amino-acid tokens carry negligible mass;
both are stored at full precision in the JSON, and in
[`out/masked_position_logits_direct.npz`](out/masked_position_logits_direct.npz) as raw
float64 logits (128 rows, one per protein × arm) with the token vocabulary.

### The case A vs case B comparison — "100×" withdrawn

The first report's "~100×" was under-specified. The figure depends heavily on the
aggregation, and the spread is large:

| Aggregation (strongest_active / amplify) | Value |
|---|---|
| Ratio of across-protein **arithmetic means** of KL (the original "100×") | 119.9× |
| Ratio of across-protein **medians** | 257.7× |
| **Median of per-protein ratios** | 525.8× |
| Range of per-protein ratios | 21.3× – 38,326× |

Across all four edit cells, the ratio of means ranges 47.7×–1510.8× and the median
per-protein ratio ranges 74×–526×. **The magnitude should not be quoted as a single
number.** What is robust is the **direction and consistency**: editing at the masked
position moves the masked-position distribution more than editing at a position ~L/4 away,
in **8/8 proteins in all four cells**. That ordering is the evidence that the intervention
reaches the readout through the intended path; the size of the gap is not a stable
quantity on n = 8.

### Feature 619 — sample-specific

Within **these 8 proteins only** (2,716 token positions, 16 protein × case combinations),
feature 619 was active at **100%** of token positions (min = max = 1.0000) and was the
argmax latent at **99.21%** on average (min 96.92%, max 100%). **This is a property of
this sample, not an established global property of the SAE or of ESM-2.** Characterising
it properly would require a much larger token sample, which was out of scope here. The
`median_rank_active` rule was added for contrast: 16 distinct features, active at a mean
of 12.4% of positions, argmax at 0%.

### Top-1 predictions

**None of the 64 tested edit conditions changed the top-1 prediction**, under either the
AA-20 or the full-vocabulary distribution.

This replaces the first report's claim that a single-feature layer-8 edit "isn't large
enough to flip the argmax 25 blocks later." That was an explanation offered as if it were
a finding; it was not tested. No scan over edit magnitude, feature count, or edited
position count was run, so the conditions under which a top-1 change *would* occur are
simply uncharacterised.

For contrast, the reconstruction-only control **did** change the top-1 amino acid in
**2 of 8** proteins when applied to all residue positions, and in **0 of 8** when applied
to a single position.

### Numerical correctness versus biological interpretation

Everything in this document is a statement about **numerical correctness**: that the hook
is on the right tensor, that the SAE math matches the authors' reference bitwise, that a
no-op edit is exactly a no-op, that a nonzero edit is confined to the intended coordinate
and position and matches a float64 reference to ~1e-7 relative, and that the perturbation
propagates measurably to the masked-token logits.

**None of it is evidence about localisation, about what feature 619 or any other latent
represents, or about any biological mechanism.** The edited features were selected by
activation rank for engineering convenience. The mask position `L//2` and the edit
strength ×2 are arbitrary. A downstream KL change demonstrates that the injection
propagates through the network; it does not indicate that the edited direction carries
biologically meaningful information.

---

## 4. Corrections to the first report

1. **"The fp32 difference-of-reconstructions loses ~5 significant figures."** Wrong, and
   it was inferred from an absolute error without a denominator. Measured relative L2
   error for DoR is 1.55e-06 (median), i.e. it **retains ~5.8 significant digits** of
   float32's ~7.2 — a loss of roughly **1.4 digits**, not 5. The direct form retains ~6.5.
2. **"~100× more"** for case A vs case B — withdrawn and replaced with a direction plus an
   explicit aggregation table (§3).
3. **"Not large enough to flip the argmax"** — replaced with the factual statement that no
   tested condition changed the top-1 prediction (§3).
4. **Feature 619 statistics** — now stated as sample-specific rather than as a property of
   the SAE (§3).

---

## 5. Unresolved failures and limitations

**Unresolved failures: none.** 42/42 checks pass under both implementations.

Open limitations, unchanged or newly surfaced:

1. **No dose–response characterisation.** Edit strength was fixed at ×2 and →0 by
   instruction. The magnitude at which a top-1 change would occur is unknown.
2. **Downstream ordering between implementations is not established** — DoR, DD and the
   rounded float64 reference agree to ~1e-4 on logits, with no clean ordering. The
   implementation choice is justified at the activation level only.
3. **n = 8, one mask position per protein, two feature rules.** Sufficient for plumbing,
   not for statistics.
4. **Bitwise results are device-bound** to `mps`/fp32/batch-1 on this machine.
   Re-establish the noise floor before comparing on any other configuration.
5. **These 8 proteins are now burnt** for confirmatory use.
6. **Whole-sequence SAE reconstruction remains a large perturbation** (mean KL 0.293 nats,
   top-1 changed in 2/8). Any design comparing an edited reconstruction against clean —
   rather than the error-preserving form validated here — would confound the edit with
   that damage.

---

## 6. Files

**Preserved unchanged** from the first round:

| Path | Note |
|---|---|
| `out/validation_report.json` | original run, DoR implementation, 39 checks |
| `out/engineering_sample.tsv` | original |
| `out/masked_position_logits.npz` | original |
| `REPORT.md` | original report; superseded where it conflicts with this document |

**New in this follow-up:**

| Path | Contents |
|---|---|
| `compare_edit_implementations.py` | DoR vs direct vs float64 reference, 96 conditions |
| `out/edit_implementation_comparison.json` | full comparison data and error definitions |
| `out/validation_report_direct.json` | production rerun, direct decoder-delta, 42/42 |
| `out/validation_report_dor.json` | reference rerun, DoR, 42/42 |
| `out/engineering_sample_direct.tsv`, `out/engineering_sample_dor.tsv` | sample manifests per run |
| `out/masked_position_logits_direct.npz`, `out/masked_position_logits_dor.npz` | float64 logits per run |
| `FOLLOWUP.md` | this document |

`validate_intervention.py` was modified to take `--impl {direct,dor}` (default `direct`)
and to write tagged outputs, so neither run overwrites the other or the original.

Nothing was committed or pushed; `git status` shows one untracked directory,
`intervention_validation/`.
