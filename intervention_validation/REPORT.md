# Layer-8 SAE intervention validation — engineering check

**Date:** 2026-11-20 · **Status:** 39/39 checks pass · **Runtime:** 32.3 s

**Question.** Can we edit a layer-8 SAE feature inside frozen ESM-2, continue the
remaining transformer computation, and correctly measure masked-amino-acid outputs?

**Answer.** Yes. The machinery is verified. **This establishes plumbing only. It is not
evidence for any localisation mechanism, and the feature edited here was chosen for
engineering convenience, not biological relevance.**

P0 was not read. HPA was not read. No existing file was modified; everything new is under
`intervention_validation/`.

---

## 1. Pass/fail table

| # | Check | Result | Evidence |
|---|---|---|---|
| 1 | SAE checkpoint SHA-256 matches pinned value | **PASS** | `1172503b…e56ebd` |
| 2 | No required encoder/LM-head parameter missing or newly initialised | **PASS** | `missing_keys=[]`, `unexpected_keys=[]`, `mismatched_keys=[]` |
| 3 | `lm_head.decoder.weight` tied to input embedding | **PASS** | tensors identical; `tie_word_embeddings=True` |
| 4 | `lm_head.dense` carries trained weights | **PASS** | mean\|w\|=0.07271, std=0.09215 |
| 5 | `lm_head.bias` non-zero (not fresh init) | **PASS** | max\|b\|=14.59375 |
| 6 | `EsmForMaskedLM` reproduces `EsmModel` `hidden_states[8]` | **PASS** | **bitwise equal**, max diff 0.000e+00 |
| 7 | Exactly one encoder layer output equals `hidden_states[8]` | **PASS** | matching indices `[7]` — unique |
| 8 | `hidden_states[8]` is the output of `encoder.layer[7]` | **PASS** | hook_layer=7 |
| 9 | All remaining transformer blocks run after the hook | **PASS** | `layer[8..32]` = 25 blocks |
| 10 | Hooked layer output structure determined empirically | **PASS** | returns a bare `Tensor`, **not** a tuple |
| 11 | Every eligible protein in development P1–4 (none from P0) | **PASS** | 56 eligible of 64 |
| 12 | Every selected sequence matches recorded SHA-256 | **PASS** | 8/8 |
| 13 | No selected protein truncated at 1022 | **PASS** | 8/8 |
| 14 | Decoder orientation `(d_hidden, d_model)` | **PASS** | `w_dec (4096, 1280)` |
| 15 | Encoder orientation `(d_model, d_hidden)` | **PASS** | `w_enc (1280, 4096)` |
| 16 | Latents match reference `get_acts` **elementwise** | **PASS** | **bitwise equal**, 0.000e+00 |
| 17 | Active-feature support identical | **PASS** | 30,720 active of 983,040 |
| 18 | Exactly k=128 features active per token | **PASS** | min=128, max=128 |
| 19 | Reconstruction matches reference `forward_val` **elementwise** | **PASS** | **bitwise equal**, 0.000e+00 |
| 20 | Pre-activations match reference `encode` elementwise | **PASS** | 0.000e+00 |
| 21 | Per-token µ/σ match reference `encode` | **PASS** | exact |
| 22 | Reference `decode(encode(x))` == `forward_val(x)` | **PASS** | 0.000e+00 |
| 23 | Reconstruction restored to original activation scale | **PASS** | x std 2.8492 → recon std 2.6874 |
| 24 | **NEG. CONTROL** biased σ gives different reconstruction | **PASS** | diff 1.259e-03 ≠ 0 |
| 25 | **NEG. CONTROL** dropping `b_pre` gives different latents | **PASS** | diff 1.304e+00 ≠ 0 |
| 26 | Repeat clean runs reproducible | **PASS** | 3 runs **bitwise identical**; noise floor = 0.000e+00 |
| 27 | Capture-only hook preserves logits bitwise | **PASS** | 8/8 proteins |
| 28 | Identity hook preserves hooked tensor bitwise | **PASS** | 8/8 |
| 29 | Identity hook preserves downstream logits bitwise | **PASS** | 8/8 |
| 30 | Zero-strength edit preserves hooked tensor bitwise | **PASS** | 32/32 arms |
| 31 | Zero-strength edit gives exactly zero activation delta | **PASS** | max = 0.000e+00 |
| 32 | Zero-strength edit preserves downstream logits bitwise | **PASS** | 32/32 arms |
| 33 | Nonzero edits change exactly one latent coordinate | **PASS** | 64/64 arms |
| 34 | No **unintended** top-K eviction or insertion | **PASS** | unintended = 0; intended = 32 (see §5) |
| 35 | Decoded difference == closed form to **double** precision | **PASS** | max rel residual 1.098e-13 |
| 36 | fp32 closed-form residual consistent with rounding | **PASS** | max rel 6.604e-05 (see §5) |
| 37 | Activation delta exactly zero outside edited position | **PASS** | 64/64 arms |
| 38 | Nonzero edits non-trivially perturb the activation | **PASS** | min delta L2 = 1.3574 |
| 39 | Nonzero edits move logits above the noise floor | **PASS** | min 2.447e-03 vs floor 0 |

---

## 2. Environment, model paths and commands

| Field | Value |
|---|---|
| ESM-2 checkpoint | `facebook/esm2_t33_650M_UR50D` |
| Cached revision | `08e4846e537177426273712802403f7ba8261b6c` |
| SAE | `interprot_ref/raw/esm2_plm1280_l8_sae4096.safetensors` |
| SAE SHA-256 | `1172503b90c950651e787a5625d91ff1c702a04ff3d851e8e95b700201e56ebd` (verified at run time) |
| SAE config | `SparseAutoencoder(d_model=1280, d_hidden=4096, k=128)` |
| torch / transformers / numpy | 2.8.0 / 4.57.6 / 2.0.2 |
| Python / platform | 3.9.6 / macOS-26.5.1-arm64 |
| Device / dtype / batch | `mps` / float32 / 1 |
| SAE execution | CPU, float32 (as in the original pilot); float64 twin used only for cross-checks |

```bash
cd "<repo root>"
python3 intervention_validation/validate_intervention.py
```

### Layer indexing, derived from the model rather than assumed

The model has **33 encoder layers** (`encoder.layer[0..32]`) and **34 hidden states**.
HF's `hidden_states[0]` is the embedding output, so `hidden_states[i]` is the output of
`encoder.layer[i-1]`. Sweeping a hook over all 33 layers and comparing each captured
tensor to `hidden_states[8]` found **exactly one match, at `encoder.layer[7]`**.

> **Correction to the earlier audit.** It described the continuation as "blocks 9–32".
> The correct 0-indexed continuation is **`encoder.layer[8..32]`, 25 blocks**, followed by
> `esm.encoder.emb_layer_norm_after` and `lm_head`. The hook site is `layer[7]`,
> not `layer[8]`. Hooking `layer[8]` would have intervened on `hidden_states[9]` —
> a silent off-by-one that no downstream metric would have revealed.

### Loading diagnostics

`EsmForMaskedLM.from_pretrained(..., output_loading_info=True)` returned **empty**
`missing_keys`, `unexpected_keys` and `mismatched_keys`. The LM head is pretrained, not
freshly initialised; `lm_head.decoder.weight` is tied to `esm.embeddings.word_embeddings.weight`
(verified by tensor identity, `tie_word_embeddings=True`), which is the one legitimate
weight-sharing difference. `esm.contact_head` is present and unused by this validation.

### A structural assumption that was wrong

In transformers 4.57.6 `EsmLayer.forward` returns a **bare `Tensor`, not a tuple**. The
hook branches on `isinstance(out, tuple)` and handles both. An implementation that
assumed a tuple and wrote `out[0]` would silently index the batch dimension and corrupt
the tensor shape. The check now records the observed structure rather than asserting a
guess.

---

## 3. Engineering sample

Drawn from `DeepLoc2/sae_pilot/sae_pilot_sample.tsv` (the existing 64-protein pilot:
seed 42, label × length-quartile stratified). **Selection rule:** keep records with
`split_role == development_P1_4` and original `seq_len <= 1022`, sort by accession
ascending, take the first 8. No label, prediction or model output was consulted. Every
sequence was re-verified against its recorded SHA-256; every record was asserted to be
outside partition 0.

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

Mask position = `L//2`, case-B edit position = `L//4`, both as 0-based residue indices;
token index = residue index + 1 (for `<cls>`). `<cls>`, `<eos>` and padding are excluded
from all edits by construction. **Caveat:** the 64-protein source sample is
label-stratified, so these 8 inherit that stratification; only the selection *within* it
is label-free.

> **These 8 accessions are designated ENGINEERING-ONLY and are reserved from future
> confirmatory intervention evaluation.** No new global dataset split was created.

---

## 4. Reference SAE agreement

Compared against the authors' `SparseAutoencoder` on an identical real layer-8 tensor
`(240, 1280)`, elementwise rather than by aggregate error:

| Quantity | Agreement |
|---|---|
| Latent activations vs `get_acts` | **bitwise identical**, max diff 0.000e+00 |
| Pre-activations vs `encode` | **bitwise identical**, 0.000e+00 |
| Per-token µ, σ vs `encode` | **bitwise identical** |
| Reconstruction vs `forward_val` | **bitwise identical**, 0.000e+00 |
| `decode(encode(x))` vs `forward_val(x)` | **bitwise identical** |
| Active support | identical, exactly 128/token |

Two **negative controls** confirm these comparisons can actually fail: substituting a
biased standard deviation changes the reconstruction (1.259e-03), and dropping `b_pre`
changes the latents (1.304e+00). Without these, the agreement above would not be
informative.

Reconstruction error **on this tensor** was 0.106 relative to zero and 0.637 relative to a
mean-vector baseline. Per instruction these were **not** used as a pass criterion and are
**not** compared to the 64-protein pilot aggregates (0.1109 / 0.5178), whose sample,
inputs (unmasked) and aggregation differ.

---

## 5. No-op behaviour and the two corrected criteria

**Noise floor first.** Three identical clean runs were **bitwise identical** — the
measured floor on this device is exactly **0.000e+00**. Every "exact" claim below is
therefore literal equality, not a tolerance, and the floor-crossing check in row 39 is
correspondingly weak (any nonzero change passes it); the magnitudes in §7 carry the real
information.

| No-op arm | Hooked tensor | Downstream logits |
|---|---|---|
| Capture-only hook | — | **bitwise equal** to clean, 8/8 |
| Identity hook | **bitwise preserved**, 8/8 | **bitwise equal**, 8/8 |
| Zero-strength error-preserving edit | **bitwise preserved**, 32/32 | **bitwise equal**, 32/32 |

The zero-strength result is exact by construction: with `edit(a) == a` the two `decode`
calls receive identical inputs through an identical code path, so their difference is
literally `0.0` and `h' = h + 0.0 = h`. This is a weaker test than it appears — it
validates plumbing, not SAE numerics, which is why checks 14–25 exist separately.

### Two first-run failures, both diagnosed rather than tolerated

**(a) "16 top-K evictions."** All 16 were `ablate` arms, exactly one each: the ablated
feature leaving its own top-K support. That *is* ablation. The check was miscounting an
intended effect as an unintended one. It now masks out the edited coordinate and counts
only changes to *other* features — **0 unintended evictions, 0 insertions** across 64
arms, with 32 intended. Amplification caused no support change at all, as expected
(128 non-zeros before and after).

**(b) "Closed-form mismatch, 3.8e-05."** The decoded difference was compared against
`((topK(a_ed) − topK(a)) @ w_dec) * std`, in which `b_pre` and µ cancel. Recomputing the
identical quantities in **float64** gave a max relative residual of **1.098e-13** — the
math is exact; the fp32 gap is catastrophic cancellation from subtracting two large
reconstructions. The criterion is now float64 exactness (≤1e-11), with the fp32 residual
reported as a measured quantity (max relative 6.604e-05) rather than silently tolerated.

**Implication worth carrying forward:** computing the edit as a difference of two
reconstructions loses ~5 significant figures in fp32. The closed form
`delta = ((topK(a_ed) − topK(a)) @ w_dec) * std` is numerically better and should be
preferred when edit magnitudes get small.

---

## 6. Feature selection, and a property worth flagging

The spec asked for one deterministic live feature. The obvious rule — strongest active
latent at the edit position — selected **feature 619 in all 16 protein × case
combinations**. Feature 619 is active at **100% of tokens** and is the argmax at
**99.2%**: a near-always-on, high-magnitude latent. Validating only on it would test one
atypical direction.

A second rule was therefore added: **median-rank active** (rank 64 of 128 by activation).
It selected **16 distinct features**, active at a mean of **12.4%** of tokens and argmax at
**0%** — ordinary features. All checks pass under both rules. Neither feature set carries
any claim of localisation relevance.

---

## 7. Measured effects

### Reconstruction-only controls (no feature edit), n = 8

| Arm | Activation Δ (rel L2) | KL(clean‖arm) AA-20 | Top-1 AA unchanged | max\|Δ logit\| |
|---|---|---|---|---|
| Replace **all residue positions** | mean **0.341** (0.287–0.419) | mean **0.2933** (3.0e-05 – **1.537**) | **6/8** | mean 5.56 (2.34–10.58) |
| Replace **one position** | mean 0.0243 | mean 1.25e-05 | 8/8 | mean 0.266 |

This is the empirical answer to the question left open earlier: **whole-sequence SAE
reconstruction is a large downstream perturbation** — it changes the top-1 amino-acid
prediction at the masked position in 2 of 8 proteins and reaches KL 1.54 nats in the worst
case. It had to be measured; the activation-space figure alone did not establish it.
Single-position replacement is comparatively mild.

### Feature edits (KL measured at the masked position)

| Case | Rule | Arm | Activation Δ L2 | Δ rel to ‖h‖ | KL AA-20 | Top-1 AA unchanged |
|---|---|---|---|---|---|---|
| A (edit at mask) | strongest | amplify ×2 | 33.04 | 0.0234 | **2.15e-03** | 8/8 |
| A | strongest | ablate | 33.04 | 0.0234 | **2.34e-03** | 8/8 |
| A | median-rank | amplify ×2 | 2.01 | 0.00132 | 3.39e-05 | 8/8 |
| A | median-rank | ablate | 2.01 | 0.00132 | 3.24e-05 | 8/8 |
| B (edit elsewhere) | strongest | amplify ×2 | 37.09 | 0.0262 | 1.79e-05 | 8/8 |
| B | strongest | ablate | 37.09 | 0.0262 | 4.91e-05 | 8/8 |
| B | median-rank | amplify ×2 | 2.02 | 0.00145 | 2.24e-08 | 8/8 |
| B | median-rank | ablate | 2.02 | 0.00145 | 2.20e-08 | 8/8 |

Three internal consistency signals, none of which was asserted in advance:

- **Amplify and ablate produce identical perturbation norms** within each cell. Expected:
  amplifying by ×2 adds `+a·w_dec` and ablating adds `−a·w_dec`, equal magnitude,
  opposite sign.
- **Case A moves the masked-position distribution ~100× more than case B** at equal
  perturbation norm (2.15e-03 vs 1.79e-05 for the strongest rule). Editing the masked
  position itself affects its own prediction more than editing a position ~L/4 away —
  the expected ordering, and evidence the intervention is reaching the readout through
  the intended path.
- **Effect scales with perturbation size**: the strongest-feature edit (Δ ≈ 2.3–2.6% of
  ‖h‖) moves the distribution 1–3 orders of magnitude more than the median-rank edit
  (Δ ≈ 0.13%).

**The top-1 amino acid never changed, in any of the 64 edit arms.** Per the brief, this is
reported as a verified injection with a small downstream effect, **not** a plumbing
failure: the injection is confirmed exact against the closed form in float64, confined to
the intended position, confined to one latent coordinate, and it moves the logits by up to
1.33 — far above a noise floor of exactly zero. A single-feature, single-position edit at
layer 8 simply is not large enough to flip the argmax 25 blocks later.

### Metric vocabulary

- **AA-20** (headline): logits restricted to the 20 standard amino acids
  `L A G V S E R T I D P K Q N F Y M H W C` (token ids 4–23), re-softmaxed.
- **Full vocab**: softmax over all 33 ESM-2 tokens, including `<cls> <pad> <eos> <unk>
  X B U Z O . - <null_1> <mask>`.
- KL is **KL(clean ‖ intervened)** in nats, computed in float64 on CPU (MPS has no
  float64). Full-precision values for both vocabularies are in the JSON; the two agree
  closely here because non-amino-acid tokens carry negligible mass.

---

## 8. Unresolved assumptions and limitations

1. **Device-bound exactness.** Bitwise results hold for `mps`/fp32/batch-1 on this
   machine. Changing device, dtype or batch composition will break bitwise equality
   legitimately; re-establish the noise floor before comparing anything there.
2. **n = 8, one mask position per protein.** Sufficient for plumbing, not for any
   statistical claim.
3. **Mask position `L//2` is arbitrary** and not chosen for biological relevance.
4. **Feature 619 is atypical** (always-on). The median-rank rule mitigates but does not
   eliminate the narrowness of the feature sample.
5. **Amplification factor ×2 is arbitrary.** No scan over edit magnitude was run, so the
   dose–response relationship is uncharacterised.
6. **Whole-sequence reconstruction damage is large** (§7). Any future design comparing an
   edited reconstruction against clean — rather than using the error-preserving form
   validated here — would confound the edit with that damage.
7. **The fp32 difference-of-reconstructions loses ~5 significant figures**; prefer the
   closed form for small edits.
8. **No biological claim is supported by anything in this report.**

---

## 9. Files created

All new, none overwritten:

| Path | Contents |
|---|---|
| `intervention_validation/validate_intervention.py` | The validation (single entry point) |
| `intervention_validation/out/validation_report.json` | Full-precision results, all 39 checks, per-protein per-arm records |
| `intervention_validation/out/engineering_sample.tsv` | The 8 engineering-only accessions with partitions and hashes |
| `intervention_validation/out/masked_position_logits.npz` | 128 rows of float64 masked-position logits (one per protein × arm) plus the token vocabulary |
| `intervention_validation/REPORT.md` | This document |

`git status` shows one untracked directory, `intervention_validation/`. Nothing was
committed or pushed.
