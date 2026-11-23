---
title: "Localisation in ESM-2"
subtitle: "Project check-in"
date: "23 November 2026"
geometry: "margin=2.4cm"
fontsize: 10pt
linkcolor: "black"
urlcolor: "black"
colorlinks: true
header-includes: |
  \usepackage{float}
  \makeatletter\renewcommand*{\fps@figure}{tbp}\makeatother
  \usepackage[section]{placeins}
  \usepackage{needspace}
  \usepackage{flafter}
---

# 1. Does the model recognise where a protein goes, or clues associated with that destination?

We compared proteins labelled extracellular (outside the cell) with proteins labelled cytoplasmic
(in the cytoplasm). One clue that often distinguishes them is a signal peptide, a short segment
near the start of a protein that can direct it into the secretory pathway. A model that
recognises that clue could predict many labels correctly. What happens when the clue and the
label point in different directions?

## How we approached the question

**First, we compared the layers.** ESM-2 is a protein language model that processes each
amino-acid sequence through 33 successive layers. Each layer produces a numerical description of
the sequence. We kept ESM-2 fixed and trained a separate classifier, called a linear probe, on the
output of each layer. Each probe gave every protein a single score: higher scores favoured the
extracellular label and lower scores the cytoplasmic label.

To choose a layer, we used probes that also received amino-acid composition. On a separate tuning
set, we measured how often each probe ranked an extracellular protein above a cytoplasmic one (a
measure called AUROC). Layer 26 scored highest, with several layers close behind. On the
evaluation set, a probe using only layer 26's output, without composition inputs, ranked the
proteins better than a baseline using only amino-acid frequencies and length.

**We also examined individual patterns.** We started at layer 8 because earlier probes could
already distinguish the labels there and we had validated a sparse autoencoder (SAE) for that
layer. It breaks a layer's output into individual patterns whose activity we can measure. One of
these, feature 1446, was much more active in extracellular than in cytoplasmic proteins (AUROC
0.974). Among proteins with a recorded signal peptide, it was about equally active on average in
both groups.

**Finally, we connected the two measurements.** A feature that tracks the label shows an
association, not that the prediction relies on it. To start testing that, we removed feature
1446's contribution at one position in layer 8, in a pilot of 24 tuning-set proteins, and read
the probe that the sweep selected at layer 26. The probe's score barely changed and no protein
changed its predicted label (section 5).

## How this builds on earlier work

InterPLM showed that SAE features can identify a broader set of biological concepts than
individual model neurons \[1\]. InterProt used linear probes on ESM-2 SAE features to recover known
sequence patterns associated with localisation, including signals near the start of proteins
\[2\]. Both show that such patterns can be found in the model. We ask two further questions.
Does the probe still beat amino-acid composition when the proteins being compared share the same
obvious clue, or come from families it never saw? And does the prediction change when one feature
linked to such a clue is removed?

# 2. Which proteins did we compare?

DeepLoc 2.0 provides Swiss-Prot records with localisation labels and five predefined partitions
(groups of proteins), P0 through P4 \[3\]. The partitions were constructed to limit sequence
similarity between them. We kept proteins carrying exactly one of our two labels, extracellular or
cytoplasmic. A kept protein could still carry other localisation labels, such as nucleus. We used
the dataset's recorded labels as targets.

\FloatBarrier

\needspace{9\baselineskip}

| Role | Partitions | Proteins | Extracellular | Cytoplasmic |
|:--------------------------|:----------|-----:|-----------:|-----------:|
| Train probes, find features | P1 + P2 | 5,190 | 1,336 | 3,854 |
| Compare layers | P3 | 2,548 | 584 | 1,964 |
| Evaluate selected probes | P4 | 2,366 | 571 | 1,795 |
| **Total analysed** | **P1--P4** | **10,104** | **2,491** | **7,613** |

We reserved P0 from fitting, selection and scoring. Of the 10,139 proteins in P1--P4, we excluded
35: eight had been used for technical validation, and 27 had another version of the same protein
(a different isoform or the canonical form) in P0. This left the 10,104 proteins above. P4 was held
out of this sweep's training and selection, although it had been used in earlier project analyses.

## What each source contributed

\FloatBarrier

| Source | How we used it |
|:-------------|:-----------------------------------------------|
| Frozen ESM-2 \[4\] | For each protein, we averaged the layer outputs over all residues seen by the model (up to 1,022) and, separately, over the first 30. ESM-only probes used these two summaries. |
| Composition baseline | For the same two windows, we calculated the frequencies of the 20 standard amino acids, and added one length feature. The composition baseline used only these 41 inputs. The combined probe used them together with the ESM-2 summaries. |
| Pfam and MMseqs2 \[5\] | Pfam supplied recorded family identifiers, which we used to find P4 proteins from families the probe never saw. MMseqs2 grouped proteins by sequence similarity so related proteins could be resampled together. |
| N-terminal chemistry | "N-terminal" refers to the start of the protein. Over the first 30 residues we computed charge balance (lysines plus arginines minus aspartates and glutamates, divided by the window length) and mean hydrophobicity (how strongly residues tend to avoid water). For each property, we split proteins into five bands, with boundaries set on P1+P2 only. |
| UniProt and SignalP 6 | UniProt annotations indicated where a signal peptide had been recorded. SignalP 6 is a separate signal-peptide predictor. We have not run it, because it requires a licensed local installation. |
| WoLF PSORT \[6\] | This external localisation predictor scored generated sequences in a separate pilot (section 6). DeepLoc supplied the original labels and WoLF PSORT supplied the pilot scores. |

## How much sequence similarity crossed the groups?

We audited all 10,139 P1--P4 proteins, before the 35 exclusions. At 30% sequence identity and 80%
coverage in both directions, 22 of 2,382 P4 proteins and 5 of 2,556 P3 proteins shared a cluster
of similar sequences with P1+P2. At 50% identity the audit found no such overlap. Relatedness within
P4 was more substantial: one cluster contained 16.2% of P4's extracellular proteins. We therefore
resampled whole clusters, rather than single proteins, for the intervals in Figures 2 and 4. This
stops one large family from counting as many independent proteins, but clusters that are related
to each other are still treated as independent.

# 3. Which layers predict the label best?

AUROC is the fraction of extracellular--cytoplasmic pairs in which a probe scored the
extracellular protein higher; 0.5 is chance and 1.0 is perfect ranking. At each of the 33 layers,
we fitted probes on P1+P2 and compared them on P3. The selection rule used the combined probe,
which received the 41 composition inputs together with the ESM-2 summaries. Its highest P3 AUROC
selected **layer 26** (0.988768), ahead of layer 8 (0.983810). **The choice was fixed before any of
this sweep's P4 results were computed.** Seventeen of the 33 layers lay within 0.005 AUROC of
layer 26, so it sits on a broad plateau rather than a sharp peak.

![Probe AUROC at each of ESM-2's 33 layers on the P3 tuning set (584 extracellular, 1,964
cytoplasmic). Probes were fitted on P1+P2. The orange combined-probe curve was used to select
layer 26; the blue curve is the ESM-only probe. The dashed line is the composition baseline
(0.958956), which does not use ESM-2 and so is the same at every layer. No intervals were computed
for this curve.](figs/f1_layer_curve.png)

On the P4 evaluation set, the layer-26 ESM-only probe reached an AUROC of 0.982902 against 0.969604
for the composition baseline: a difference of +0.0133, with a 95% interval of [+0.0047, +0.0237]
from resampling whole sequence clusters. We compare the ESM-only probe here, rather than the
combined probe used for selection, so that ESM-2's information is set directly against
composition. Adding the composition inputs to layer 26 barely changed the result: the combined
probe's difference was +0.013308 against +0.013298 for ESM-only.

Comparing layer 26 with layer 8 directly, using the combined probe, gave a difference of +0.020102
with an interval of [-0.003503, +0.063002]. That interval includes zero, so the two layers are
unresolved; it is not evidence that they are equivalent. Each interval in this section was
computed on its own, with no adjustment for making several comparisons.

![Layer-26 ESM-only AUROC minus composition-baseline AUROC in four sets of proteins or protein
pairs. Intervals (95%) come from resampling whole sequence clusters, with the same draw for both
models. The orange interval includes zero. The rows share AUROC units but use different proteins
or pairs, and charge and hydrophobicity were narrowed in separate
comparisons.](figs/f2_p4_comparisons.png)

\FloatBarrier

\needspace{9\baselineskip}

| Row | What it answers | $\Delta$ AUROC | 95% interval |
|:----------------|:--------------------------|-----------:|:----------------|
| All P4 proteins | Does ESM-2 beat composition overall? | **+0.0133** | [+0.0047, +0.0237] |
| Unseen recorded family | Does the advantage hold for families the probe never saw? | **+0.0115** | [-0.0054, +0.0372] |
| Same charge band | Does it hold among proteins with similar N-terminal charge? | **+0.0177** | [+0.0081, +0.0293] |
| Same hydrophobicity band | Does it hold among proteins with similar N-terminal hydrophobicity? | **+0.0664** | [+0.0391, +0.0995] |

**Defining "unseen".** A P4 protein counts as unseen-family when it has at least one Pfam family
assignment that applies to its sequence, and none of its Pfam family assignments appears in P1, P2
or P3. We used only Pfam entries classified as families, not those classified as domains. That
gives 443 proteins (183 extracellular, 260 cytoplasmic). The row is uneven: 104 of its 183 extracellular
proteins (57%) belong to one family, cobra-type snake toxins (Pfam PF21947), so it largely
tests that family against a mix of cytoplasmic families. On these proteins the AUROCs were 0.9865
for the ESM-only probe and 0.9750 for the composition baseline. Their difference of +0.0115 has an
interval that includes zero, so the advantage is not resolved for this group. Family status was
unknown for 1,633 of the 2,366 P4 proteins, mostly proteins with Pfam domain assignments but no
family assignment, so this row covers only proteins whose novelty could be determined.

**The chemistry rows.** In these rows, each extracellular protein is compared only with
cytoplasmic proteins from the same band, so a probe cannot win simply by reading charge or
hydrophobicity. The advantage held in both. Charge and hydrophobicity were examined in separate
comparisons, and each uses a different set of protein pairs. So the larger value in the
hydrophobicity row does not mean the probe relies more on hydrophobicity than on charge.

# 4. Which features distinguish the labels?

The layer sweep tested each layer's whole output. To look inside one layer, we used the InterProt
SAE we had already validated for layer 8. Each feature has an identifying number and an activation
that measures how strongly its learned pattern appears. For the 5,190 P1+P2 discovery proteins, we
averaged each feature's activation over the first 30 residues.

We kept the 2,867 features that were active in at least 5% of discovery proteins. For each, we
computed AUROC using its average activation as the score, and ranked features by how far that AUROC
lay from 0.5 in either direction. With extracellular as the positive class, an AUROC near 1 means
the feature is more active in extracellular proteins, and an AUROC near 0, like feature 3900's,
means it is more active in cytoplasmic proteins. Rankings computed separately on P1 and on P2
shared 19 of their top 20 features.

\FloatBarrier

\needspace{9\baselineskip}

| Rank | SAE feature | More active in | Discovery AUROC |
|---:|---:|:--------------|-----------:|
| 1 | 3900 | Cytoplasmic | 0.023011 |
| 2 | 642 | Extracellular | 0.974461 |
| 3 | 1446 | Extracellular | 0.974143 |
| 4 | 243 | Extracellular | 0.971749 |
| 5 | 963 | Extracellular | 0.971066 |

Why follow feature 1446? Its AUROC almost tied 642's, but its activity was far more concentrated
inside the recorded signal-peptide region. Among the 1,306 discovery proteins whose recorded
signal-peptide position could be mapped onto the sequence, 1446's mean activity was 5.587 inside
that region and 0.123 outside it, and higher inside for 1,305 of the 1,306. For 642 the means were
2.145 and 0.292, higher inside for 1,289. Because 1446 is so tightly tied to a known clue, we know
what we are removing when we remove it. We therefore compared 1446 across groups that share an
annotation but differ in their localisation label.

![Feature 1446 in four groups of discovery proteins, split by localisation label (colour) and by
whether UniProt records a signal peptide ("annotated"). Left: mean activity over the first 30
residues. Right: percentage of proteins with any activity there. The groups cover the 5,009
discovery proteins whose UniProt record matched the sequence we used.](figs/f3_feature1446_groups.png)

**How to read Figure 3.** Across all discovery proteins, mean 1446 activity was 4.357 for the
extracellular group and 0.158 for the cytoplasmic group. The figure asks whether that gap follows
the label or the signal peptide. The first two bars share a recorded signal peptide but differ in
label, and 1446 behaved almost identically in them: mean 4.48 and 4.45, active in 99.8% and 98.3%
of proteins. The first and third bars share the extracellular label but differ in annotation, and
without a recorded signal peptide 1446 fell to a mean of 0.95, active in 50.9%. So 1446 follows
the signal peptide more closely than the label. Among proteins without a recorded signal peptide,
extracellular ones still scored above cytoplasmic ones (0.95 against 0.09; 50.9% against 21.5%),
possibly because some carry a signal peptide that UniProt does not record; we have not checked
this. The second bar rests on a small group whose annotations are mostly predicted: of its 59
signal peptides, 53 came from computational prediction, four from experiments and two from
similarity to other proteins. Whether the localisation probe uses the signal-peptide clue is a
separate question, which section 5 starts to test.

# 5. Does removing 1446 change the probe's prediction?

We reused the 24 P3 proteins from an earlier pilot, 12 from each label, chosen in an arbitrary but
reproducible order from proteins 60--1,022 residues long with only standard amino acids. For each
protein we set feature 1446's activity to zero in layer 8 at one position (10, 20 or 30, near the
start of the protein where signal peptides lie). This is the ablation. We did this three times,
once per position, and averaged the three results. After the remaining 18 layers processed the
change, we read the layer-26 ESM-only probe's score. A higher score, measured in log-odds, means
the probe favours the extracellular label more strongly.

We compared each ablated score with two references: the clean score, with no edit, and a random
control, an edit to the layer-8 output of the same size as the ablation but in a random direction.
Running the editing code without changing anything left every score exactly unchanged, so the code
itself adds nothing.

Feature 1446 was active at the tested positions in only 13 proteins (12 extracellular, one
cytoplasmic); in the other 11 there was nothing to remove, so their scores did not change. In all
12 extracellular proteins where it was active, removing 1446 raised the score slightly, by at most
0.039 log-odds; in the one cytoplasmic protein it lowered the score by 0.002. Across all 24, the mean
rise was 0.0045 log-odds. A random edit of the same size had a similarly tiny effect: the mean
ablation-minus-random difference was 0.0060, with a median of zero. Every protein was already
classified with more than 99% confidence before the edit, so changes this small left the predicted
probabilities essentially unchanged, and no protein changed its predicted label.

**Answer: in this pilot, no.** If the probe relied on 1446 to call proteins extracellular, removing
it would have lowered the score. Instead the score barely moved, slightly upward, and no more than
after a random edit. Three things limit what this shows. The test effectively covers 13 proteins,
almost all extracellular, and all far from the probe's decision boundary. The edit was small, one feature at one position, and other features carry
the same clue (642 and 243 are also concentrated inside the signal-peptide region), so a near-zero
change does not show that the probe ignores that clue. And no other SAE feature matched 1446's
activation frequency and strength, so the random-direction edit was the only control available; we
could not compare removing 1446 with removing a comparably active feature.

Section 5 measured only one thing: the score of the probe we trained. A feature edit could change
the model in ways this probe does not pick up. Section 6 reports two other measurements, from pilots
run before this one: whether the edits changed which amino acids ESM-2 predicts, and whether they
changed the localisation score that an independent tool, WoLF PSORT, gives to sequences ESM-2
generates.

# 6. What changed beyond the probe score?

## Did ESM-2 change its amino-acid predictions?

The first pilot used 400 P4 proteins, 200 from each label. In each protein we hid one residue near
the start and asked ESM-2 which amino acid belonged there, as probabilities over the 20 amino acids.
We then edited a layer-8 feature at that residue, removing or doubling it, and measured how far the
probabilities moved from the unedited ones, using Kullback--Leibler (KL) divergence in nats: 0 means
no change, and larger values mean a bigger shift.

Any edit moves the probabilities somewhat, so each feature was compared with control features: other
layer-8 features with similar activation frequency and strength, edited by the same amount. For each
protein we took the feature's shift minus its controls' shift; a positive difference means the
feature moved ESM-2's predictions more than its controls did. Three features were planned. Two came
from the top-20 list in section 4: 3900, ranked first and more active in cytoplasmic proteins, and
2564, ranked 14th and more active in extracellular proteins. The third, 1446, had no matching control
feature and could not be compared; this is the missing control mentioned in section 5. Each
comparison used only proteins where the feature was active at the hidden residue and a matched
control edit could be built, so the counts in Figure 4 range from 61 to 301 of the 400.

That gave four comparisons (two features, removed or doubled). Because we made four, each uses a
98.75% interval rather than 95%, which keeps the overall chance of a false alarm near 5%. Intervals
resample clusters of related sequences, since related proteins are not independent. All four
intervals include zero, and even the largest average difference was below $10^{-4}$ nats. For
removing 2564, an earlier analysis that resampled proteins individually gave an interval excluding
zero. But only 61 proteins qualified, and 22 came from one cluster of related sequences; resampling
clusters gave [$-0.19$, $+3.5$] $\times$ 10$^{-5}$ nats, which includes zero.

![Each row is one feature, removed (×0) or doubled (×2). The dot is the average, across proteins, of
the feature's KL shift minus its controls' shift; positive means the feature moved ESM-2's
amino-acid predictions more than its controls did. The line is a 98.75% interval. Units are
$10^{-5}$ nats.](figs/f5_masked_token.png)

## Did an independent predictor change its localisation score?

The second pilot used 40 P3 animal proteins, 20 from each label, so that WoLF PSORT \[6\], a tool
that predicts localisation from sequence, could use one organism setting throughout. We hid ten
residues near the start of each protein and let ESM-2 fill them in by sampling from its predicted
probabilities, producing a new sequence. We did this with no edit, with 1446 doubled, with 1446
removed, and with random-direction edits of the same size, four times per protein (160 pairs per
comparison).

WoLF PSORT's margin is its extracellular score minus its cytoplasmic score. Refilling the residues
with no edit already raised the margin by 0.93 on average, so we measured 1446's effect against
that: doubled or removed versus no edit, and versus a random edit.

All four comparisons had a median difference of exactly zero. The main reason is that the edit
rarely changed what ESM-2 generated: in 132--139 of 160 pairs, the two runs produced identical
sequences. Five of the 40 original proteins are also in WoLF PSORT's own reference set, so it may
already have seen them.

**Answer: nothing detectable in either pilot.** Both were small, and the second had little to detect
because the edits rarely changed the generated sequence. Neither measured where any protein actually
goes.

# 7. Findings to date and next steps

**What is complete.** The 33-layer probe sweep is finished. Layer 26 was selected by a prespecified
rule from a broad plateau, and its ESM-only probe ranked P4 proteins better than the composition
baseline (+0.0133 AUROC). That advantage persisted within charge bands and, in a separate
comparison, within hydrophobicity bands. It was not resolved for the 443 proteins with unseen
recorded families, and layer 26 against layer 8 remains unresolved. Individual SAE features have
been screened at layer 8: feature 3900 pointed most strongly toward the cytoplasmic label and 642
toward the extracellular label, with 1446 close behind.

**What the feature and pilot results show.** Feature 1446 was about equally active, on average, in
annotated extracellular and annotated cytoplasmic proteins (4.48 and 4.45), so it is linked to the
signal-peptide annotation as well as to the label. Removing it at one position moved the probe's
mean score slightly upward, with some proteins moving up and some down, and no change to any
predicted class.

**What remains open.** Whether ESM-2 represents localisation directly, or represents clues that
accompany it, is unresolved. Everything here is correlational, or measured on a classifier fitted
within the project, and none of it measures where a protein actually goes.

**Next research question.** We will extend the contrastive feature search to other released SAE
layers. The InterProt release on Hugging Face contains checkpoints for nine distinct layers (4, 8,
12, 16, 20, 24, 28, 32 and 33) across ten files, because layer 24 has two versions. There is no
layer-26 checkpoint in this release, so the layer our sweep selected cannot be studied feature by
feature with a released SAE. We have validated and screened only layer 8 so far. We will ask
whether the features that best separate the labels at other depths show the same link to annotated
signal-peptide regions as at layer 8, or whether different patterns emerge. We will use the same
5,190 discovery proteins, the same first-30-residue summary and the same ranking rule. Each SAE
learns its own features, so feature numbers do not carry over between layers; feature 1446 belongs
to the layer-8 SAE only. We will compare activation patterns and annotation associations instead.

**Still outstanding.** SignalP 6 predictions, which require a licensed local installation. Also
outstanding are concept probes: probes trained at each layer to predict the clues themselves, such
as signal-peptide annotation or charge. The family, charge and hydrophobicity comparisons in
section 3 are complete, but they re-score the existing localisation probe within subsets rather
than training new probes.

# References

\[1\] Simon E, Zou J. InterPLM: discovering interpretable features in protein language models via
sparse autoencoders. *Nature Methods* 22:2107--2117, 2025. doi:10.1038/s41592-025-02836-7.

\[2\] Adams E, Bai L, Lee M, Yu Y, AlQuraishi M. From mechanistic interpretability to mechanistic
biology: training, evaluating, and interpreting sparse autoencoders on protein language models.
*Proceedings of the 42nd International Conference on Machine Learning*, PMLR 267, 2025.
https://proceedings.mlr.press/v267/adams25a.html. SAE weights:
https://huggingface.co/liambai/InterProt-ESM2-SAEs.

\[3\] Thumuluri V, Almagro Armenteros JJ, Johansen AR, Nielsen H, Winther O. DeepLoc 2.0:
multi-label subcellular localization prediction using protein language models. *Nucleic Acids
Research* 50(W1):W228--W234, 2022. doi:10.1093/nar/gkac278.

\[4\] Lin Z, Akin H, Rao R, et al. Evolutionary-scale prediction of atomic-level protein structure
with a language model. *Science* 379:1123--1130, 2023. doi:10.1126/science.ade2574.

\[5\] Mistry J, Chuguransky S, Williams L, et al. Pfam: the protein families database in 2021.
*Nucleic Acids Research* 49(D1):D412--D419, 2021. doi:10.1093/nar/gkaa913. Steinegger M, Söding J.
MMseqs2 enables sensitive protein sequence searching for the analysis of massive data sets.
*Nature Biotechnology* 35:1026--1028, 2017. doi:10.1038/nbt.3988.

\[6\] Horton P, Park K-J, Obayashi T, et al. WoLF PSORT: protein localization predictor. *Nucleic
Acids Research* 35:W585--W587, 2007. doi:10.1093/nar/gkm259.
