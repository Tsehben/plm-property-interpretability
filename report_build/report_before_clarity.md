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

To choose a layer, we used probes that also received amino-acid composition. We asked how often
each one ranked an extracellular protein above a cytoplasmic protein on a separate tuning set.
Layer 26 scored highest, with several layers close behind. On the evaluation set, a probe using
layer 26 alone ranked the proteins better than a baseline using amino-acid frequencies and length.

**We also examined individual patterns.** We started at layer 8 because earlier probes could
already distinguish the labels there and we had validated a sparse autoencoder (SAE) for that layer. It breaks a layer's output into individual patterns whose activity we can measure. One of these, feature 1446, was much more active in extracellular than in cytoplasmic proteins (AUROC 0.974). Among proteins with a
recorded signal peptide, it was about equally active on average in both groups.

**Finally, we connected the two measurements.** A feature that tracks the label shows an
association, not that the prediction relies on it. To start testing that, we removed feature
1446's contribution at one position in layer 8, in a pilot of 24 tuning-set proteins, and read
the probe that the sweep selected at layer 26. The probe's score barely changed and no protein
changed its predicted label (section 5).

## How this builds on earlier work

InterPLM showed that SAE features can identify a broader set of biological concepts than
individual model neurons. InterProt used linear probes on ESM-2 SAE features to recover known
sequence patterns associated with localisation, including signals near the start of proteins.
Our experiments ask how well prediction survives when obvious clues are similar, and whether a
later prediction changes when one such feature is removed. \[1, 2\]

# 2. Which proteins did we compare?

DeepLoc 2.0 provides Swiss-Prot records with localisation labels and five predefined groups of
proteins, P0 through P4. The groups were constructed to limit sequence similarity across them. We
kept proteins carrying exactly one of the two labels in our comparison. A protein could also have
a third localisation label. We used the dataset's recorded labels as targets. \[3\]

\FloatBarrier

\needspace{9\baselineskip}

| Role | Partitions | Proteins | Extracellular | Cytoplasmic |
|:--------------------------|:----------|-----:|-----------:|-----------:|
| Train probes, find features | P1 + P2 | 5,190 | 1,336 | 3,854 |
| Compare layers | P3 | 2,548 | 584 | 1,964 |
| Evaluate selected probes | P4 | 2,366 | 571 | 1,795 |
| **Total analysed** | **P1--P4** | **10,104** | **2,491** | **7,613** |

We reserved P0 from fitting, selection and scoring. Of the 10,139 proteins in P1--P4, eight had
been used for technical validation and 27 shared a base accession with P0 records. This left the
10,104 proteins above. P4 was held out of this sweep's training and selection, although it had
been used in earlier project analyses.

## What each source contributed

\FloatBarrier

| Source | How we used it |
|:-------------|:-----------------------------------------------|
| Frozen ESM-2 \[4\] | For each protein, we averaged the layer outputs over all residues seen by the model (up to 1,022) and, separately, over the first 30. ESM-only probes used these two summaries. |
| Composition baseline | For the same two windows, we calculated the frequencies of the 20 standard amino acids and added one length feature. The combined probe used these 41 inputs together with the ESM-2 summaries. |
| Pfam and MMseqs2 \[5\] | Pfam supplied recorded family identifiers. MMseqs2 grouped proteins by sequence similarity so related proteins could be resampled together. |
| N-terminal chemistry | "N-terminal" refers to the start of the protein. We grouped proteins by charge balance and, in a separate analysis, by hydrophobicity (how strongly residues tend to avoid water) in the first 30 residues. |
| UniProt and SignalP 6 | UniProt annotations indicated where a signal peptide had been recorded. SignalP 6 is a separate predictor. It still requires a licensed local installation. |
| WoLF PSORT \[6\] | This external predictor scored generated sequences in a separate pilot. DeepLoc supplied the original labels and WoLF PSORT supplied the pilot scores. |

## How much sequence similarity crossed the groups?

We audited the original cohort before exclusions. At 30% sequence identity and 80% coverage in both
directions, 22 of 2,382 P4 proteins and 5 of 2,556 P3 proteins shared a cluster with P1+P2. At
50% identity the audit found no such overlap. Relatedness within P4 was more substantial: one
cluster contained 16.2% of its extracellular proteins. We therefore resampled whole clusters in
Figures 2 and 5. This accounts for the measured groups while leaving possible relationships
between clusters unresolved.

# 3. Which layers carry the label?

AUROC is the fraction of extracellular--cytoplasmic pairs in which a probe scored the
extracellular protein higher; 0.5 is chance and 1.0 is perfect ranking. Each of the 33 probes was
fitted on P1+P2 and compared on P3. The selection rule used the combined probe, which received the
41 composition inputs together with the ESM-2 summaries. Its highest P3 AUROC selected **layer
26** (0.988768), ahead of layer 8 (0.983810). **The selection was written to disk before any of
this sweep's P4 predictions were computed.** Seventeen of the 33 layers lay within 0.005 AUROC of
layer 26, so it sits on a broad plateau rather than a sharp peak.

![Probe AUROC at each of ESM-2's 33 layers on the P3 tuning set (584 extracellular, 1,964
cytoplasmic). Probes were fitted on P1+P2. The orange combined-probe curve was used to select
layer 26. The dashed line is the layer-independent composition baseline (0.958956). No intervals
were computed for this curve.](figs/f1_layer_curve.png)

On the P4 evaluation set, the layer-26 ESM-only probe reached an AUROC of 0.982902 against 0.969604
for the composition baseline: a difference of +0.0133 with a 95% interval of [+0.0047, +0.0237].
Comparing layer 26 with layer 8 directly, using the combined probe, gave a difference of +0.020102
with an interval of [-0.003503, +0.063002]. That interval includes zero, so the two layers are
unresolved; it is not evidence that they are equivalent. All intervals in this section are
pointwise, with no correction for multiplicity.

![Layer-26 ESM-only AUROC minus composition-baseline AUROC in four sets of proteins or protein
pairs. Intervals (95%) come from resampling whole sequence clusters, with the same draw for both
models. The orange interval includes zero. The rows share AUROC units but use different proteins
or pairs, and charge and hydrophobicity were narrowed in separate
comparisons.](figs/f2_p4_comparisons.png)

\FloatBarrier

\needspace{9\baselineskip}

| Row | What it answers | $\Delta$ AUROC | 95% interval |
|:----------------|:--------------------------|-----------:|:----------------|
| All P4 proteins | Does ESM beat composition overall? | **+0.0133** | [+0.0047, +0.0237] |
| Unseen recorded family | Does it survive when families are new? | **+0.0115** | [-0.0054, +0.0372] |
| Same charge band | Does it survive fixing N-terminal charge? | **+0.0177** | [+0.0081, +0.0293] |
| Same hydrophobicity band | Does it survive fixing N-terminal hydrophobicity? | **+0.0664** | [+0.0391, +0.0995] |

**Defining "unseen".** A P4 protein is in the unseen-family row when it carries at least one
applicable Pfam signature whose own entry type is *family*, and none of its recorded family-type
Pfam IDs appears in P1, P2 or P3. That is 443 proteins (183 extracellular, 260 cytoplasmic). On
these proteins the absolute AUROCs were 0.9865 for the ESM-only probe and 0.9750 for the
composition baseline; their difference of +0.0115 has the interval above, which includes zero.
Family status was unknown for 1,633 of the 2,366 P4 proteins, mostly proteins with Pfam domain
signatures but no family-type signature, so the unseen-family row describes only the proteins for
which novelty could be determined.

**The two chemistry rows answer different questions.** Charge and hydrophobicity were examined in
separate comparisons, and each compares a different set of protein pairs. The difference in size
between those rows cannot be attributed to either property.

# 4. Which features distinguish the labels?

The layer sweep tested each layer's whole output. To look inside one layer, we used the InterProt
SAE we had already validated for layer 8. Each feature has an identifying number and an activation
that measures how strongly its learned pattern appears. For the 5,190 P1+P2 discovery proteins, we
averaged each feature's activation over the first 30 residues.

We ranked 2,867 eligible features by how far their AUROC lay from 0.5 in either direction. With
extracellular as the positive class, feature 3900's low AUROC means that greater activity points
toward the cytoplasmic label. The separate P1 and P2 rankings shared 19 of their top 20 features.

\FloatBarrier

\needspace{9\baselineskip}

| Rank | SAE feature | Higher mean activity | Discovery AUROC |
|---:|---:|:--------------|-----------:|
| 1 | 3900 | Cytoplasmic | 0.023011 |
| 2 | 642 | Extracellular | 0.974461 |
| 3 | 1446 | Extracellular | 0.974143 |
| 4 | 243 | Extracellular | 0.971749 |
| 5 | 963 | Extracellular | 0.971066 |

Why follow feature 1446? It almost tied 642 on class separation and showed a sharper contrast at
the annotated signal peptide. For 1446, mean activity was 5.587 inside the recorded interval and
0.123 outside it, with higher activity inside for 1,305 of 1,306 proteins. For 642 the means were
2.145 and 0.292. We therefore compared 1446 across groups that share an annotation but differ in
their localisation label.

![Of the 5,009 discovery proteins with verifiable annotation status, the two annotated groups had
similar mean activity (left) and similar percentages with any activity in the first 30 residues
(right). Proteins without an annotation may still have a signal
peptide.](figs/f3_feature1446_groups.png)

Across all discovery proteins, mean 1446 activity was 4.357 for the extracellular group and 0.158
for the cytoplasmic group. Among annotated proteins the means were 4.48 (1,252 extracellular) and
4.45 (59 cytoplasmic). Of those 59 annotations, 53 were predictor-derived, four had experimental
evidence and two were assigned by similarity. This pattern links 1446 to a sequence clue. Whether
the localisation probe uses that clue is a separate question.

# 5. Does removing 1446 change the probe's prediction?

We tested 24 P3 proteins, 12 from each class. With each sequence held fixed, we set feature 1446's
activity to zero in layer 8 at one position at a time (10, 20 or 30). This is the ablation. After
the intervening 18 layers processed the change, we read the localisation probe at layer 26. A
higher score, measured in log-odds, means the probe favours the extracellular label more strongly.

\FloatBarrier

\needspace{10\baselineskip}

| Condition | What the comparison tells us |
|:------------|:-------------------------------------------|
| Clean | The starting score without any intervention. Clean probabilities reproduced the saved predictions to within 0.000051. |
| No-op | We ran the editing procedure while retaining the original activation. The score changed by exactly zero. |
| Ablation | Subtracting the clean score from the score after ablation gives the change associated with removing 1446. |
| Random control | We also applied a perturbation at the same position and of the same size in a random direction. We compared its score with the ablation score. |

![Each bar compares the ablation with the random control for one protein, averaged over the tested
positions. Diamonds compare ablation with clean. Of the 24 differences, nine were positive, four
negative and eleven zero. In those 11 cases, feature 1446 was inactive at all three
positions.](figs/f4_ablation_24.png)

The feature was active in 13 proteins (12 extracellular, one cytoplasmic). Across all 24, removing
it raised the mean probe score by +0.004534 log-odds versus clean. The mean ablation-minus-random
difference was +0.005968, with a median of zero. AUROC was 1.0 before and after ablation. At a
predicted-probability threshold of 0.5, 0 of 24 proteins changed class.

If the feature always pushed the probe toward the extracellular label, removing it would lower the
score. The pilot average moved slightly upward. We found no control feature with comparable
activation frequency and strength, so the relative effect of removing a similarly active feature
remains open.

# 6. What changed beyond the probe score?

## Did ESM-2 change its amino-acid predictions?

In a separate experiment on 400 P4 proteins, we hid one residue and asked ESM-2 which amino acid
belonged there. Two other layer-8 features, 2564 and 3900, had usable control features. We removed
or doubled each candidate and compared the change in the model's predicted probabilities with the
change after an equally large control edit. We measured the shift from the unedited distribution
using Kullback--Leibler (KL) divergence, in nats.

![Positive values mean the model's probabilities shifted more for the candidate than for its
matched controls. All four 98.75% intervals include zero after resampling sequence clusters and
accounting for four comparisons. Units are $10^{-5}$ nats.](figs/f5_masked_token.png)

For feature 2564 ablation, the adjusted interval originally excluded zero when proteins were
resampled individually. Only 61 proteins had a feasible match, and 22 belonged to one cluster. The
paired mean was +2.10 $\times$ 10$^{-5}$ nats. Resampling clusters gave a 98.75% interval of
[-0.19, +3.5] $\times$ 10$^{-5}$. All four contrasts remain inconclusive under this retrospective
analysis.

## Did an external predictor change its localisation score?

For 40 P3 animal proteins, we sampled replacement residues from ESM-2 and scored the resulting
sequences with WoLF PSORT \[6\]. We compared sampling with feature 1446 increased or removed
against unmodified sampling and equal-size random-direction interventions. We used four paired
replicates per protein (160 pairs per comparison).

The WoLF PSORT margin is its extracellular score minus its cytoplasmic score. We measured the
change from each original protein to each generated variant, then compared those changes across
paired conditions. Ordinary resampling with the unmodified model raised the margin by +0.93 on
average. The additional effect of a feature edit is measured by the paired comparison.

All four median paired differences across proteins were zero. Depending on the comparison, 132--139
of 160 pairs produced identical sequences in the two conditions. Five original proteins exactly
matched sequences in the predictor's reference set. The pilot measured how this predictor
responded to the sampled sequences. Where the proteins actually went was never measured.

# 7. Findings to date and next steps

**What is complete.** The 33-layer probe sweep is finished. Layer 26 was selected by a prespecified
rule from a broad plateau, and its ESM-only probe ranked P4 proteins better than the composition
baseline (+0.0133 AUROC). That advantage persisted within charge bands and, in a separate
comparison, within hydrophobicity bands. It was not resolved for the 443 proteins with unseen
recorded families, and layer 26 against layer 8 remains unresolved. Individual SAE features have
been screened at layer 8: feature 3900 pointed most strongly toward the cytoplasmic label and 642
toward the extracellular label, with 1446 close behind.

**What the feature and pilot results show.** Feature 1446 was about equally active, on average, in
annotated extracellular and annotated cytoplasmic proteins (4.48 and 4.45), so it is linked to a
sequence clue as well as to the label. Removing it at one position moved the probe's mean score
slightly upward, by a small and mixed-sign amount, with no change to any predicted class.

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

**Still outstanding.** SignalP 6 predictions, which require a licensed local installation, and
separate layer-wise concept probes. The within-property comparisons above are complete, but they
re-score an existing localisation probe within subsets; fitting separate probes for annotated
concepts at each layer has not been done.

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
