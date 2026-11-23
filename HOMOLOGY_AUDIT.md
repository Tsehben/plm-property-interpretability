# Homology Audit — DeepLoc Localisation Development Cohort (P1–P4)

**Status:** complete. Read-only with respect to every earlier result. No P0 access, no label
or partition change, no model fitting, no embedding extraction, no intervention, not
committed, not pushed.

**Purpose.** [`LOCALISATION_EXPERIMENT_SPEC.md`](LOCALISATION_EXPERIMENT_SPEC.md) §2.6 records
that cluster identifiers, pairwise identities and any measurement of residual homology between
partitions were **not available**, and refuses the family-disjointness claim outright. §3.3
names the work that would convert that assumption into a measurement. This audit is that work.
It is the precondition for the §5.6 cluster bootstrap and for any statement about
generalisation in a layer sweep.

---

## 1. Input and tool

### 1.1 Audit population vs. scientific-analysis population

These are deliberately different, and the distinction matters.

| Population | n | What it is | Used here for |
|---|---|---|---|
| **Audit population** | **10,139** | Every P1–P4 development protein | All clustering |
| Scientific-analysis population (`in_core`) | 10,104 | Audit population minus E1/E2 exclusions | Reported separately, as a check |
| Annotated population (`in_annotated`) | 9,745 | Core minus E3 (coordinates not transferable) | Not re-clustered; flag carried |

**All 10,139 are clustered.** The E1/E2/E3 flags are **carried as columns, not applied**. A
protein excluded from scientific analysis is still a sequence that can be homologous to an
included one; dropping it would make cross-partition relatedness look smaller than it is.
Every output table therefore reports both the full 10,139 cohort and the 10,104 core.

Exclusion flags carried (from [`discovery_analysis/v1/inclusion_manifest.tsv`](discovery_analysis/v1/inclusion_manifest.tsv)):
8 × `E1_engineering_protein`, 27 × `E2_p0_base_collision`, 359 × `E3_coordinates_not_transferable`.

### 1.2 Sequences

- **Original full-length DeepLoc sequences**, read from
  [`DeepLoc2/raw/Swissprot_Train_Validation_dataset.csv`](DeepLoc2/raw/Swissprot_Train_Validation_dataset.csv)
  (SHA-256 `29a07b29…8384a2`). **Not** the 1,022-residue ESM-truncated inputs. Lengths run
  40 – 18,141, median 414; 1,111 of the 10,139 exceed 1,022 and would have been truncated for
  ESM-2. Clustering truncated sequences would have fabricated similarity between long proteins
  that share only an N-terminal region.
- **Every sequence verified** byte-for-byte against the frozen `sha256` in the inclusion
  manifest before being written. 10,139/10,139 matched. **0 dropped, 0 substituted.**
- **Non-standard residues:** 22 proteins carry them — U ×59, X ×18, B ×1. Written **verbatim**.
  MMseqs2 maps unrecognised residue codes to `X` at scoring time; we did not pre-edit,
  substitute or remove any residue, and no protein was dropped for carrying one.
- **Isoforms:** 268 accessions carry a `-N` suffix. Each is resolved by **its own accession**,
  exactly as it appears in the DeepLoc CSV. The canonical sequence is **never** substituted for
  an isoform. Isoforms of the same base accession are allowed to cluster together or not on
  their sequences alone, with no special handling.

### 1.3 Tool

No sequence-clustering tool was present. MMseqs2 was installed **locally**; **no sequence left
this machine** — MMseqs2 clustering is entirely local and makes no network call.

```
brew install mmseqs2          # Homebrew bottle, arm64_tahoe
mmseqs version                # 19-ab57c
```

Binary `/opt/homebrew/bin/mmseqs`; machine Darwin arm64, 12 cores. Exact clustering commands are
in [`homology_audit/run_clustering.sh`](homology_audit/run_clustering.sh) and the full console
log in [`homology_audit/v1/clustering_run.log`](homology_audit/v1/clustering_run.log).

---

## 2. Clustering specification — frozen before any cluster output was examined

Written to [`homology_audit/v1/frozen_specification.json`](homology_audit/v1/frozen_specification.json)
(SHA-256 `79b0985c…f208440`) **before** the clustering was run. Nothing below was chosen after
seeing a cluster.

| Item | Setting |
|---|---|
| Algorithm | `mmseqs easy-cluster`, cascaded clustering, `--cluster-mode 0` (Set-Cover, greedy) |
| Reassignment | `--cluster-reassign 1` |
| Sensitivity | `-s 7.5` (sensitive) — 30 % identity is remote; the default `-s 4.0` would under-recall distant pairs |
| E-value | `-e 0.001` |
| Alignment mode | `--alignment-mode 3` |
| **Identity definition** | `--seq-id-mode 0` — identical aligned residue pairs **÷ alignment length** (not ÷ shorter or full sequence) |
| **Coverage requirement** | `-c 0.8 --cov-mode 0` — **bidirectional**: ≥ 80 % of the query **and** ≥ 80 % of the target |
| **Primary grouping** | **30 % identity, ≥ 80 % bidirectional coverage** |
| Secondary | 50 % and 80 % identity, **same** coverage requirement |
| Threads | 12 |

The same settings are used at all three identity thresholds so the three results are comparable.

### 2.1 What these clusters are, and are not

> These are **operational sequence-similarity clusters**. They are **not verified protein
> families** and they are **not proof of evolutionary independence**.

Specific limitations, all of which bear on how §3 should be read:

1. **Bidirectional 80 % coverage demands near-global similarity.** It will **miss** shared
   domains, partial homology, a multi-domain protein that shares one domain with another, and
   homologous regions embedded in proteins of very different length. Two proteins in different
   clusters may still be substantially related.
2. **Set-Cover is greedy and representative-based.** Membership is transitive *through the
   representative*, so a cluster is not guaranteed to be a clique: two members may be less
   similar to each other than each is to the representative.
3. **The clusterings are not strictly hierarchical.** 55 of the id50 clusters and 17 of the id80
   clusters are not nested inside a single id30 cluster. This is expected behaviour for greedy
   cascaded clustering, not an error, but it means "the 30 % clusters" and "the 50 % clusters"
   are two different partitions of the cohort, not a refinement of one another.
4. **Counts are artifacts of this configuration**, not intrinsic properties of the cohort. A
   looser coverage mode or a domain-aware method would give different, larger clusters.
5. **Different clusters ≠ no sequence similarity.** Everything in §3 bounds relatedness *under
   this definition*. It does not establish independence.

---

## 3. Results

### 3.1 Cluster structure, full cohort (10,139)

| | **30 % id** | 50 % id | 80 % id |
|---|---|---|---|
| Clusters | **6,486** | 7,378 | 8,724 |
| Singletons | **4,804** (47.4 % of clusters) | 5,916 | 7,705 |
| Proteins in a cluster of ≥ 2 | **5,335 (52.6 %)** | 4,223 (41.7 %) | 2,434 (24.0 %) |
| Size 2–4 / 5–9 / 10–49 / ≥ 50 | 1,484 / 156 / 41 / 1 | 1,329 / 104 / 28 / 1 | 979 / 36 / 4 / 0 |
| Largest cluster | **93** | 52 | 16 |
| Mean cluster size | 1.563 | 1.374 | 1.162 |

Core (10,104) is nearly identical — 6,464 / 7,354 / 8,694 clusters, same largest cluster,
52.6 % in clusters ≥ 2. **The 35 excluded proteins make no material difference to the homology
picture**, which is itself a useful result: the E1/E2 exclusions are not concentrated in
redundant families.

### 3.2 Cross-role and cross-partition membership — the headline

Roles: **discovery** = P1 + P2 (5,201), **tuning** = P3 (2,556), **evaluation** = P4 (2,382).

| | **30 % id** | 50 % id | 80 % id |
|---|---|---|---|
| Clusters confined to one role | 6,462 | 7,378 | 8,724 |
| Clusters spanning **two** roles | **24** | **0** | **0** |
| Clusters spanning three roles | **0** | 0 | 0 |
| **P3 proteins sharing a cluster with discovery** | **5 / 2,556 = 0.20 %** | **0** | **0** |
| **P4 proteins sharing a cluster with discovery** | **22 / 2,382 = 0.92 %** | **0** | **0** |

Pairwise, at 30 % identity ([`cross_partition_overlap.tsv`](homology_audit/v1/cross_partition_overlap.tsv)):

| Pair | Shared clusters | Proteins involved (A / B) | % of A / % of B |
|---|---|---|---|
| P1–P2 | 7 | 10 / 12 | 0.40 / 0.44 |
| P1–P3 | 2 | 3 / 4 | 0.12 / 0.16 |
| P1–P4 | 5 | 5 / 6 | 0.20 / 0.25 |
| P2–P3 | 1 | 6 / 1 | 0.22 / 0.04 |
| P2–P4 | 9 | 20 / 16 | 0.74 / 0.67 |
| P3–P4 | 7 | 7 / 9 | 0.27 / 0.38 |

At 50 % and 80 % identity, **every pair is zero**.

**Reading.** Under this operational definition, DeepLoc's published homology partitioning holds
up well: cross-partition relatedness at 30 % identity / 80 % bidirectional coverage touches
under 1 % of proteins in any partition, and vanishes entirely at 50 %. §2.6's concern about
*cross-partition* leakage is, at this definition, measured and small.

**This is overlap under this clustering definition. It is not a bound on homology.** The
figures above count proteins that share a 30 %-identity / 80 %-bidirectional-coverage cluster
and nothing else. That criterion cannot see shared-domain or partial homology (§2.1), so
relatedness it does not detect is simply unmeasured here — the numbers do not cap it, and a
domain-aware or lower-coverage method would be expected to report more. DeepLoc's own
partitioning also used a 30 % similarity criterion, so a 30 % clustering partly re-measures
their construction rather than testing it independently. The statement these numbers support
is: *overlap under this definition is small*. They do not support *the partitions are
family-disjoint*, and they do not support any statement about total homology.

### 3.3 The real finding: within-partition redundancy

Cross-partition overlap is small. **Within-partition redundancy is not**, and it is
concentrated asymmetrically by class ([`class_group_counts.tsv`](homology_audit/v1/class_group_counts.tsv), 30 % identity):

| Role | Class | n | Distinct clusters | Singletons | Largest | % in largest | % in clusters ≥ 2 |
|---|---|---|---|---|---|---|---|
| Discovery | Extracellular | 1,339 | 645 | 425 | 45 | 3.4 % | **68.3 %** |
| Discovery | Cytoplasm | 3,862 | 2,688 | 2,077 | 25 | 0.6 % | 46.2 % |
| Tuning | Extracellular | 584 | 336 | 241 | 21 | 3.6 % | 58.7 % |
| Tuning | Cytoplasm | 1,972 | 1,315 | 978 | 27 | 1.4 % | 50.4 % |
| **Evaluation** | **Extracellular** | **575** | **312** | 232 | **93** | **16.2 %** | 59.7 % |
| Evaluation | Cytoplasm | 1,807 | 1,256 | 933 | 12 | 0.7 % | 48.4 % |

Two things stand out.

**(a) Extracellular proteins are consistently more redundant than cytoplasmic ones.** In
discovery, 68.3 % of extracellular proteins share a cluster with another, against 46.2 % of
cytoplasmic. The minority class is also the more correlated class, so a protein-level bootstrap
overstates precision *more* for the class that has fewer proteins to begin with.

**(b) One cluster is 16.2 % of the entire P4 extracellular class.** The largest cluster in the
cohort — representative `Q2VBP4`, 93 proteins — is **entirely extracellular, entirely Metazoa,
and entirely in P4**. Cross-referencing
[`annotation_readiness/out/annotation_primary.tsv`](annotation_readiness/out/annotation_primary.tsv):
84 of the 93 carry Pfam **PF21947** (snake-toxin-like / three-finger toxin), and the organisms
are elapid snakes (*Ophiophagus hannah* ×12, *Micrurus altirostris* ×11, *Bungarus multicinctus*
×10, *Drysdalia coronoides* ×9, *Naja atra* ×7, …). It is a single venom-toxin family that
happens to sit wholly inside the evaluation partition. At 50 % identity it splits into 10
clusters, at 80 % into 54.

### 3.4 Class composition and mixed-class clusters

| | 30 % id | 50 % id | 80 % id |
|---|---|---|---|
| Extracellular-only clusters | 1,244 | 1,508 | 2,005 |
| Cytoplasm-only clusters | 5,199 | 5,849 | 6,708 |
| **Mixed-class clusters** | **43 (0.66 %)** | 21 (0.28 %) | 11 (0.13 %) |
| Proteins in mixed-class clusters | 214 (2.1 %) | 69 (0.7 %) | 30 (0.3 %) |

Mixed-class clusters are rare. Sequence-similarity groups at this definition are
overwhelmingly label-homogeneous — which is exactly why cluster membership is a plausible
confound for any feature that predicts the label, and why cluster-level resampling is the
honest unit.

### 3.5 Completed experiments — cluster membership only

Reported from the frozen samples as they stand. **No outcome was recomputed and no new sample
was drawn.** ([`completed_sample_clusters.tsv`](homology_audit/v1/completed_sample_clusters.tsv))

**400-protein P4 evaluation sample** ([`evaluation_p4/v1/eval_sample.tsv`](evaluation_p4/v1/eval_sample.tsv)), 200 Extracellular / 200 Cytoplasm:

| | 30 % id | 50 % id | 80 % id |
|---|---|---|---|
| Distinct clusters spanned | 322 | 349 | 385 |
| Proteins sharing a cluster with another sample member | 108 | 79 | 30 |
| Largest within-sample cluster | **31** | 11 | 2 |
| Sharing a cluster with **any discovery** protein | **5** | **0** | **0** |
| Singletons within the full cohort | 180 | 228 | 301 |

By class at 30 %: Extracellular 200 proteins → 134 clusters, largest **31**; Cytoplasm 200
proteins → 193 clusters, largest 2. **That 31-member block is part of the `Q2VBP4`
three-finger-toxin family** — 15.5 % of the sample's extracellular half is one venom family.
Stated as what it is: the 200 extracellular proteins fall into **134 clusters at this
definition**, and the largest holds 31 of them. This is a cluster count, **not an effective
sample size** — translating 200 proteins into "134 independent observations" would assume the
134 clusters are independent of each other, which this audit does not establish (§2.1).

Separately, 5 of the 400 share a 30 % cluster with a discovery protein. That is small, but it
is **not zero**, and those five are the proteins for which the evaluation is least independent
of feature selection.

**40-protein sequence-variant pilot** ([`sequence_variant_pilot/v1/pilot_sample.tsv`](sequence_variant_pilot/v1/pilot_sample.tsv)), all P3, 20 Extracellular / 20 Cytoplasm:

| | 30 % id | 50 % id | 80 % id |
|---|---|---|---|
| Distinct clusters spanned | 38 | 39 | 40 |
| Within-sample sharing | 2 pairs (`Q90YK5`, `P61364`) | 1 pair | none |
| Sharing a cluster with discovery | **1** | **0** | **0** |

This sample is close to cluster-independent: 38 clusters for 40 proteins at the primary
threshold, 1 protein linked to discovery.

> **Do not read "different clusters" as "no sequence similarity."** Every figure in §3.5 bounds
> relatedness under 30 %/80 %-bidirectional. Proteins sharing a single domain, or homologous
> over less than 80 % of their length, are counted as unrelated here and may not be.

---

## 4. Recommendation for the layer sweep

Using cluster membership and class counts only. P0 stays reserved throughout; nothing below
touches it.

### 4.1 Keep the existing partitions, and report overlap explicitly — **primary**

Keep discovery = P1+P2, tuning = P3, evaluation = P4 exactly as they are. The measurement now
supports this rather than merely assuming it: **0.92 % of P4 and 0.20 % of P3 share a 30 %
cluster with discovery, and 0 % at 50 %.** That is a reportable number, not an assumption.

What must accompany it:

1. State the overlap figures above in any generalisation claim, with the coverage caveat.
2. **Report the 5 P3 and 22 P4 proteins that do share a discovery cluster**, and give the
   evaluation statistic with and without them. 22 of 2,382 will not move a result; showing that
   it does not is cheap and closes the question.
3. **Use the §5.6 cluster bootstrap, resampling 30 % clusters, as the primary interval.** It is
   now runnable — [`cluster_membership.tsv`](homology_audit/v1/cluster_membership.tsv) supplies
   the resampling unit. The protein-level interval becomes the secondary. §5.6's fallback
   language ("report the protein-level interval with the explicit statement that it assumes an
   independence we have not verified") no longer applies and should not be used.
4. For any P4 evaluation, **report the largest-cluster share of the extracellular class
   (16.2 % cohort-wide, 15.5 % in the completed 400-protein sample)**. A layer sweep that
   reports per-layer extracellular performance on P4 without this is reporting, in part, how
   well each layer handles snake venom toxins.

### 4.2 Supplementary cluster-held-out development evaluation — **secondary, and labelled as such**

**This does not replace the existing splits.** It is an additional, clearly-labelled analysis
run alongside them.

Construct a development-only split that keeps **whole 30 % clusters together**: assign each of
the 6,486 clusters, not each protein, to a fold; stratify assignment on class; keep P0 out.
This targets within-partition family redundancy, which §3.3 shows is the larger effect, rather
than cross-partition overlap, which §3.2 shows is already small under this definition.

> **Cluster-held-out splitting and cluster-balanced weighting are two different fixes, and
> only one of them is this.** Assigning whole clusters to folds prevents a cluster from
> appearing on both sides of the split — it removes *overlap* between fitting and evaluation.
> It does **not** stop a large cluster dominating the metric computed *within* a fold: the
> 93-member toxin family contributes 93 of the protein-weighted observations in whatever fold
> it lands in, exactly as it does now. Removing that domination is a separate, explicit choice
> — weighting each cluster equally rather than each protein — and it answers a different
> question ("how well does this work per family?" rather than "per protein?"). Neither
> substitutes for the other, and any analysis using one should say which.

**It must be labelled as a development-set evaluation**, for a reason that cannot be engineered
away: **every one of these 10,139 proteins has already contributed to earlier analyses** —
discovery feature selection, tuning, the 400-protein P4 experiment, the pilots. Re-partitioning
them does not restore naivety. The result is a *robustness check on redundancy*, not a fresh
held-out estimate, and P0 remains the only genuinely untouched data.

### 4.3 Are there enough clusters to resample? — **yes for both classes, with one caveat**

Distinct 30 % clusters available, by class. **These are cluster counts under the §2 grouping,
not effective sample sizes and not counts of independent observations.**

| Class | Discovery | Tuning | Evaluation | **Cohort-wide** |
|---|---|---|---|---|
| Extracellular | 645 | 336 | 312 | ~1,244 exclusive + shared |
| Cytoplasm | 2,688 | 1,315 | 1,256 | ~5,199 exclusive + shared |

- **Cluster bootstrap on existing partitions (§4.1):** ample. 312 extracellular and 1,256
  cytoplasmic clusters in P4 is more than enough resampling units for a stable BCa interval.
- **Cluster-held-out development evaluation (§4.2):** ample. 1,244 extracellular-only and 5,199
  cytoplasm-only clusters support 5 folds with ≈ 249 and ≈ 1,040 clusters per fold.
- **Caveat — the minority class governs.** Extracellular has roughly one fifth as many
  clusters as cytoplasm at every role, is the more redundant class (§3.3a), and carries the one
  93-member block. Any per-layer difference in the sweep should be reported alongside the
  **extracellular cluster count** as well as the protein count: 575 proteins in 312 clusters in
  P4. Neither number is an effective sample size.
- **Not useful at 80 % identity:** that clustering is too fine to be a useful resampling unit
  (24 % of proteins in clusters ≥ 2, largest cluster 16) — it is close to resampling proteins.
  Use 30 % as the resampling unit and cite 50 % / 80 % only as sensitivity.

> **What cluster-based uncertainty does and does not buy.** Resampling clusters instead of
> proteins makes the interval account for **the grouping defined in §2** — it stops 93 members
> of one toxin family entering a resample as 93 free draws. It does **not** establish that the
> clusters are independent of one another. Two clusters may share a domain, or be related
> below 80 % coverage, and the cluster bootstrap will treat them as independent anyway. A
> cluster-resampled interval *accounts for the specified cluster grouping, without
> establishing independence between clusters*. It is not guaranteed to be wider than a
> protein-level interval.

### 4.4 What this audit does not license

- It does not show the partitions are family-disjoint. It shows no cross-partition relatedness
  is **detectable at 30 % identity with 80 % bidirectional coverage**.
- Shared-domain and partial homology remain **unmeasured**. A domain-level audit against
  InterPro/Pfam — the manifest already carries `n_pfam` / `n_interpro` for 9,745 proteins —
  would answer a different and complementary question, and is not done here.
- Nothing here justifies a layer sweep. [`LOCALISATION_EXPERIMENT_SPEC.md`](LOCALISATION_EXPERIMENT_SPEC.md)
  §6 makes a sweep conditional on the layer-8 experiment first yielding a candidate with a
  reproducible, specific effect, and recommends a single layer-24 check over a sweep if one
  cross-layer check is ever made. That condition is unchanged by this audit.

---

## 5. Deliverables and integrity

All under [`homology_audit/`](homology_audit/), version `v1`.

| File | Contents |
|---|---|
| [`build_input.py`](homology_audit/build_input.py) | Builds and verifies the frozen FASTA |
| [`run_clustering.sh`](homology_audit/run_clustering.sh) | Exact MMseqs2 commands, all three thresholds |
| [`analyse_clusters.py`](homology_audit/analyse_clusters.py) | All tables and integrity checks |
| [`v1/frozen_specification.json`](homology_audit/v1/frozen_specification.json) | Configuration, frozen pre-run |
| [`v1/input_report.json`](homology_audit/v1/input_report.json) | Input composition and hashes |
| [`v1/development_10139.fasta`](homology_audit/v1/development_10139.fasta) | The clustered sequences |
| [`v1/clusters/id{30,50,80}_cluster.tsv`](homology_audit/v1/clusters/) | Raw MMseqs2 output |
| [`v1/cluster_membership.tsv`](homology_audit/v1/cluster_membership.tsv) | **Per-protein table**: cluster id + size at all three thresholds, with role, class, kingdom, length, and every exclusion flag |
| [`v1/cluster_summary.tsv`](homology_audit/v1/cluster_summary.tsv) | 6 rows: {cohort, core} × 3 thresholds |
| [`v1/largest_clusters.tsv`](homology_audit/v1/largest_clusters.tsv) | Top 10 per threshold, with class and partition breakdown |
| [`v1/cross_partition_overlap.tsv`](homology_audit/v1/cross_partition_overlap.tsv) | All 6 partition pairs × 3 thresholds |
| [`v1/class_group_counts.tsv`](homology_audit/v1/class_group_counts.tsv) | Independent groups per role × class |
| [`v1/completed_sample_clusters.tsv`](homology_audit/v1/completed_sample_clusters.tsv) | Per-protein membership for the 400- and 40-protein samples |
| [`v1/integrity_checks.json`](homology_audit/v1/integrity_checks.json) | Checks below, plus all file hashes |
| [`v1/clustering_run.log`](homology_audit/v1/clustering_run.log) | Full MMseqs2 console output |

### Input hashes

| File | SHA-256 |
|---|---|
| `Swissprot_Train_Validation_dataset.csv` | `29a07b293fed2994a966b70bdcd6bacc59915b8b01fa200cb2b07d8db18384a2` |
| `discovery_analysis/v1/inclusion_manifest.tsv` | `783931349700b3082cac37edc5a2729e9b439cdb98bb7824a8c30e206ab2adb1` |
| `v1/development_10139.fasta` | `2fd9d8cc08c2745245ea9ab9dcbc213ed225e8a86ad87de8dba9fcf644497121` |
| `v1/frozen_specification.json` | `79b0985c7354f23a54b58b4901f23990ad59af10e3c3b296c38864573f208440` |

### Integrity checks — all pass

| Check | Result |
|---|---|
| Every protein assigned to exactly one cluster, all three thresholds | **Pass** — 10,139 each, no double assignment |
| Cluster membership set equals the manifest accession set | **Pass**, all three |
| Every representative is a member of its own cluster | **Pass**, all three |
| Cluster sizes sum to 10,139 | **Pass**, all three |
| All 10,139 sequences match their frozen `sha256` | **Pass** — 0 dropped, 0 substituted |
| Accessions unique in the manifest | **Pass** — 10,139 unique |
| Sequences unique | **Pass** — 10,139 distinct |
| Completed-sample accessions all inside the development manifest | **Pass** — 400/400 and 40/40 |
| Specification frozen before cluster output existed | **Pass** — written and hashed pre-run |

**Noted, not a failure:** 55 id50 clusters and 17 id80 clusters are not nested within a single
id30 cluster. Greedy cascaded clustering is not strictly hierarchical; the three thresholds are
three partitions, not a refinement chain. Recorded in
[`integrity_checks.json`](homology_audit/v1/integrity_checks.json) and in §2.1.
