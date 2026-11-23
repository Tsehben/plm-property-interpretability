# Discovery analysis inputs and controlled-comparison feasibility — view `v1`

Saved data only. No network requests, no ESM or SAE extraction, no model fitting, no
intervention runs, nothing committed. Source data is read and never written: the frozen
manifest, the raw annotation cache, `annotation_readiness/out/`, `out_reconciled/`, and the
intervention-validation results are all unchanged.

Everything below comes from **discovery data (P1+P2) only**, except role-level integrity
counts, which are permitted. No P3 or P4 feature association, model fit or performance
figure was inspected, and no group or threshold was chosen on a model outcome.

---

## 1. Analysis view — paths and reproducible rules

| artefact | path |
|---|---|
| builder (reproducible, idempotent) | `discovery_analysis/build_analysis_view.py` |
| inclusion manifest, all 10,139 proteins with exclusion reasons | `discovery_analysis/v1/inclusion_manifest.tsv` |
| discovery analysis view (P1+P2, core sample) | `discovery_analysis/v1/analysis_view_discovery.tsv` |
| sample counts by role and class | `discovery_analysis/v1/sample_counts.tsv` |
| discovery counts by class × signal status × evidence | `discovery_analysis/v1/discovery_signal_evidence_counts.tsv` |
| per-group discovery counts (Pfam + InterPro) | `discovery_analysis/v1/discovery_group_counts.tsv` |
| eligibility at 5 / 10 / 20 per class per group | `discovery_analysis/v1/discovery_group_eligibility.tsv` |
| domain-annotation applicability | `discovery_analysis/v1/domain_applicability.tsv` |
| integrity checks and input hashes | `discovery_analysis/v1/integrity_checks.json` |

### Inclusion rules

Applied in fixed precedence, so each excluded protein is attributed to exactly one reason
and overlaps are reported rather than summed.

| id | rule | scope |
|---|---|---|
| **E1** | `engineering_protein == 1` | excluded from every selection and evaluation sample |
| **E2** | `p0_base_collision != 'none'` | excluded from every selection and evaluation sample |
| **E3** | `coordinates_transferable == 0` | excluded from **annotation-conditioned** samples only |

Samples:

- **`cohort`** — all 10,139 development proteins. Never altered; E1/E2/E3 proteins are
  retained here with their flags.
- **`core`** — cohort minus E1 ∪ E2. Label-conditioned analyses (feature discovery against
  the frozen DeepLoc binary label) use this.
- **`annotated`** — core minus E3. The only sample in which a signal-peptide annotation may
  be conditioned on, because it is the only one where the annotation is known to apply to
  the frozen extraction-input sequence.

Fixed rules enforced in code, each as an assertion that fails the build:

- P1+P2 = discovery, P3 = tuning, P4 = evaluation, from the manifest partition column.
- **No P0 row enters any sample.** The partition column is read only to assert this.
- **Targets come from the frozen DeepLoc binary label and nothing else.**
  `subcellular_location_text` is not in the view's field list and is asserted absent from
  every output row. It is never a concept input and never a source of a target.
- `not_annotated` is carried through as its own state and is never merged with absence.

---

## 2. Before / after counts and integrity

### Exclusion accounting, without double-counting

| criterion | n (cohort-wide) |
|---|---|
| E1 engineering | 8 |
| E2 p0_base_collision | 27 |
| E3 coordinates not transferable | 381 |
| E1 ∩ E2 | 0 |
| E1 ∩ E3 | 1 (`O60462-6`) |
| E2 ∩ E3 | 21 |
| E1 ∩ E2 ∩ E3 | 0 |
| **union E1 ∪ E2** | **35** (naive sum 35) |
| **union E1 ∪ E2 ∪ E3** | **394** (naive sum 416 — 22 would be double-counted) |

### Sample sizes

| sample | role | n | Extracellular | Cytoplasm | prevalence |
|---|---|---|---|---|---|
| cohort | discovery (P1+P2) | 5201 | 1339 | 3862 | 0.2575 |
| cohort | tuning (P3) | 2556 | 584 | 1972 | 0.2285 |
| cohort | evaluation (P4) | 2382 | 575 | 1807 | 0.2414 |
| **core** | **discovery** | **5190** | **1336** | **3854** | **0.2574** |
| core | tuning | 2548 | 584 | 1964 | 0.2292 |
| core | evaluation | 2366 | 571 | 1795 | 0.2413 |
| **annotated** | **discovery** | **5009** | **1305** | **3704** | **0.2605** |
| annotated | tuning | 2453 | 573 | 1880 | 0.2336 |
| annotated | evaluation | 2283 | 554 | 1729 | 0.2427 |

Discovery loses 11 proteins to E1 ∪ E2 (3 engineering, 8 collisions) and a
further 181 to E3 when annotations are conditioned on. Class prevalence moves by less than
0.004 at every step, so no exclusion materially reweights the classes.

### Integrity checks

| check | result |
|---|---|
| rows in inclusion manifest | 10,139 = expected |
| unique accessions | 10,139 |
| P0 rows in any sample | **0** |
| labels and partitions vs frozen manifest | identical, asserted row by row |
| target source | frozen DeepLoc binary label only |
| `subcellular_location_text` present in any view output | **no** (asserted) |
| core sample = cohort − \|E1 ∪ E2\| | 10,104 = 10,139 − 35 ✓ |
| annotated sample = cohort − \|E1 ∪ E2 ∪ E3\| | 9,745 = 10,139 − 394 ✓ |
| network requests | 0 |
| source files modified | none |

Input hashes are recorded in `v1/integrity_checks.json`.

---

## 3. Discovery-only feasibility

### 3.1 Usable counts by class, signal status and evidence

Discovery core, n = 5,190. `n_conditioned` is the count surviving E3, i.e. the count usable
when the annotation must apply to the frozen sequence.

| class | signal status | evidence | n (core) | n (annotation-conditioned) |
|---|---|---|---|---|
| Extracellular | annotated_present | sequence_analysis_predictor | 954 | 954 |
| Extracellular | annotated_present | experimental | 236 | 236 |
| Extracellular | annotated_present | by_similarity | 33 | 33 |
| Extracellular | annotated_present | unspecified | 14 | 14 |
| Extracellular | annotated_present | curator_inference | 14 | 14 |
| Extracellular | annotated_present | author_statement | 1 | 1 |
| Extracellular | **not_annotated** | — | 53 | 53 |
| Extracellular | **unverifiable** | — | 31 | **0** |
| Cytoplasm | annotated_present | sequence_analysis_predictor | 53 | 53 |
| Cytoplasm | annotated_present | experimental | 4 | 4 |
| Cytoplasm | annotated_present | by_similarity | 2 | 2 |
| Cytoplasm | **not_annotated** | — | 3645 | 3645 |
| Cytoplasm | **unverifiable** | — | 150 | **0** |

`not_annotated` means a matching record was retrieved and carries no signal-peptide
annotation. It is **not** biological absence, and it is not pooled with `unverifiable`
anywhere above.

The binding constraint for any evidence-stratified analysis is visible here: the
low-dependency subset is **240 experimentally evidenced annotated proteins in discovery**
(236 Extracellular + 4 Cytoplasm), against 1,007 with computational evidence.

### 3.2 Group counts, multiple memberships retained

A protein appears in every Pfam and InterPro group it carries; group sizes therefore sum to
more than the number of proteins, and no protein is assigned one family.

| | Pfam | InterPro |
|---|---|---|
| groups present in discovery core | 2,933 | 6,175 |
| of which singletons (1 protein) | 1,519 | 3,063 |
| of which mixed-class (both classes present) | **103** | **272** |

Per-group counts, with the P1/P2 split of each class and the applicability breakdown, are in
`v1/discovery_group_counts.tsv`.

**InterPro entry type is not available from saved data.** The UniProt cross-reference we
retrieved carries only the InterPro accession and entry name; it does not carry the entry
*type*, so "Family", "Domain" and "Homologous superfamily" cannot be distinguished here.
Entry names follow conventions (`_sf`, `_dom`) but a naming convention is not a type field,
and typing groups from it would be inventing data. Every InterPro group below is therefore
**untyped**, and several are homologous superfamilies rather than families. Obtaining the
types requires an InterPro fetch, which is outside this step's scope.

### 3.3 Mixed-class groups and their P1/P2 distribution

The largest mixed-class groups with ≥10 per class, discovery core:

| db | group | name | n | Ext (P1/P2) | Cyt (P1/P2) | applies to frozen seq |
|---|---|---|---|---|---|---|
| InterPro | IPR013320 | ConA-like_dom_sf | 92 | 59 (41/18) | 33 (14/19) | 86 |
| InterPro | IPR013783 | Ig-like_fold | 80 | 34 (17/17) | 46 (23/23) | 65 |
| InterPro | IPR029058 | AB_hydrolase_fold | 76 | 42 (32/10) | 34 (16/18) | 74 |
| InterPro | IPR018247 | EF_Hand_1_Ca_BS | 74 | 10 (4/6) | 64 (44/20) | 69 |
| InterPro | IPR032675 | LRR_dom_sf | 65 | 10 (4/6) | 55 (20/35) | 63 |
| InterPro | IPR036179 | Ig-like_dom_sf | 36 | 21 (8/13) | 15 (9/6) | 23 |
| InterPro | IPR014756 | Ig_E-set | 33 | 11 (7/4) | 22 (10/12) | 33 |
| InterPro | IPR003599 | Ig_sub | 32 | 18 (7/11) | 14 (9/5) | 19 |
| InterPro | IPR007110 | Ig-like_dom | 32 | 17 (6/11) | 15 (9/6) | 21 |
| InterPro | IPR003598 | Ig_sub2 | 28 | 14 (5/9) | 14 (8/6) | 17 |
| InterPro | IPR016166 | FAD-bd_PCMH | 28 | 13 (9/4) | 15 (**0**/15) | 28 |
| InterPro | IPR016169 | FAD-bd_PCMH_sub2 | 28 | 13 (9/4) | 15 (**0**/15) | 28 |
| InterPro | IPR036318 | FAD-bd_PCMH-like_sf | 28 | 13 (9/4) | 15 (**0**/15) | 28 |

All seven eligible **Pfam** groups at ≥5 per class:

| group | name | n | Ext (P1/P2) | Cyt (P1/P2) | applies to frozen seq |
|---|---|---|---|---|---|
| PF00622 | SPRY | 24 | 5 (5/**0**) | 19 (9/10) | 24 |
| PF07679 | I-set | 19 | 5 (3/2) | 14 (9/5) | 16 |
| PF00041 | fn3 | 16 | 7 (7/**0**) | 9 (6/3) | 15 |
| PF00061 | Lipocalin | 16 | 8 (2/6) | 8 (**0**/8) | 16 |
| PF13855 | LRR_8 | 16 | 8 (3/5) | 8 (4/4) | 16 |
| PF13765 | PRY | 15 | 6 (6/**0**) | 9 (3/6) | 15 |
| PF00149 | Metallophos | 14 | 7 (**0**/7) | 7 (3/4) | 14 |

Two structural facts matter more than the group sizes. **Groups are nested and redundant** —
the three FAD-bd_PCMH entries are the same 28 proteins, and the Ig entries overlap heavily —
so a count of eligible groups is not a count of independent comparisons. And **several
groups have an empty partition cell**, which makes a within-group P1-versus-P2 check
impossible for them regardless of total size.

### 3.4 Eligibility, descriptive, at 5 / 10 / 20 per class per group

No threshold is selected here, and none was chosen on a model outcome.

| db | ≥ n/class | eligible groups | distinct member sets | nested in another | ≥2 per class in **both** P1 and P2 | with an empty partition cell | distinct proteins covered (Ext / Cyt) |
|---|---|---|---|---|---|---|---|
| Pfam | 5 | **7** | 7 | 0 | **2** | 5 | 93 (38 / 55) |
| Pfam | 10 | **0** | 0 | 0 | 0 | 0 | 0 |
| Pfam | 20 | **0** | 0 | 0 | 0 | 0 | 0 |
| InterPro | 5 | **38** | 35 | 17 | **21** | 16 | 704 (333 / 371) |
| InterPro | 10 | **13** | 11 | 4 | **10** | 3 | 423 (170 / 253) |
| InterPro | 20 | **3** | 3 | 0 | **3** | 0 | 246 (133 / 113) |

Coverage is of the 5,190-protein discovery core, deduplicated across overlapping groups:
13.6% at InterPro ≥5, 8.1% at ≥10, 4.7% at ≥20; 1.8% at Pfam ≥5.

### 3.5 Which annotations apply to the frozen sequence

Pfam and InterPro cross-references in the retrieved records are **entry-level**: they are not
residue-mapped, and for an isoform accession they describe the parent entry.

| applies to | n (discovery core) | share |
|---|---|---|
| entry matches the frozen sequence | 5,008 | 96.5% |
| parent-entry annotation, **unresolved isoform applicability** | 115 | 2.2% |
| parent-entry annotation, retrieved sequence differs from the frozen input | 67 | 1.3% |

So 182 discovery proteins (3.5%) carry domain annotations whose applicability to the exact
sequence we extract from is unestablished. Per-group counts of this split are in
`discovery_group_counts.tsv` (`n_applies_to_frozen`, `n_parent_entry_isoform`,
`n_parent_entry_seq_differs`) — note that for three of the seven eligible Pfam groups, and
for Ig-like_fold, a non-trivial minority of members fall in this category.

Signal-peptide annotations are handled more strictly: the `annotated` sample already removes
every protein whose coordinates do not transfer, so all 1,311 annotation-conditioned
discovery signal-peptide annotations apply to the frozen sequence by construction.

### 3.6 Homology clusters

**None exist.** No cluster identifier column appears in any saved table, and no clustering
artefact (MMseqs2, CD-HIT or otherwise) exists anywhere in the project tree. Consequences,
stated rather than assumed away:

- Independent homology-group counts **cannot be reported** for this cohort.
- Family-controlled generalisation **remains unverified**.
- Pfam or InterPro co-membership is not a homology cluster, and the exclusions in §1 do not
  establish homology independence. Preserving DeepLoc's published partitions does not make
  the partitions family-disjoint.

---

## 4. Feasibility table

| comparison | status | basis |
|---|---|---|
| Candidate discovery against the frozen binary label on P1+P2 | **supported** | 5,190 core proteins, 1,336 / 3,854, no annotation needed |
| P1-versus-P2 stability of a discovery statistic, cohort-wide | **supported** | 2,484 P1 and 2,706 P2 core proteins; both classes well populated |
| Signal-peptide-conditioned description of candidates | **supported** | 5,009 annotation-conditioned discovery proteins, coordinates verified |
| The same, restricted to experimental evidence | **supported but small** | 240 proteins (236 Ext / 4 Cyt). Enough to describe, not to stratify within groups |
| "Is the feature more than a signal-peptide detector?" via `not_annotated` Extracellular proteins | **weakly supported** | 53 Extracellular proteins only, and `not_annotated` is not absence — descriptive at best |
| Within-**Pfam**-group comparison, both classes, ≥10/class | **not supported** | 0 eligible groups |
| Within-**Pfam**-group comparison, ≥5/class | **marginal** | 7 groups, 93 proteins, only 2 with both partitions populated |
| Within-**InterPro**-group comparison, ≥10/class | **supported for a bounded set** | 13 groups (11 distinct member sets, 4 nested), 423 proteins, 10 with both partitions populated |
| Within-group comparison typed to *families* specifically | **needs more data** | InterPro entry type is not in the saved records; requires an InterPro fetch |
| Family-controlled generalisation claim | **not supported** | no clustering exists |
| Independent homology-group counts | **not supported** | no clustering exists |
| Residual-homology quantification between partitions | **not supported** | requires clustering |
| Domain-conditioned analysis on the exact extracted sequence for all members | **supported with a stated caveat** | 3.5% of discovery members carry parent-entry annotations of unresolved applicability; the column is in the view |

Two items would move several rows: an MMseqs2 `easy-cluster` run at 30% identity over the
frozen sequences (bounded, local, no network), and an InterPro type lookup (needs network).
Neither is performed here.

---

## 5. One recommended bounded candidate-discovery run

Proposed, not executed. Discovery data only; nothing here consults P3 or P4.

| element | proposal |
|---|---|
| **Sample** | `analysis_view_discovery.tsv`, n = 5,190 (core). Label = frozen DeepLoc binary target |
| **Activations** | layer 8, InterProt SAE, top-k k = 128, as already validated |
| **Pooling** | five summaries per feature per protein: `mean_whole`, `mean_nterm`, `frac_active_nterm`, `frac_active_whole`, `max_act_nterm`. N-terminal window = first 30 residues after the start token, fixed before the run |
| **Eligibility** | feature active in ≥5% of discovery proteins, measured on our own activations, never on the SAE's training-time dead list |
| **Selection statistic** | \|AUROC − 0.5\| of the pooled summary against the binary label, computed on P1+P2 together. One statistic, fixed in advance |
| **Candidate limit** | top **20** features by the statistic, across all five poolings combined (not 20 per pooling) |
| **Stability** | recompute the statistic **independently on P1 and on P2**, rank within each, and require a candidate to appear in the top 20 of **both**. Report the overlap count; the pre-registered bar is ≥10 of 20. Report Spearman rank correlation between the P1 and P2 statistics over eligible features as a descriptive diagnostic |
| **Frozen before any tuning or evaluation** | candidate list, eligibility set, pooling definitions and window, written to disk and SHA-256 hashed. Control families L and M (§4.4a/§4.4b of the specification) are drawn only after this hash exists |
| **Cost** | one extraction pass over 5,190 proteins. At the measured 4.11 seq/s this is ≈21 minutes, plus SAE encoding; no intervention passes |
| **What it does not do** | it does not select strengths, does not touch P3 or P4, does not use annotations to select anything, and does not establish that any candidate encodes localisation |

Annotations enter only **after** the candidate list is hashed, and only descriptively: for
each candidate, report its signal-status and evidence-category breakdown on the 5,009
annotation-conditioned discovery proteins, split by evidence category, with `not_annotated`
kept distinct from absence.

---

## 6. Control specification, as revised

The specification now distinguishes two control families rather than one (§4.4a, §4.4b,
§5.3 of `LOCALISATION_EXPERIMENT_SPEC.md`):

- **Family L — low-association controls.** Activation frequency and mean active magnitude
  matched to ±20%, *plus* a constraint that discovery \|AUROC − 0.5\| ≤ 0.02. Because that
  third constraint selects on a measured statistic, family L is a conditioned subset of the
  dictionary and **is not an unbiased random-feature null**. It is never described as "a
  random feature".
- **Family M — activation/norm-matched controls, no association filter.** The same two
  matching constraints, no AUROC condition, several controls per candidate from a recorded
  seed.

They answer different questions and are reported in separate rows, never pooled. A candidate
that beats one but not the other is an informative result, reported as such. **No controls
are drawn and no candidate SAE features are selected at this stage.**
