# Biological annotation: protocol and reproducibility

This is an additive annotation layer for **13,880 unchanged development proteins**
(9,880 train; 4,000 validation). It preserves heterogeneous A, including single-class
families, no-Pfam-hit proteins and missing annotations. Neither labels, centres nor
model performance determine annotation methods or parameters. Label/source strata
are used only for the operational pilot and descriptive summaries. No probe, SAE,
intervention, feature selection, ESM fine-tuning or independent-test access occurs.

## Inputs and evidence

`annotate_biology.py` verifies `results/cohort_manifest.tsv` against the completed
`pfam_annotation/annotated_cohort.tsv`, selecting the already-frozen `clean & cached`
records. It recovers full sequences from `esm2_probe/uesolds_pairs.tsv` by legacy
cache key, then checks length and SHA-256 for every original `raw_uid`. Do not join
only on protein names: identifiers can represent different constructs.

| Block | Evidence category | Important dependency |
|---|---|---|
| Composition, entropy, hydropathy, estimated charge/pI | Direct sequence-derived | Defined residue scales and ionization assumptions; not measurements |
| Aggregation | Biophysical/sequence-propensity predictor | Published AGGRESCAN scale and algorithm, independently reimplemented |
| Disorder | Pretrained-model prediction | metapredict V3, consensus-disorder and AlphaFold2-derived targets |
| TM topology | Pretrained-model prediction | DeepTMHMM2, ESM-2-650M final-layer representations |
| Secondary structure/accessibility, awaiting server output | Pretrained-model prediction | NetSurfP-3.0 uses ESM-1b |
| Signal peptides, awaiting server output | Pretrained-model prediction | SignalP-6.0 uses a protein language model and structured prediction |
| Pfam, taxonomy | External biological annotation | Profile models or exact matches to an archived reference sequence |
| Homology | Sequence-similarity-derived grouping | Existing MMseqs2 result; included with the external/grouping block |

There is **no new experimental measurement**. In particular, neither an exact
reference match nor a structural prediction verifies UESolDS solubility labels.
No pretrained solubility predictor is an explanatory feature.

## Deterministic biochemical descriptors

`biological_features.py` version 1.0 implements the following rules, also recorded
in `biological_annotation/protocol.json` and the feature dictionary.

- Canonical alphabet: `ACDEFGHIKLMNPQRSTVWY`. All 20 fractions divide by original
  sequence length. Entropy is Shannon composition entropy in bits. It is not a
  disorder predictor. Empty inputs and unsupported residues are explicit.
- Groups: positive `KRH`, negative `DE`, charged `DEKRH`, polar uncharged `STNQCY`,
  nonpolar `AGILMFPVW`, aromatic `FWY`, hydrophobic `AVILMFWY`, small `AGST`,
  proline `P`, cysteine `C`. The first three broad categories charged/polar
  uncharged/nonpolar partition the 20 residues; other groups intentionally overlap.
  “Positive” here is group membership, not a fully ionized histidine assumption.
- Hydropathy uses the Kyte–Doolittle scale already used in curation. Prespecified
  9- and 19-residue complete sliding windows produce means, minima/maxima and
  their coordinates. Earliest coordinates break tied extrema. Short sequences
  have missing window statistics, not zero. Hydrophobic patches are not TM calls.
- Charge uses Henderson–Hasselbalch estimates at pH **7.0 and 7.4**, including
  free termini. [Biopython 1.86's Bjellqvist implementation](https://github.com/biopython/biopython/blob/biopython-186/Bio/SeqUtils/IsoelectricPoint.py)
  supplies pKas: K 10.0, R 12.0, H 5.98; D 4.05, E 4.45, C 9.0, Y 10.0;
  default N-terminus 7.5 and C-terminus 3.55. N-terminal exceptions: A 7.59,
  M 7.0, S 6.93, P 8.36, T 6.82, V 7.44, E 7.7. C-terminal exceptions:
  D 4.55 and E 4.75. pI uses bisection over **[0,14]**, preventing default-bound
  clipping of very acidic proteins. These estimates ignore folded environments,
  disulfides, covalent modifications, terminal blocking and assay conditions.
- Local charge windows use **side-chain fractional charge per residue** without
  invented free termini at window boundaries; positive/negative group-density
  extrema are also retained. `(K+R−D−E)/length` remains a separately named proxy.
- For unsupported residues, whole-sequence hydropathy/charge/pI/entropy are
  missing; local windows touching unknown residues are missing. Known residue
  fractions remain observed fractions with an explicit unsupported-input status.
  No unknown residue is silently made neutral or assigned a standard amino acid.

## Aggregation: implementation and discrepancy check

The [2007 AGGRESCAN paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC1828741/)
and [published residue scale](https://ddd.uab.cat/pub/artpub/2007/112404/bmcbioinformatics_a2007m2v8n65addfin1.pdf)
define the implementation. Window widths are 5/7/9/11 for lengths ≤75/≤175/≤300/>300.
An APR is ≥5 consecutive residues with averaged propensity **>−0.02**, excluding
proline. Virtual N/C terminal residues use mean K/R and D/E propensities,
respectively; uncentred edge residues copy the nearest complete-window value.
Integer half-milliunit summation avoids floating-point threshold artefacts.
Outputs include mean/max propensity, APR count, residue fraction, longest region,
scores and intervals. Area metrics are not claimed.

**This is a published-algorithm reimplementation, not the original executable.**
The [supplementary SH3 example](https://ddd.uab.cat/pub/artpub/2007/112404/bmcbioinformatics_a2007m2v8n65addfin3.pdf)
matches interior scores, three/four APRs and −13.9/−11.9 normalized means for
wild type/N47A. Its terminal values instead match the reversed terminal convention.
Both conventions are implemented: prose-defined is primary, example-defined is
a recorded sensitivity. The discrepancy is unresolved, not hidden by changing
the primary rule. Inputs >2,000 residues are flagged as extrapolations beyond
the original server limit; the published >300-residue window rule still applies.
Inputs <5 or noncanonical are unsupported. Treat APR annotations as exploratory.

## Disorder

Use **metapredict package 3.1.1, network V3**, CPU, one thread, seed 42,
`torch.use_deterministic_algorithms(True)`, `normalized=True`, `round_values=False`.
The package internally reports `3.1.01`; installed distribution version is 3.1.1.
The active checkpoint `backend/networks/epoch099_val_loss1107.08.ckpt` has SHA-256
`3209c42d91a77dc5828618c7835ebc5dd912f56d43d0111d3c3c36138d9739a8`.
All checkpoint hashes, PyTorch version and runtime settings are recorded.

The [authors' documentation](https://metapredict.readthedocs.io/en/stable/faq.html)
describes a bidirectional sequence network trained with consensus-disorder and
AlphaFold2-derived targets. V3 uses actual AlphaFold2 pLDDT targets smoothed over
25 residues. It does not consume ESM embeddings at inference, but its outputs are
not independent experimental disorder ground truth. Metapredict is MIT licensed.
The [IUPred3 download](https://iupred3.elte.hu/download_new) requires academic
registration; it was not available locally. Metapredict is an explicitly different
validated method, not an implementation claimed to reproduce IUPred3 scores.

Disorder calls use **score ≥0.5**, consecutive segments of minimum one residue,
no gap closure or minimum-domain postprocessing. Mean score, fraction and longest
segment are descriptive outputs. Scores should not be read as calibrated physical
probabilities. The 71-record pilot spans split/source/label cells, length deciles,
extrema and terminal-motif flags; repeated pilot inference was identical.
Full-cohort inference is single-sequence, preventing batch-padding variation.

## TM topology, structure and targeting

[DeepTMHMM2](https://github.com/fteufel/DeepTMHMM2) is pinned to commit
`526f6c60f93c0868c7315413da8e669d12373cc3`, with its ten checkpoint hashes recorded
in `transmembrane_pilot.json`. A separate temporary environment installs the
upstream specified dependencies, including PyTorch 2.10.0, fair-esm 2.0.0 and
Lightning 2.6.1. Use CPU, four threads, batch size one, `--simplify-io`, default
upstream ensemble/thresholds. ESM backbone and companion file are outside the repo;
the backbone checksum is recorded. This predictor requires new annotation-time
ESM final-layer inference: existing mean-pooled layer-8 caches cannot substitute.
It does not train or modify ESM or run an ESM solubility probe.

The complete 71-record pilot succeeded, including long sequences. The user
authorized the full CPU run after its ~10–11-hour throughput estimate. The batch
runner checkpoints each 128 records and resumes completed shards after verifying
their input/settings/output hashes. Do not use MPS: upstream explicitly reports
numerical inconsistencies. Upstream chunks sequences beyond its context limit;
long-construct topology therefore remains exploratory. TM helix and beta-barrel
strand segments are distinct from interfacial/reentrant regions and generic
hydropathy. Signal/transit calls from DeepTMHMM2 remain separately named and do
not substitute for SignalP types or calibrated cleavage-site confidence.

[NetSurfP-3.0](https://services.healthtech.dtu.dk/services/NetSurfP-3.0/) and
[SignalP-6.0](https://services.healthtech.dtu.dk/services/SignalP-6.0/) require
institutional registration for local packages; the user has none and will upload
FASTA batches to their servers. No secondary-structure, accessibility or SignalP
values are invented. SignalP settings are derived only from exact species-consensus
lineage: Eukaryota → Eukarya, Bacteria/Archaea → Other; everything else unresolved.
Use SignalP **Slow**, **Long output**; retain full JSON/GFF3 and job metadata.
Server model/checkpoint hashes are not exposed: archive outputs/settings/date and
report this weaker reproducibility honestly. Do not claim pinned local checkpoints
for future server results.

## Existing Pfam, homology and exact taxonomy

Pfam **37.0**, PyHMMER **0.11.2**, profile-specific sequence/domain gathering
thresholds and all original MMseqs2 groups are reused; see [PFAM_METHODS](PFAM_METHODS.md).
`pfam_domains.tsv` preserves raw scores, accessions/versions, clans, alignment/HMM/
envelope coordinates and competing-clan overlap flags. Architectures order every
retained domain by alignment start, end, versioned accession, retaining repeats.
They are candidate ordered architectures, not resolved nonoverlapping domain maps.
No-hit is distinct from failed search. Coverage uses the union of aligned positions.
Class/family counts per cluster are cohort-level descriptive metadata: a multi-Pfam
protein itself can produce a cluster with multiple families.

`annotate_taxonomy.py` queries the [UniParc archive](https://www.uniprot.org/help/uniparc)
uniformly by full-sequence MD5, then requires **actual returned sequence equality**
and SHA-256 agreement. Queries are batched at 40 records, three workers, three
attempts, 60-second request timeout; errors are distinct from no match. Successful
raw responses are immutable gzip snapshots. Reference release is **UniProt 2026_04**
(04-November-2026), recorded in response headers; future-release downloads are
rejected by the script. Reproduction uses the archived responses because this
live API cannot promise historical query results indefinitely.

Only one archived common taxon at species rank or below, with verified lineage,
gets `assigned_exact_species_consensus`. Higher common ancestors, multiple taxa,
missing taxonomy IDs and unresolved lineages do not get a guessed organism.
Synthetic/artificial lineages are explicit engineered-taxon matches. UniParc IDs,
candidate taxa, lineage, identity/coverage and matched-sequence hashes are retained.
This is archive taxonomic consensus, potentially including historical cross-references,
not independent verification of a native construct. No high-identity approximate
mapping, tag removal, centre/host inference or legacy organism hint fills gaps.

## Coordinates, visibility and missingness

Intervals are **1-based inclusive**. Ragged NPZ offsets are **0-based half-open**,
with ordered IDs and sequence hashes. Disorder shards and aggregation profiles
retain every residue. `sequence_profiles.npz` stores exact ASCII residues; together
with fixed scales and code it losslessly reconstructs per-residue KD/charge values
and complete-window profiles without storing redundant arrays.

`full_*` summaries describe originals. `esm_*` descriptors use positions 1–1,022
with complete windows inside that prefix. Regional predictions are calculated on
the **full sequence**, then clipped: a qualifying region may intersect the visible
prefix by fewer than its original minimum length. No artificial terminal is added
to a clipped regional prediction. Prefix pI/net-charge estimates are explicitly
named `esm_hypothetical_fragment_*`; `esm_original_context_charge_*` instead sums
the original molecule's per-residue charges in the prefix. Full-context disorder
predictions can depend on unseen tail residues; clipping does not erase that
dependency. The preserved full predictions enable a later explicit comparison.

`ok`, `zero_features`, `no_hit`, `unsupported_input`, `failed`, `not_run_compute`,
and `unavailable` distinguish successful zeros from absent evidence. Numeric missing
values use `NA`; missing annotations never remove records or become biological zeros.
Original Pfam status and pilot-organism metadata are preserved separately.

## Reproduction and audit

Use a separate Python 3.13 environment with `requirements-biological.txt`.
Commands below run from the repository root; outputs default to the additive
`curation/biological_annotation/` directory. Use `--out` for isolated reproduction.

```sh
R=Violet969-PLM_Sol-3589a95/curation
python "$R/annotate_biology.py" prepare
python -m unittest discover -s "$R" -p test_biological_annotation.py -v
python "$R/annotate_biology.py" disorder-pilot
python "$R/annotate_biology.py" disorder
python "$R/annotate_taxonomy.py" --pilot
python "$R/annotate_taxonomy.py"
python "$R/annotate_biology.py" build
python "$R/summarize_biology.py"
python "$R/render_biological_report.py"
python "$R/export_annotation_fastas.py" --destination /tmp/uesolds-server-upload
```

For taxonomy replay, copy `taxonomy_responses/` into the new output directory first;
matching URLs reuse snapshots without network access. Copy no cohort outputs from
other experiments. Sequence/manifest/Pfam/homology hashes must match the frozen inputs.
Large HMMs, model weights, environments and temporary FASTAs remain outside git.

For the authorized TM job, install the pinned GitHub source in a separate environment
and download its ESM files under `TORCH_HOME=/tmp/uesolds-tm-models`. Run
`pilot_transmembrane.py`, then `annotate_transmembrane.py` (or resume the latter).
Options expose executable/checkpoint/work directories, CPU/CUDA and thread count.
The current job and log location are recorded in `transmembrane_job.json`;
`transmembrane_progress.json` reports completed-shard coverage. A watcher refreshes
the unified table, QC/EDA and report when all 13,880 records finish. The computer
must remain running; after interruption, rerun the batch command to resume.

QC checks every record's identity, hash, original split/label/cluster and retained
Pfam metadata, all domain/region coordinates, score preservation and scalar/profile
agreement. Unit tests cover published external aggregation fixtures, the terminal
discrepancy, unknown/short inputs, pI boundaries, termini, clipping and proline breaks.
EDA uses fixed concept representatives, all 20 AA fractions, source summaries and
families with ≥5 records per split for variation summaries. Descriptive label
contrasts require ≥3 records of each label in the reported stratum. Multidomain
memberships overlap; strata and contrasts are not independent hypothesis tests.
No association strength selects an annotation or changes primary A.
