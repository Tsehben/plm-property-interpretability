# UESolDS development curation

This additive analysis implements the 2026-11-08 review. Read `EDA_REPORT.md` for conclusions. The proposal, `FINDINGS.md`, `CHECKIN_1.md`, cohort builder, annotation pilot, baseline and embedding code were reviewed first. No original data, cached embeddings or old results are overwritten. The independent test file is never opened, including for hashing or homology searches.

## Reproduction

From the repository root, with Python 3.13 and MMseqs2 installed:

```sh
python3 -m venv /tmp/uesolds-venv
/tmp/uesolds-venv/bin/pip install -r Violet969-PLM_Sol-3589a95/curation/requirements-lock.txt
/tmp/uesolds-venv/bin/python Violet969-PLM_Sol-3589a95/curation/pipeline.py prepare
/tmp/uesolds-venv/bin/python Violet969-PLM_Sol-3589a95/curation/pipeline.py homology --mmseqs /path/to/mmseqs
/tmp/uesolds-venv/bin/python Violet969-PLM_Sol-3589a95/curation/pipeline.py finish
/tmp/uesolds-venv/bin/python Violet969-PLM_Sol-3589a95/curation/summarize.py
```

The optional report-rendering script reads the default results directory. Use the same `--out /path/to/new/results` on all three stages to reproduce separately. The run records input SHA-256 hashes, package versions, MMseqs2 binary version/hash and exact commands. This run downloaded the official macOS universal binary from https://mmseqs.com/latest/mmseqs-osx-universal.tar.gz; use the recorded revision for exact reproduction. MMseqs temporary databases live outside the repo. The original source parser is reused without calling its main routine. ESM tensors are loaded directly from trusted NPZ shards without loading torch or any model weights.

## Frozen policy and evaluation

1. Audit the complete published train and validation FASTAs, before exclusions. `raw_uid` is split plus actual zero-based FASTA position. `legacy_uid` preserves the old post-exclusion index used by caches. Join caches through verified ID, label, sequence hash and sequence; never silently reinterpret old keys.
2. Preserve the already documented `PDB_230301_28958` training exclusion without reopening test. Invalid syntax/empty sequences and ambiguous residues B/X/Z/J are excluded from cleaned cohorts, with reasons. U/O are recognized rare residues, retained. Descriptor hydropathy uses zero for unsupported residues; affected records are flagged. Conflicting identical sequences are quarantined; same-label duplicates retain the first occurrence within split. Training copies matching validation are removed. Original labels are never changed.
3. Source eligibility: SG centres with at least 100 examples of **each** label in the full eligible training pool. This threshold is a pragmatic support criterion, not optimized against AUROC. All their eligible validation records are retained, including small classes. No resampling, class weighting or prevalence balancing. Macro scores include only defined per-centre AUROCs; single-class scores remain undefined.
4. MMseqs2 direct validation-to-full-training search at >=30% identity, >=80% coverage of both sequences (`--cov-mode 0`), E <=0.001, sensitivity 7.5, actual identity (`--alignment-mode 3`), maximum 10,000 candidates. Also cluster all development sequences at the same thresholds using connected-component mode. Remove training records with direct hits, membership in a validation-containing cluster, or a validation-shared identifier. Validation is unchanged by these training exclusions. This is conservative, may remove whole related families, and does not establish absence of remote or partial-domain homology. Clustering and heuristic search are approximate, not exhaustive all-pairs alignment.
5. Preserve tags and variants. Flag terminal His runs (>=6 in first/last 30 residues), cleavage motifs ENLYFQ/LEVLFQ/LVPRGS and GS-rich runs (>=8 in first/last 40). Non-M starts and repeated IDs are **suspicions**, not evidence of truncation. Removing tags would alter an experimentally tested construct and invalidate cached embeddings. Do not do so without matched reference sequences and experimental metadata.
6. Strict-label inclusion requires construct-matched experimental outcomes. None are recoverable from the local release and annotations; the strict manifest is empty with an explicit reason. DNASU-like identifiers are candidates for future recovery, not verified evidence.
7. Refit source-only, composition plus log length, and frozen mean-pooled ESM-2 layer-8 logistic probes with fixed C=1, seed 42, no tuning, training-only scaling/encoding. Composition sees the same first 1,022 residues as ESM. All models use identical cached rows within each cohort. The full pool is audited, but ESM comparison training remains restricted to the existing 10,000-record sample and its retained subset. Full-training historical baselines remain in their existing directory. Report pooled, every centre and macro AUROC plus per-group AP; macro `n_eval` denotes the number of centres, not proteins. These are development measurements, not independent validation of preprocessing choices.

## Outputs

- `records.tsv`: complete original-record audit; `cache_mapping.tsv`: cache correspondence.
- `cohort_manifest.tsv`: all decisions, exclusions and cluster memberships; individual cohort manifests select rows. Manifests reference immutable original FASTA records rather than copying/rewriting their sequences.
- `source_counts.tsv`, `descriptor_summary.tsv`, `within_source_effects.tsv`, `artifact_summary.tsv`: EDA by split, source and label. Descriptor and artifact summaries have two header rows (statistic under feature).
- `annotation_coverage.tsv`: coverage over all records; `annotation_attempted_only.tsv`: success conditional on an attempted pilot annotation. Unqueried records must not be interpreted as biologically lacking an organism or family.
- `cross_split_hits.tsv`, `homology_pairs_annotated.tsv`, `diversity_cluster.tsv`: homology evidence and diversity proxy.
- `cohort_statistics.tsv`, cohort source counts, `baseline_results.tsv`, `predictions.tsv`: before/after comparison and reproducible scores.
- PNG figures and `summary.json`; input hashes and environment/MMseqs provenance.

The workflow asserts input preservation, cache alignment, unique embedding keys, finite features, model convergence and cleaned train–validation separation by exact sequence, detected homology, cluster and identifier. Search against the locked independent test is deliberately deferred to the final evaluation protocol.
