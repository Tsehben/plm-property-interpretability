# Experimental solubility dataset comparison

Start with [DATASET_RECOMMENDATION.md](DATASET_RECOMMENDATION.md). The audit covers source validity, exact outcome counts, biological diversity, within-family feasibility and cross-release overlap. It contains no trained models or new ESM embeddings.

Reports:

- [EXPERIMENTAL_LABEL_AUDIT.md](EXPERIMENTAL_LABEL_AUDIT.md)
- [DATASET_COUNTS.md](DATASET_COUNTS.md)
- [FAMILY_DIVERSITY.md](FAMILY_DIVERSITY.md)
- [DATASET_OVERLAP.md](DATASET_OVERLAP.md)
- [DATASET_RECOMMENDATION.md](DATASET_RECOMMENDATION.md)

## Inputs and preservation

`source_manifest.json` records download URLs, Git commit, retrieval date, exact byte counts and SHA256. Public raw data are retained under `raw/`; papers/supplements are local research reference copies. Git's new directory-specific ignore file excludes those documents and temporary outputs from accidental commits. `fetch_sources.py` downloads missing sources and rejects checksum changes. The SoDoPe original reconstruction uses author nucleotide exports and class mappings, rather than loading untrusted pickles. `PSI_reconstruction_exclusions.tsv` records the cleanup stages. Author benchmark FASTAs are retained as reference downloads, not silently substituted for the 12,216-protein reconstruction.

Existing UESolDS inputs are explicitly limited to:

- `curation/results/cohort_manifest.tsv` (train/validation only)
- `esm2_probe/uesolds_pairs.tsv` (existing cached development sequences)
- completed `curation/pfam_annotation/` data and Pfam clan names
- previous compact DNASU recovery evidence, referenced without repeating web requests

Only 13,880 `clean & cached` proteins enter the homology/annotation union. The larger train/validation manifest is used only for a separate exact-hash comparison. No independent UESolDS test path is used. Original sequences, labels, split fields, source assignments, cluster IDs and annotation outputs remain unchanged. New `audit_cluster` and `audit_fold` columns are explicitly separate from original roles.

## Reproduce

Python 3.13.7 was used; package versions are in `requirements.txt`. Local environment: `/tmp/uesolds-venv/bin/python`. From the repository root:

```sh
python dataset_comparison/fetch_sources.py
python dataset_comparison/audit.py prepare
python dataset_comparison/audit.py annotate
python dataset_comparison/audit.py homology
python dataset_comparison/analyze.py
python dataset_comparison/supplement.py
MPLCONFIGDIR=/tmp/uesolds-matplotlib python dataset_comparison/plot_results.py
python dataset_comparison/render_reports.py
python -m unittest discover -s dataset_comparison -p 'test_audit.py' -v
```

The existing Pfam reference is expected at `/tmp/uesolds-pfam37/`; MMseqs2 at `/tmp/mmseqs/bin/mmseqs`. These paths are constants in `audit.py` and can be edited for another machine without changing thresholds. The union FASTA and temporary MMseqs databases are generated outside the checkout in `/tmp/uesolds-dataset-comparison/`. Original downloaded sequence files are retained as audit evidence; no ESM-input sequences are overwritten.

### Pfam protocol

Reuse completed full-sequence hits by exact SHA256 for 13,880 sequences. Search the other 30,277 unique sequences against all 21,979 Pfam 37.0 profiles with PyHMMER 0.11.2; sequence and domain GA cutoffs must both pass. Same implementation as `curation/annotate_pfam.py`: seed 42, six query-parallel CPUs, F1=.02/F2=.001/F3=.00001, bias filter and null2 enabled. New-search Z=30,277; reused E-values retain their original Z=13,880 provenance. Selection uses profile GA scores rather than a cohort-dependent E-value threshold. Exact reference SHA256 and settings are in `pfam_run.json`.

Reference URL: https://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam37.0/

Pfam-A.hmm.gz SHA256: `66f57012fe5825c82f2e31447d7a8026367e1ca9e5e91b4b716c594c268db69b`.

Coordinates are 1-based inclusive. Full domains, scores, cutoffs, HMM coverage, clan metadata, search origin and competing-overlap flags are saved. Three unique sequences contain X; HMMER receives X unchanged, with no amino-acid imputation. Physicochemical summaries for unknown-residue sequences remain missing. All sequences are accounted for in `annotation_status.tsv.gz`; a searched/no-hit record is not an annotation failure.

### Homology protocol

Direct all-vs-all MMseqs2 `easy-search`, not just cluster-representative matching: ≥30% identity, ≥80% bilateral coverage (`--cov-mode 0`), E≤.001, sensitivity7.5, alignment-mode3, maximum10,000 targets, four threads. This uses the existing UESolDS threshold policy and binary; the exact version/hash/command are in `homology_run.json`. Search hit counts are checked against the maximum-target cap. Union connected components use the lexicographically smallest SHA256 as identifier. Direct ≥90%/≥95% overlap is also summarized.

The audit search can join original UESolDS groups through newly included sequences. It does not revise the original cluster column. All detected homology edges stay within the proposed component-hash folds. Remote homology, short shared domains and approximate-search false negatives remain possible; Pfam memberships are never substituted for these groups.

### Outcomes and statistics

Canonical `records.tsv.gz` retains version-specific IDs, exact original sequence and label, raw E/S, original split/partition/source/cluster fields, and SHA256. `policy_membership.tsv.gz` maps every proposed outcome policy to those records, with no mutation of published labels. `policy_groups.tsv.gz` adds audit groups/folds. `additional_eligibility.tsv` records conservative conflict/deduplication and external-evaluation views. Threshold variants were chosen for interpretation, not model performance.

Family support counts records once per family; eligible-protein totals are unions. Duplicate records may increase raw protein counts, so group thresholds and consistent-deduplicated sensitivities are reported separately. Class-exclusive group thresholds are stricter than per-class group-presence thresholds. Representative distributions use longest domain span, bit score, accession and coordinate; they sum to 100%, but intentionally discard architecture information. Main plots instead show all-membership protein prevalence.

Direct biochemical descriptors are length, all standard amino-acid fractions, `(K+R-D-E)/length` (a proxy, not pH-specific charge), mean Kyte–Doolittle hydropathy, terminal His-run flags and >1,022-residue flags. None are fitted to labels.

## Key supporting files

| File | Purpose |
|---|---|
| `results/dataset_counts.tsv` | Exact counts and explicit policy denominators |
| `results/NESG_E_by_S.tsv` | Complete raw-score contingency table, including fractional E |
| `results/NESG_exclusion_reasons.tsv` | Disjoint exclusion reasons by score policy |
| `results/NESG_score_anomalies.tsv` | Unresolved fractional/zero-expression observations |
| `results/PSI_construct_mapping.tsv` | Original ↔ curated sequence offsets and label agreement |
| `results/NetSolP_NESG_score_checks.tsv` | Released-score validation, including three rule exceptions |
| `results/family_diversity.tsv` | All policies, GA and annotation-quality sensitivities |
| `results/family_support.tsv.gz` | Family-specific protein/group/evaluation support |
| `results/top20_named_families_clans.tsv` | Actual named biological categories by class |
| `results/five_fold_feasibility.tsv` | Fixed-split sensitivity; no best-fold selection |
| `results/PSI_family_vector_support.tsv` | Controlled same-vector contrasts |
| `results/dataset_overlap.tsv` | Exact and direct homology overlap |
| `results/duplicate_records.tsv` | Within-release duplicate/conflicting binary outcomes |
| `results/exact_conflicts.tsv` | Between-release exact-sequence binary disagreements |
| `results/biochemical_summary.tsv` | Sequence-descriptor distributions by label |

Reports are descriptive. New diagnostic partitions are not historically untouched holdouts. No observed AUROC enters the recommendation.
