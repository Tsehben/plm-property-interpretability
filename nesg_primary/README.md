# NESG dataset validation

Start with the [executed notebook](NESG_Dataset_Validation.ipynb) and [validation summary](NESG_VALIDATION_SUMMARY.md). The notebook contains the analysis, tables and ten embedded figures; saved outputs are visible without a running kernel.

The strict cohort contains **3,762 proteins: 2,614 soluble and 1,148 insoluble**. It supports heterogeneous ESM analysis, with limited support for controlled within-family comparisons.

## ESM analyses

See the [ESM workflow index](esm_layer_probing/README.md) for the executed extraction, solubility and uncontrolled-property notebooks, and the prepared within-property evaluation notebook.

## Results

- [Validation summary](NESG_VALIDATION_SUMMARY.md): phenotype, cohort, split design and limitations.
- [Mixed-label family tables](SPLIT_MIXED_LABEL_FAMILIES.md): complete training, validation and test tables, using ≥1 protein per class independently in each split.
- [Family-distribution figures](figures/split_family_distributions/): sizes, counts and soluble/insoluble proportions for those families.
- [Test-label exposure record](TEST_LABEL_EXPOSURE.md): the subsequent authorized inspection of held-out labels by family.

## Reproduce

Use Python 3.13.7 and the versions in `requirements.txt`. The original environment is recorded in `environment_freeze.txt` and `results/execution_provenance.json`. From the repository root:

```sh
python nesg_primary/execute_notebook.py
python nesg_primary/write_summary.py
python -m unittest discover -s nesg_primary -p 'test_validation.py' -v
python nesg_primary/plot_split_family_distributions.py
```

The runner uses the invoking Python, saves executed outputs and fails on cell errors or missing figures. It needs local Jupyter sockets but no network access or new annotations.

Edit notebook content in `build_notebook.py`. Running that builder clears saved outputs; execute the notebook afterward. `analysis_helpers.py` handles data loading, family counts, hashes and split finalization. The analysis itself remains in notebook cells.

`report_split_mixed_families.py` regenerates the full family tables and **reopens reserved-test labels**. It is a separate inspection operation, not part of routine notebook reproduction. The plotting script uses the existing tables only.

## Frozen data and test access

The initial run checked the cohort against the pinned NESG archive, created splits once and wrote `results/test_lock.json`. Ordinary reruns display frozen aggregates in Sections 1–4 and analyze training/validation records in Sections 5–6. A kernel audit hook blocks restricted test records and upstream project inputs; `results/notebook_project_file_access.json` records access. This is an access guard, not encryption.

**Test labels were subsequently inspected by family at the user's request.** The unchanged lock records the original freeze; it does not mean the test remains uninspected. See [TEST_LABEL_EXPOSURE.md](TEST_LABEL_EXPOSURE.md). Do not use those reports to revise the cohort, splits, model selection or hypotheses. The UESolDS independent test was not accessed.

Keep the lock and frozen manifests. Reproduce the initial full-cohort audit in a separate copy if needed, with prior exposure disclosed.

## Definitions

- **Cohort:** 9,703 recovered records → 9,272 singleton-sequence records → 3,762 with integer E=4/5 and S=0/4/5. All repeated observations remain in the exclusion ledger. The historical 9,644-record membership is unresolved.
- **Pfam:** release 37.0, PyHMMER 0.11.2, sequence and domain gathering thresholds. The conservative sensitivity removes overlapping different-family hits within the same clan. Proteins with no hits remain in the cohort.
- **Homology:** existing MMseqs2 connected components at ≥30% identity, ≥80% coverage of both sequences and E≤0.001. Other audited proteins can bridge NESG components; no reclustering is performed.
- **Split:** existing audit fold 0 is test, 1 validation, and 2–4 training. Training-only CV uses SHA256(`nesg-primary-inner-v1:` + group ID), first 16 hex digits modulo 5. No class balancing or split search.
- **Family-held-out analysis:** a separate training-only design joins proteins sharing any GA family or homology group. It excludes shared Pfam accessions across folds, but cannot exclude every remote relationship or shared clan.

## Data files

| Path under `results/` | Contents |
|---|---|
| `development_records.tsv.gz`, `development_domains.tsv.gz` | Training/validation records and Pfam hits |
| `training_manifest.tsv`, `validation_manifest.tsv` | Frozen roles and labels |
| `training_cv_manifest.tsv` | Training-only CV assignments |
| `split_assignments.tsv` | All record/hash/group assignments, without test phenotypes |
| `reserved_test/` | Restricted test records and full-cohort manifest |
| `restricted_prelock/raw_record_decisions.tsv.gz` | Restricted original observations and exclusion reasons |
| `within_family_feasibility.tsv`, `controlled_family_candidates.tsv` | Training-selected family support and validation counts |
| `development_family_support.tsv` | Development-family counts under both annotation policies |
| `training_family_heldout_feasibility_manifest.tsv` | Separate training-only family folds |
| `split_mixed_label_families/` | Later split-local tables, including inspected test labels |
| `integrity_checks.json`, `final_validation.json` | Creation-time and development-only checks |
| `test_lock.json`, `input_checksums.json` | Frozen artifact and input hashes |

Pfam and MMseqs2 run metadata remain in `dataset_comparison/results/`. Notebook plots are also saved under `figures/`.

## AGGRESCAN server automation

`run_aggrescan_server.py` submits the prepared FASTA batches to the authors' server at `http://andromeda.uab.cat/bioinf/aggrescan/`. It processes one job at a time, saves the received HTML and spreadsheet, and matches every returned amino-acid sequence to the exported sequence hash. Test sequences were explicitly authorized for this fixed annotation procedure; no solubility labels enter submissions or parsing.

```sh
python nesg_primary/run_aggrescan_server.py --validate-pilot
python -m unittest discover -s nesg_primary -p 'test_aggrescan_server.py' -v
python nesg_primary/run_aggrescan_server.py --run
python nesg_primary/validate_aggrescan_results.py
```

Network access is needed only for `--run`. Completed batches are verified and reused; an uncertain request stops the run rather than submitting a duplicate automatically. A run without `--run` rebuilds the combined tables from saved jobs. `--limit 1` submits at most one new batch; `--split train`, `validation` or `reserved_test` restricts submission to one split.

Results are under `results/aggrescan/`: `protein_annotations.tsv`, `annotation_status.tsv`, `residue_scores.tsv.gz`, `hotspots.tsv.gz` and `run_summary.json`. Original responses, job URLs, request times and checksums are under `aggregation_server/jobs/`. The server spreadsheet has an `.xls` extension but contains tab-separated text.

The legacy form requires CRLF line endings and converts the colon in `NESG_raw:100` to `3A`. The parser accounts for that observed name transformation, verifies the entire returned sequence, and checks that spreadsheet and HTML scores agree. Hotspot intervals come from the server's highlighted residues, avoiding threshold decisions based on rounded scores. Coordinates are 1-based inclusive.

These are predicted aggregation annotations. The three X-containing development sequences remain in the cohort with `unsupported_input`; they are never assigned zero aggregation. Proteins successfully annotated without hotspots use `zero_features`. An unfinished job remains `not_run`. Server defaults are used without fitting parameters to NESG labels. The server page does not expose a version identifier, so received bytes, dates, endpoint and settings are retained as provenance.

Completed server run: **3,759/3,762 proteins annotated** across 39 batches (training 2,359/2,361; validation 690/691; test 710/710). The three unsupported records contain X. Sixteen training proteins have valid predictions with zero hotspots. The saved outputs contain 925,761 residue rows and 33,092 hotspot intervals. All returned sequences, coordinates and summary/profile consistency checks passed; see `results/aggrescan/quality_checks.json`. Five parser tests also passed. Original cohort data and splits are unchanged.
