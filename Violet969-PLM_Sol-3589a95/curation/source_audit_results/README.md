# Source audit reproduction

All new results are under this directory. Earlier analyses and the frozen 13,880-protein cohort are read-only inputs. The independent test is never opened. Run from `curation/` with the environment pinned in `../requirements-source-audit.txt` (the project currently uses `/tmp/uesolds-venv/bin/python`).

```sh
python source_validity_audit.py --stage all
python recover_targettrack_sample.py --xml /tmp/uesolds-source-audit/tt.xml.gz
python audit_dnasu_sample.py --pages source_audit_results/dnasu_raw_pages.zip
python source_audit_supplement.py
python source_audit_uncertainty.py
python source_audit_evidence.py
python render_source_audit.py
python -m unittest test_source_validity_audit.py -v
```

`source_validity_audit.py` reuses completed annotations and cached ESM vectors. It performs lightweight logistic fits only. Stage `diversity` skips probes; stage `probes` skips diversity. The full run took about 8.5 minutes on this machine. There is no hyperparameter tuning, SAE fitting, ESM extraction, fine-tuning or intervention.

The original public archive URL/checksums are in `targettrack_archive.json`. Download it outside the checkout, verify its published MD5, and extract only `TargetTrack-1Jul2017/TargetTrack XML files/tt.xml.gz` to the path above. The parsed gzip SHA-256 is recorded in `targettrack_recovery_run.json`. Compact sampled XML and referenced protocols are kept in `targettrack_sample/`; DNASU raw pages are archived for offline reproduction. Research evidence was retrieved on 2026-11-12. The source evidence table is manually curated from the linked primary publications and recovered records; rerendering does not query live websites.

The fixed recovery samples deliberately cover split, label and length; their success rates must not be generalized as population rates. No experimental outcome is inferred from an unmatched or missing record. `source_audit_supplement.py` annotates evidence flags on the new canonical policy manifest, without editing original labels or files. Rerunning only the diversity stage resets those evidence fields; run the supplement afterwards to restore them.

Key result tables:

- `source_provenance.tsv`: programme facts versus record evidence and unresolved questions.
- `source_diversity.tsv`, `source_family_clan_distribution.tsv`, `family_source_support.tsv`: source contribution and support; family memberships can overlap.
- `leave_source_out_diversity.tsv`, `source_pair_overlap.tsv`, `source_biochemical_neighbour_loss.tsv`: marginal losses and redundancy.
- `source_label_descriptors.tsv`, `taxonomy_coverage.tsv`: source/split/label distributions and missingness.
- `source_prediction_metrics.tsv`, `source_prediction_per_centre.tsv`, `source_paired_differences.tsv`, `source_confusion.tsv`: multiclass source probes and uncertainty.
- `source_within_family*.tsv`: conditional source prediction; overlapping family strata are not independent observations.
- `solubility_metrics.tsv`, `solubility_paired_differences.tsv`: fixed-population source-removal sensitivities. Supplemental cross-source point estimates in the former have empty CI fields where not prespecified for bootstrapping; these are not zero-width intervals.
- `cross_centre_uncertainty.tsv`: all single-centre training transfer results with group-bootstrap intervals.
- `solubility_family_strata.tsv`: per-family and family/source AUROC, AP, Brier, log loss and intervals.
- `policy_comparison.tsv`, `policy_source_counts.tsv`, `policy_manifest.tsv`: candidate policies and evidence flags, without changing dataset membership.

All probability tables are saved at 17-digit precision. Reused historical solubility predictions retain the precision of their existing artifact. Bootstrap intervals are conditional on the fitted model, use 200 group resamples, and do not incorporate label or training uncertainty. SSGCID has only one negative validation record: its bootstrap intervals are not reliable population-level evidence.

Core-plus-extension eligibility is a scientific decision flag, not automatic thresholding of a probe. There is presently no certified comparable binary assay core. Two sampled MCSG soluble stages and twelve matching DNASU status fields are exposed separately as recovered evidence, not silently promoted to verified labels.
