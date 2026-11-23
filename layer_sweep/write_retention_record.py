#!/usr/bin/env python3
"""Record exactly which pooled embeddings were retained and which were not.

Usage: write_retention_record.py <v_dir>
"""
import os, sys, json
V = sys.argv[1]
prov = json.load(open(os.path.join(V, 'pooled_embedding_provenance.json')))
sel = json.load(open(os.path.join(V, 'layer_selection.json')))
EVAL = sorted({int(sel['selected_layer']), 8})
ALL = list(range(1, 34))
notret = [l for l in ALL if l not in EVAL]
f = os.path.join(V, 'pooled_evaluated_layers.npz')
rec = {
 "layers_extracted": ALL,
 "layers_whose_pooled_vectors_ARE_retained": EVAL,
 "retained_file": "layer_sweep/v1/pooled_evaluated_layers.npz",
 "retained_file_bytes": os.path.getsize(f) if os.path.exists(f) else None,
 "layers_whose_pooled_vectors_are_NOT_retained": notret,
 "n_layers_not_retained": len(notret),
 "WHAT_THIS_MEANS": (
   "The pooled embeddings for the 31 non-evaluated layers were held in RAM for the duration "
   "of the sweep and were NOT written to disk. They cannot be reloaded. Any future analysis "
   "that needs them - refitting those probes, a different pooling, a different classifier, a "
   "cluster-held-out re-split at a non-evaluated layer - REQUIRES RE-EXTRACTION."),
 "why_not_retained": {
   "reason": "measured disk constraint",
   "needed_gb_float16_all_33_layers": 1.71,
   "free_gb_at_decision": 1.34,
   "swap_used_gb": 16.4, "swap_total_gb": 17.4,
   "note": "swap could not shrink because the volume was ~99% full, so the free space was not "
           "recoverable. The cohort (10,104) and the layer range (1-33) were NOT reduced."},
 "what_IS_retained_for_every_layer": [
   "per-layer SHA-256 of both pooled arrays (pooled_embedding_provenance.json)",
   "P3 predictions for all 33 layers x 3 feature sets (p3_predictions.npz)",
   "per-layer metrics (p3_layer_curves.tsv)",
   "fitted coefficients and intercepts (probe_parameters_all_layers.npz)",
   "training-fitted StandardScaler mean_ and scale_ (probe_parameters_all_layers.npz)",
   "convergence diagnostics (probe_convergence_diagnostics.tsv)",
   "extraction provenance (pooled_embedding_provenance.json)"],
 "re_extraction_command": (
   "python3 layer_sweep/run_full_sweep.py . <out_dir> --batch-tokens 6144    "
   "# same checkpoint, layer indexing, truncation, pooling and compute dtype; "
   "per-layer SHA-256 in pooled_embedding_provenance.json verifies an identical re-extraction"),
 "extraction_settings_needed_to_reproduce": {
   k: prov[k] for k in ('checkpoint', 'implementation', 'transformers_version',
                        'torch_version', 'device', 'compute_dtype', 'pooling_dtype',
                        'stored_dtype', 'pooling', 'max_residues', 'nterm_positions',
                        'batch_tokens', 'layers_extracted', 'accession_order_sha256')},
}
json.dump(rec, open(os.path.join(V, 'embedding_retention.json'), 'w'), indent=1)
print(json.dumps({k: rec[k] for k in ('layers_whose_pooled_vectors_ARE_retained',
                                      'n_layers_not_retained', 'WHAT_THIS_MEANS')}, indent=1))
