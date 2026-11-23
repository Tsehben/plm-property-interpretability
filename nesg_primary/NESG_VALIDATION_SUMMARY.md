# NESG dataset validation summary

## Ready with limitations

**Ready for heterogeneous ESM analysis on training/validation.** Broad confirmatory within-family analysis needs further annotation and architecture validation. The [executed notebook](NESG_Dataset_Validation.ipynb) documents the evidence in six sections and ten figures.

## Final phenotype and cohort

The cohort contains **3,762 proteins: 2,614 soluble and 1,148 insoluble (69.48% soluble)**. All have E ∈ {4,5}; soluble means S ∈ {4,5}, insoluble S=0. The phenotype is recombinant soluble recovery conditional on high expression. Family proportions and missing-annotation records are preserved (§1, Figures 1–2).

The 9,703-record archive yields 9,272 singleton records; all 431 repeated-sequence observations remain in the exclusion ledger. Historical 9,644-record membership is unresolved. NetSolP's 1,323-row release derives from the same experiments and has three fractional-E threshold exceptions (§1).

## Experimental validity and technical uncertainty

Gel-based E/S scores describe expression and recovery in E. coli. Record-level vector, tag, batch, protocol and replicate metadata are unavailable. Soluble proportions are 685/973 (70.40%) for E=4 and 1,929/2,789 (69.16%) for E=5; this does not establish assay equivalence. No terminal His-run occurs in the distributed strict sequences, but tags may have been omitted from those sequences (§2, Figure 3).

Fractional E and E=0/S>0 observations remain flagged. The 437 candidate construct-ID groups are unverified. Three strict sequences contain X; they remain included with missing physicochemical summaries.

## Biological diversity and controlled contrasts

The cohort retains **1,383 Pfam families, 243 clans, 1,820 existing homology components and 390 mixed-label families**. Pfam coverage is 3,510/3,762 (93.30%); 1,337 proteins carry multiple family assignments. All sequences fit within 1,022 residues. Natural family and biochemical associations are retained (§3, Figures 4–6).

Training-selected controlled-family support (§5, Figure 8):

| Pfam policy | Groups/class in train | Train families | Validation ≥1/class | Validation ≥threshold/class | Eligible train proteins | Eligible validation proteins |
| --- | --- | --- | --- | --- | --- | --- |
| all_GA | 3 | 31 | 16 | 7 | 359 | 111 |
| all_GA | 5 | 13 | 8 | 2 | 228 | 80 |
| no_competing_overlap | 3 | 7 | 1 | 0 | 52 | 3 |
| no_competing_overlap | 5 | 1 | 0 | 0 | 0 | 0 |

A mixed-label homology group can count in both classes; the notebook also reports class-exclusive-group sensitivities. Protein counts deduplicate across families. **The single conservative ≥3 family is Radical SAM (PF04055), with three validation proteins; it fails the class-exclusive-group sensitivity.** Seven/two all-GA families meet ≥3/≥5 groups per class in validation as well. Check competing calls and architectures before confirmatory analysis.

The separate training-only family-held-out design has 695 joint family/homology components; the largest contains 153 proteins (6.48% of training). All five folds contain both labels among annotated proteins, with unequal sizes (§5, Figures 9–10).

## Final split and integrity

| Partition | Proteins | Soluble | Insoluble | Homology groups |
| --- | --- | --- | --- | --- |
| train | 2361 | 1645 | 716 | 1090 |
| validation | 691 | 491 | 200 | 364 |
| reserved_test | 710 | 478 | 232 | 366 |

Existing audit fold 0 is reserved test (18.87%), fold 1 validation, and folds 2–4 training. Five additional deterministic group-hash CV folds use **training only**. No split/seed search or performance-based selection was used (§4, Figure 7).

**Integrity checks passed:** sequence hashes and E/S agreement, unchanged membership, and no exact-sequence, homology-group or detected-MMseqs2-link overlap across splits. Training CV also keeps groups together. Manifests are hashed. The locked notebook run executed 13 code cells and ten figures without opening restricted records.

**Subsequent test exposure:** the user later authorized inspection of test labels by Pfam family. [TEST_LABEL_EXPOSURE.md](TEST_LABEL_EXPOSURE.md) records this; the test is no longer uninspected for family–label composition. Ordinary notebook reruns still use frozen aggregates and development records. The UESolDS independent test remains untouched.

## Recommended next step

Extract frozen ESM-2 activations for training/validation with explicit X handling. Tune in training group CV and validate family architectures using development data. Specify feature selection, readouts and controls before SAE interventions. Do not use the test-family report to revise model selection or hypotheses.

Reproduction instructions, dependency pins and manifest locations are in [README.md](README.md). `results/test_lock.json` records the full-manifest hash and access policy.
