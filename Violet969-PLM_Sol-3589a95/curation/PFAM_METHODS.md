# Uniform Pfam annotation: frozen protocol

All 13,880 proteins in the cached subset of `clean_all_sources` are searched:
9,880 training and all 4,000 published validation records. The script selects only
`clean & cached` from the previously frozen manifest; no source, label, annotation
or probe score enters selection. Original full sequences are recovered from the
existing pairs file and checked against their SHA-256 and length. No trimming,
truncation, relabeling, resampling, embedding extraction or test access occurs.

## Reference and search thresholds

- **Pfam 37.0**, fixed release directory:
  <https://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam37.0/>.
- 21,979 Pfam-A HMMs; profile archive MD5
  `8272c3f18626bc300ed7cc00e86bad40`. Validate all three downloaded archives
  against the release's `md5_checksums`, and record independent SHA-256 hashes.
- **PyHMMER 0.11.2** on the local CPU; search every profile against the same
  13,880 full sequences (`hmmsearch`, six query-parallel threads).
- Use **both sequence and domain gathering cutoffs** supplied by each Pfam
  profile (`bit_cutoffs='gathering'`, equivalent threshold policy to `--cut_ga`).
  A domain is retained only when both its target hit and domain are included.
  No cohort-fitted global E-value, score, length or HMM-coverage cutoff is added.
- Prespecify seed 42, sequence database size `Z=13880`, default significance-count
  domain search-space estimation, MSV/Viterbi/Forward filters
  `F1=0.02`, `F2=0.001`, `F3=0.00001`, composition-bias filtering and null2 correction.
  Record scores, conditional/independent E-values, model lengths, per-model GA
  thresholds and 1-based inclusive domain/alignment/envelope coordinates.

The [HMMER manual](https://eddylab.org/software/hmmer/Userguide.pdf) defines
model-specific GA sequence and domain thresholds. Pfam release background is
available in the [Pfam documentation](https://pfam-docs.readthedocs.io/en/dev/pfam.html).
This is a uniform HMM-based accession-membership analysis, not the complete
PfamScan post-processing pipeline. Family HMM coverage is the aligned model span
divided by model length; it is distinct from sequence coverage and from the
MMseqs2 bilateral coverage threshold.

## Family and controlled-analysis definitions

Primary annotation retains all GA-passing domains. Multiple copies of one domain
within a protein count once toward protein–Pfam membership; a multidomain protein
can belong to multiple families. Union counts deduplicate original `raw_uid`.
Single-class families and searched/no-hit proteins remain in primary A.
Pfam family accessions are not replaced with clans, MMseqs2 groups or organisms.
Holding one Pfam domain constant does not hold the complete protein architecture
constant; this limitation applies to all B claims.

Count ≥3 and ≥5 proteins **per label**, with separate support tables for training,
validation and pooled development. For B, select families using training support;
then report whether validation has both labels (minimal evaluability) and whether
it independently reaches ≥3/≥5 per label (stronger support). Apply the same rule
within a Pfam × source/centre stratum for same-source contrasts. Class-exclusive
MMseqs2 groups per class are additionally counted; meeting protein thresholds alone
is not treated as independent-group support. Original clusters and splits are
preserved; neither families nor clusters trigger a new split.

The same-source key remains a header proxy. It cannot establish matched organism,
expression host, construct, assay protocol or verified negative outcome. Pfam
annotation does not improve the reliability of a workflow-derived label.

Quality sensitivities are fixed before reading family support:

1. Require HMM coverage ≥80%, assessing partial-domain evidence; retain all A
   proteins, marking those without a surviving assignment as no usable hit.
2. Remove both members of competing **different-family, same-clan** domains with
   any aligned-residue overlap. No label-dependent tie-break is used. This is a
   deliberately conservative ambiguity sensitivity, not official PfamScan clan
   resolution; legitimate nested assignments can be lost. Preserve primary hits.
3. Require domain alignment end ≤1,022 for annotation accessible to cached ESM
   inputs. Whole-sequence primary annotation remains untruncated.
4. Previously defined source-controlled, terminal-artifact-free and ≤1,022-length
   subsets; these are diagnostic cohort restrictions, never new primary selection.

No model prediction/performance file is read. No label is changed and no
experimentally verified label is inferred from a Pfam hit. Existing pilot and
family-feasibility results remain historical snapshots.

## Reproduction

Keep the large reference library outside the checkout, e.g. `/tmp/uesolds-pfam37`.
Download `Pfam-A.hmm.gz`, `Pfam-A.clans.tsv.gz`, `Pfam.version.gz` and
`md5_checksums` from the fixed directory above. The annotation script verifies
the MD5s and release version before searching. Install the existing analysis
requirements plus the pinned package:

```sh
/tmp/uesolds-venv/bin/pip install -r Violet969-PLM_Sol-3589a95/curation/requirements-pfam.txt
/tmp/uesolds-venv/bin/python Violet969-PLM_Sol-3589a95/curation/annotate_pfam.py --reference /tmp/uesolds-pfam37 --cpus 6
/tmp/uesolds-venv/bin/python Violet969-PLM_Sol-3589a95/curation/summarize_uniform_pfam.py
```

Both scripts accept `--out /path/to/results` for isolated reproduction. Search
writes `domains.partial.tsv` and only promotes it to `domains.tsv` after all
21,979 models complete. Summarization requires a completed run and matching domain
hash. Interrupted searches are rerun; partial outputs are not interpreted as
no-hit annotations. The domain table, profile-specific thresholds, complete
searched manifest, source files for small reference metadata and run provenance
are retained. The large profile library is downloadable, not committed. No FASTA
copy is left in the working tree.

After summarization, render the report and validate the completed run:

```sh
/tmp/uesolds-venv/bin/python Violet969-PLM_Sol-3589a95/curation/render_uniform_pfam_report.py
/tmp/uesolds-venv/bin/python Violet969-PLM_Sol-3589a95/curation/validate_uniform_pfam.py --reference /tmp/uesolds-pfam37
```

Validation repeats the first 50 profiles against all 13,880 sequences with a
single CPU thread, verifies exact domain coordinates and matching scores against
the six-thread run, independently recomputes primary family thresholds using
cross-tabulations, and checks complete accounting and unchanged original inputs.
These checks do not tune annotation thresholds. The report renderer reads the
default output directory. Historical run-provenance manifests refer to their
respective snapshots, not to later documentation updates.
