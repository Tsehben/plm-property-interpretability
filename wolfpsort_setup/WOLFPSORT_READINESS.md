# Local WoLF PSORT readiness

Status: **installed, verified, usable.** Nothing was submitted to any web
server; the predictor runs in a container started with `--network none`. The
project's ESM/Python environment is untouched — WoLF PSORT lives entirely inside
a separate container image.

---

## 1. What was installed

| | |
|---|---|
| package | WoLF PSORT Package **Command Line Version 0.2** (packaged docs; no later version is distributed) |
| source | `https://github.com/fmaguire/WoLFPSort` — the maintainer's rehost after the original site went down |
| commit | `23786c0d17738c0ae1ce7e96d8131b551f50141d` (2023-04-06), tree `399d86a0b64ffd5a7118f530209f93f5d08ebaf7` |
| authors | Paul Horton; Kenta Nakai (original PSORT) |
| citation | Horton et al., *Nucleic Acids Research* 35 (suppl_2) W585–W587, 2007, doi:10.1093/nar/gkm259 |
| local copy | [src_download/wolfpsort_repo/](src_download/wolfpsort_repo/) (downloaded material, kept in its own directory; no script is run from inside it) |

### Licence — read the caveat

`LICENSE` was updated on 2023-04-04 and now states **`LICENSE: GPLv3`**. The
same file still carries its 2005 header `Copyright: All Rights Reserved` and the
original two-part academic/commercial clause, under which academics "may freely
use this package on site or redistribute it for academic use as is", with
**"this does not include providing a WoLF PSORT prediction service over the
internet"**, and "any rights not explicitly given in this license are reserved".

The file is therefore internally inconsistent: a GPLv3 grant and a
non-commercial, no-public-service restriction cannot both bind. For our use —
academic, non-profit, run locally, nothing redistributed, no service exposed —
every reading of the file permits it. **Before this tool's outputs appear in
anything redistributed or public-facing, the licence should be settled with the
author rather than assumed.**

### Dependencies

* **Perl ≥ 5.8** (fails on 5.005_03 per `INSTALL`) with `Getopt::Long`,
  `Pod::Usage`, `FileHandle`, `File::Basename`, `FindBin`. Debian's `perl-base`
  is not enough — `Pod::Usage` lives in `perl-modules`, so the full `perl`
  package is installed.
* **One compiled program**, `bin/wolfPredict` (plus `WoLFPSORTpredictAndAlign`
  and `runAddAlignmentToHtml`, used only by the HTML path).
* Optional, for the HTML path only: a `dbm` implementation.

Everything else — `runWolfPsortSummary`, `checkFastaInput.pl`,
`psortModifiedForWoLF` — is Perl.

### The one real obstacle, and how it was handled

**The distribution ships no source for its compiled program.** `bin/binByPlatform/`
contains only prebuilt binaries:

```
bin/binByPlatform/binary-i386/wolfPredict   ELF 32-bit LSB, Intel 80386, statically linked, for GNU/Linux 2.2.0
bin/binByPlatform/binary-sparc/wolfPredict  ELF 32-bit MSB, SPARC, statically linked
```

There is no C/C++ source, no `Makefile` and no `configure` anywhere in the
repository (verified across all 102 tracked files). This host is macOS 26.5.1 on
arm64, which can execute neither binary, and Rosetta 2 translates macOS x86-64,
not Linux i386 ELF. A native build is therefore **impossible, not merely
inconvenient** — there is nothing to build.

Resolution: run the published artefact unmodified under `linux/386` emulation in
an isolated container. This is a faithful installation of the distributed
program, not a substitute predictor.

| | |
|---|---|
| image | `wolfpsort:0.2`, id `sha256:2a0b57fa1dad484f46c411ac23273e2bf77760f1a1d973f4197fc12331a440e6` |
| base | `i386/debian@sha256:abd54731664e83f582c92db7581fb65a7bb81498c8b7001101844301bf6efe36` (`bookworm-slim`) |
| runtime | Docker 29.9.3, platform `linux/386` |
| binary used | `wolfPredict`, sha256 `f705cc11aa94b40c859c756424f67bf53cb91c7940aa6ccaa8508ad80e5c3171` |

`bullseye` was tried first and fails: Debian 11's i386 security pool no longer
serves the `perl` packages apt resolves to (404). `bookworm` works.

### Installation commands, exactly as run

```bash
git clone https://github.com/fmaguire/WoLFPSort.git wolfpsort_repo    # commit 23786c0d
docker build --platform linux/386 -f Dockerfile -t wolfpsort:0.2 src_download/wolfpsort_repo
```

[Dockerfile](Dockerfile) performs the package's own `INSTALL` steps: copy
`bin/binByPlatform/binary-i386/*` into `bin/` (step 1), then run
`psortModifiedForWoLF -t all.seq` to build the dbm feature tables (step 2; needed
only by the HTML path, and its failure would not affect `runWolfPsortSummary`).
Step 2 completed. `LANG=C` is set because the wrappers test `$ENV{LANG}` without
checking it is defined, which otherwise prints an uninitialised-value warning on
stderr.

### Running it

```bash
./wolfpsort.sh <animal|plant|fungi> <input.fasta> <output_dir>
```

[wolfpsort.sh](wolfpsort.sh) writes, all raw: `summary.txt` (stdout),
`summary.err` (stderr), `preflight.txt` (see §3), `run.json` (image id, input
and output hashes). [parse_wolfpsort.py](parse_wolfpsort.py) parses and
reconciles.

### Organism models

| model | reference sequences | feature rows | k | classes (single / dual) | extr | cyto |
|---|---|---|---|---|---|---|
| animal | 12,771 | 12,775 | **32** | 21 (10 / 11) | 3,130 | 1,438 |
| fungi | 2,158 | 2,158 | **27** | 14 (10 / 4) | 140 | 354 |
| plant | 2,426 | 2,426 | **14** | 24 (11 / 13) | 114 | 432 |

Class vocabularies differ by design (`chlo` only in plant, `lyso` only in
animal), and `k` differs. The documentation notes the program does not check
that the organism type matches the query's actual organism.

---

## 2. Verification — 37 checks, 37 passed

[verify_install.py](verify_install.py) · results in
[verification/verification_results.tsv](verification/verification_results.tsv)
· raw outputs preserved under [verification/raw/](verification/raw/) · fixtures
in [verification/fixtures/](verification/fixtures/).

Every query sequence comes from the package itself — `bin/testQuery.fasta` (the
author's own regression query) and `data/*.psort` (the packaged reference sets).
No project sequence was used.

| group | checks | result |
|---|---|---|
| 1. packaged regression | 6 | PASS |
| 2. single- and multi-sequence input | 3 | PASS |
| 3. identifier preservation | 3 | PASS |
| 4. animal / fungi / plant modes | 4 | PASS |
| 5. extracellular / cytoplasmic / dual parsing | 5 | PASS |
| 6. missing labels, invalid input, failed predictions | 13 | PASS |
| 7. determinism | 3 | PASS |

**Regression.** `runWolfPsortSummary` reproduces the `REGRESSION` file's recorded
output **byte for byte in all three organism modes**, with empty stderr. That
file was written by the author in 2006 against the Pentium build; the i386
binary under emulation on arm64 in 2026 gives identical numbers.

**Identifiers** are preserved verbatim up to the first whitespace, including
`sp|P12345|NAME_HUMAN`, dots, colons, hyphens and a 120-character id. Free-text
description is not absorbed. Duplicate identifiers produce two records and are
not collapsed; the parser matches them positionally.

**Modes.** All three run and report their own `k` (32 / 27 / 14) — recorded
precisely because it means scores are not comparable across models.

**Parsing.** Extracellular-, cytoplasmic- and dual-labelled reference queries all
parse. Dual classes (`cyto_nucl`, `extr_plas`, `cyto_mito`, `E.R._mito`) parse as
single tokens. *These fixtures are reference-set members, so the predictor
retrieves them as their own neighbours; nothing here is an accuracy estimate.*

**Determinism.** Three identical runs give a byte-identical `summary.txt`
(sha256 `ccad3e4098c8f151…` recorded in the results file), and a sequence scores identically
alone and inside a 10-sequence batch. The predictor is deterministic; in the
planned study all variability will come from the sequence edit, none from the
predictor.

### The finding that matters most: silent record loss

`checkFastaInput.pl` **drops records without telling the caller.** It writes
`#Warning:` lines to *stdout*, which the downstream feature computer then
consumes, so they never surface; and the exit status is `0` regardless.

In the edge-case fixture, **10 input records produced 6 output lines, exit 0,
empty stderr**. Dropped: three sequences under 30 residues and one containing
the residue `O`.

Documented rules, confirmed by probe:

| condition | behaviour |
|---|---|
| cleaned length < 30 | **dropped silently** |
| cleaned length > 10,000 | **dropped silently** |
| contains residue `O` | **dropped silently** |
| length exactly 30 | accepted (boundary verified) |
| digits, `*`, whitespace | deleted *before* the length check, then accepted |
| lowercase | upper-cased, accepted |
| does not start with `M` | warning only — accepted |
| headerless raw sequence | id becomes `queryProtein` |
| bare `>` header | ids become `queryProtein1`, `queryProtein2`, … |
| empty input | `Error: empty input`, plus a Perl array error from the feature computer — **still exit 0** |
| non-FASTA prose | letters are concatenated into a "sequence"; here it was dropped only because the text contained an `O` |

**Consequences adopted as standing rules.** Exit status is not a success signal.
Output must be reconciled against input identifiers; `wolfpsort.sh` therefore
runs `checkFastaInput.pl` separately to capture the warnings that the pipeline
swallows, and `parse_wolfpsort.py` emits one row per *input* record with an
explicit status. Every sequence entering the planned study must be pre-screened
for length ≥ 30 and absence of `O`.

---

## 3. What the numbers are — and are not

From the packaged documentation (`bin/runWolfPsortSummary`, `=head1 DESCRIPTION`):
each number "roughly corresponds to the number of nearest neighbours" carrying
that site, computed through the class-similarity matrix `data/*.wolfu`, so that
for a compound class

```
score(extr_plas) = #extr_plas + 0.5·#extr + 0.5·#plas
```

They are **similarity-weighted neighbour vote counts**. They are not
probabilities, they are not normalised, and across one line they need not sum to
`k` (fungi `RCC1_YEAST`: 12.5 + 11.5 + 9.5 + 3 = 36.5 against k = 27). The
documentation supports no probabilistic reading, so nothing in this project will
call them probabilities.

> **CORRECTED 2026-11-22 — a class absent from the printed line does NOT mean
> zero.** An earlier version of this document said it did. Verified against the
> implementation on 840 records: `runWolfPsortSummary` prints only a *selected
> subset* of the non-zero classes. Every one of those 840 records omitted at
> least one non-zero class (mean 7.6, max 15); `extr` was non-zero but omitted in
> 119 records and `cyto` in 89. Scores must therefore be recomputed for every
> class from the `--print-neighbors` list through `data/<model>.wolfu`, and that
> recomputation reproduces **every** printed score exactly (2,650 of 2,650, zero
> mismatches). See [../sequence_variant_pilot/PILOT_RESULTS.md](../sequence_variant_pilot/PILOT_RESULTS.md) §1.
> A class that recomputes to zero is a genuine zero; a record with no output line
> is missing. Those two remain distinct.

A direct consequence used below: `score("extr")` **already contains half-credit
for every `extr_*` neighbour**, so adding a compound score to the single-site
score double counts the same neighbours.

---

## 4. Reference-data audit

[audit_reference_data.py](audit_reference_data.py) · results in
[reference_audit/reference_audit.json](reference_audit/reference_audit.json) ·
per-sequence index in
[reference_audit/reference_sequence_sha256.tsv](reference_audit/reference_sequence_sha256.tsv).
Package-side only — the project cohort was not read and P0 was not touched.

**The package contains all three kinds of data, in different files:**

| file | contents | keyed by |
|---|---|---|
| `data/{animal,fungi,plant}.psort` | **full sequences** + class label, FASTA | UniProt **entry name** (`GCC1_HUMAN`) |
| `data/*.wolff` | **transformed feature matrix**, 56 features + class, **no sequence** | the same entry names |
| `data/*.wolfu`, `data/*.wolfw` | class-similarity matrix, feature weights, `k` | — |
| `bin/psortModifiedForWolfFiles/all.seq` | 1,552 full sequences + 3-letter PSORT site | UniProt **accession** (all 1,552 match accession syntax) |

Union across the three models: **17,345 distinct reference identifiers,
17,242 distinct sequences.** No sequence is shared between organism models;
83 are duplicated within a model.

`all.seq` feeds the classical-PSORT dbm tables used by the HTML path; it is *not*
consulted by `runWolfPsortSummary` at prediction time. The kNN prediction draws
on `*.wolff` / `*.psort`.

### Which overlap checks are technically possible

| check | possible | how |
|---|---|---|
| **exact sequence overlap**, cohort frozen sequence vs reference sequence | **yes** | SHA-256 of the upper-cased, whitespace-stripped residue string on both sides; covers all 17,242 distinct reference sequences |
| near-identity / homology overlap (e.g. ≥ 90% id) | **yes, in principle** | the sequences are present; needs an aligner (BLAST/DIAMOND/MMseqs2), none of which is installed here — a separate, declared step |
| **identifier overlap by UniProt accession** | **no** | the kNN files are keyed by **entry name**, not accession; mapping needs the UniProt ID-mapping service, i.e. a network query. Only `all.seq` is accession-keyed, and it does not drive the kNN prediction |
| overlap against the **feature matrix** as such | **no** | `*.wolff` carries only transformed features and no sequence. Overlap can be established only via the matching identifier in `*.psort`, never from a feature vector |
| label agreement with the frozen cohort labels | yes, but conditional | only for proteins first matched by exact sequence, and only under an explicit pre-registered mapping — the vocabularies differ (WoLF `extr`/`cyto` plus dual classes vs the binary DeepLoc labels) |

Because full sequences *are* supplied, sequence-overlap measurement is genuinely
available here; it would not have been if only `*.wolff` had shipped.

**No overlap was measured in this task.** Running it is a separate, declared
step, and it must not touch P0.

### Unresolved limitations

* The package records no provenance or date for how the reference sets were
  assembled, so their relationship to the DeepLoc 2.0 cohort cannot be
  established from the distribution alone.
* Entry names are not stable keys — they can be reassigned between UniProt
  releases — so even a successful name match would not prove record identity.
* Exact-hash overlap is a **lower bound on dependence**, not a measure of
  independence: a sequence absent from the reference set may still have close
  homologues in it. This is the same unmeasured homology dependence that limits
  the P4 evaluation.
* `data/animal.psort` holds 9 duplicated identifiers and `data/animal.wolff` has
  1 identifier with no sequence record, so the two animal files are not in exact
  bijection.
* One header in `data/animal.psort` (`>TCOF_HUMAN`) lacks the space after `>`;
  a naive `grep '^> '` undercounts the animal set by one. The audit parser
  handles it; the correct count is 12,771.

---

## 5. Proposed fixed evaluation readout

**Fixed now, before any variant exists**, and not to be revised in light of
project variant outcomes.

### Definition

For a sequence *s* that produced a prediction, from the single-site classes only:

```
E(s) = score of class "extr"
C(s) = score of class "cyto"
Δ(s) = E(s) − C(s)
```

**Extraction (corrected — see §3):** `E` and `C` are read from the full class
score vector recomputed from the `--print-neighbors` list through
`data/animal.wolfu`, **not** from the printed summary line, because the printed
line omits non-zero classes. The recomputation is asserted to reproduce every
printed score exactly, on every record.

The paired, per-protein statistic over an original *o* and its variant *v*:

```
D(o, v) = Δ(v) − Δ(o)
```

in similarity-weighted neighbour votes. Positive means the edit moved the
sequence toward the extracellular reference neighbourhood relative to the
cytoplasmic one. Summarised across proteins with the protein as the sampling
unit, exactly as in the P4 evaluation.

### Treatment of dual labels

Compound classes are **not** added into `E` or `C`. `score("extr")` already
carries half-credit for every `extr_*` neighbour, so adding `extr_plas` would
count the same neighbours twice. Each record instead carries flags —
`top_is_compound`, `top_mentions_extr`, `top_mentions_cyto`,
`top_ambiguous_for_readout` (top class names both sites) — which are reported as
coverage, and one **pre-registered** sensitivity analysis excludes proteins whose
original's top class is compound. The flags never alter the primary readout.

### Treatment of missing outputs

A protein contributes only if **both** members of its pair produced a prediction.
If either is absent, the pair is dropped, counted, and the recorded reason is
reported — **never imputed as zero**. Dropping is distinguished from a genuine
zero: an `extr` score that *recomputes* to zero means no neighbour carried that
site (a class merely missing from the printed line does not — see §3); an absent
record means no prediction exists. Pre-screening (length ≥ 30, no `O`,
length ≤ 10,000) is applied before generation so this should be rare, and any
occurrence is reported rather than silently absorbed.

### Limitations

1. **Not probabilities.** Bounded by `k` and by the composition of the reference
   set; `D` is a difference of vote counts, not of likelihoods.
2. **Not comparable across organism models.** `k` differs (32 / 27 / 14), the
   reference sets are disjoint, and the class vocabularies differ. A `D` of +3
   under `animal` and +3 under `plant` are different quantities. **The study will
   therefore fix one organism model — `animal` — for every protein, and restrict
   the sample to Metazoa so that the chosen model matches the query organism.**
   No cross-model pooling, no per-protein model switching.
3. **Not a localisation call for a single sequence.** The verification run shows
   two reference proteins labelled `extr_plas` whose Δ is *negative* because
   `cyto` outscores the single-site `extr`. `D` is meaningful as a paired
   difference, not as an absolute classification.
4. **Predictor, not experiment.** WoLF PSORT is itself a kNN predictor over
   UniProt-derived annotations. Agreement between an ESM-2 intervention and WoLF
   PSORT is agreement between two sequence-based models. It is not experimental
   validation, and both may share annotation and sampling biases.
5. **Possible reference contamination, currently unmeasured.** See §4.
6. **Coarse granularity.** Scores are sums of (half-)integer votes, so `D`
   changes in steps of 0.5; small edits may well produce `D = 0` exactly, and a
   zero is a real tie, not a missing value.

---

## 6. Files

| path | purpose |
|---|---|
| [Dockerfile](Dockerfile) | the isolated `linux/386` image |
| [wolfpsort.sh](wolfpsort.sh) | run wrapper; preserves raw stdout, stderr, preflight, provenance |
| [parse_wolfpsort.py](parse_wolfpsort.py) | parser, input/output reconciliation, the §5 readout |
| [verify_install.py](verify_install.py) | the 37-check suite |
| [audit_reference_data.py](audit_reference_data.py) | the §4 audit |
| [verification/](verification/) | results, fixtures, preserved raw outputs |
| [reference_audit/](reference_audit/) | audit JSON and the reference sequence SHA-256 index |
| [src_download/wolfpsort_repo/](src_download/wolfpsort_repo/) | the downloaded package, unmodified |
