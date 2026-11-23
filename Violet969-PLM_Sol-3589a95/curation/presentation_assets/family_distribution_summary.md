# Named Pfam composition of the UESolDS operational cohort

Verified: **13,880 proteins = 9,880 training + 4,000 validation**; soluble (label 1) **7,996**, insoluble (label 0) **5,884**.
Both candidates use exactly this cohort and completed uniform Pfam 37.0 gathering-threshold annotations. No sequences, labels, splits, clusters or membership were changed; no independent-test files were accessed.

[Family PNG](pfam_family_distribution.png) · [Family SVG](pfam_family_distribution.svg) · [Clan PNG](pfam_clan_distribution.png) · [Clan SVG](pfam_clan_distribution.svg)
[Statistics](family_distribution_statistics.tsv) · [Protein assignments](family_distribution_assignments.tsv) · [Script](family_distribution_plot.py)

## Protein-level assignment and category selection

For each protein, combine all aligned residue intervals assigned to each Pfam accession (1-based inclusive); rank families by total union coverage, then greatest single aligned-domain span, then ascending accession. The top family is the representative. Repeated/overlapping hits to the same family do not double-count residues. Coverage favors large domains or repeated domains: this is a reproducible display convention, not a biological claim that one domain defines the entire protein.
Map that same representative family to its Pfam 37.0 clan for the broader candidate. A family without a clan remains explicitly “Representative family has no clan,” even if another domain has a clan. Clans are evolutionary groupings of domain families, not MMseqs clusters or complete protein architectures.
For each granularity, select seven named categories with the most representative assignments in the **combined** cohort, breaking count ties by accession. This selection and the representative rule use no solubility labels, sources or performance. Every protein contributes exactly once to each stacked bar; residual categories are exhaustive and disjoint. “No Pfam assignment” means no GA-passing hit, not absence of biological domains.
Display names are shortened biological descriptions, with accessions retained to distinguish similarly named families. Exact release descriptions are in the TSV.

## Family candidate: exact plotted counts

| Category | Soluble (7,996) | Insoluble (5,884) | Combined (13,880) |
|---|---:|---:|---:|
| Protein kinase (PF00069) | 27 (0.337669%) | 64 (1.087695%) | 91 (0.655620%) |
| Methyltransferase (PF13489) | 28 (0.350175%) | 50 (0.849762%) | 78 (0.561960%) |
| RNA recognition motif (PF00076) | 35 (0.437719%) | 24 (0.407886%) | 59 (0.425072%) |
| Radical SAM (PF04055) | 29 (0.362681%) | 25 (0.424881%) | 54 (0.389049%) |
| Extracellular solute-binding (PF13416) | 31 (0.387694%) | 20 (0.339905%) | 51 (0.367435%) |
| GNAT acetyltransferase (PF00583) | 25 (0.312656%) | 22 (0.373895%) | 47 (0.338617%) |
| LysR substrate-binding (PF03466) | 24 (0.300150%) | 20 (0.339905%) | 44 (0.317003%) |
| Other families | 6,115 (76.475738%) | 3,928 (66.757308%) | 10,043 (72.355908%) |
| No Pfam assignment | 1,682 (21.035518%) | 1,731 (29.418763%) | 3,413 (24.589337%) |

The seven named family categories cover **424/13,880 (3.05%)** of the cohort.

## Clan candidate: exact plotted counts

| Category | Soluble (7,996) | Insoluble (5,884) | Combined (13,880) |
|---|---:|---:|---:|
| FAD/NAD(P)-binding Rossmann fold (CL0063) | 239 (2.988994%) | 273 (4.639701%) | 512 (3.688761%) |
| Helix-turn-helix (CL0123) | 274 (3.426713%) | 163 (2.770224%) | 437 (3.148415%) |
| P-loop NTP hydrolase (CL0023) | 212 (2.651326%) | 211 (3.585996%) | 423 (3.047550%) |
| Phosphate-binding TIM barrel (CL0036) | 182 (2.276138%) | 149 (2.532291%) | 331 (2.384726%) |
| Protein kinase (CL0016) | 64 (0.800400%) | 119 (2.022434%) | 183 (1.318444%) |
| Tetratricopeptide repeat (CL0020) | 104 (1.300650%) | 68 (1.155676%) | 172 (1.239193%) |
| Periplasmic binding protein (CL0177) | 92 (1.150575%) | 72 (1.223657%) | 164 (1.181556%) |
| Other clans | 3,622 (45.297649%) | 2,225 (37.814412%) | 5,847 (42.125360%) |
| Representative family has no clan | 1,525 (19.072036%) | 873 (14.836846%) | 2,398 (17.276657%) |
| No Pfam assignment | 1,682 (21.035518%) | 1,731 (29.418763%) | 3,413 (24.589337%) |

The seven named clan categories cover **2,222/13,880 (16.01%)** of the cohort.

## Additional purple rows

| Quantity | Soluble | Insoluble |
|---|---:|---:|
| In a mixed-label family | 3,827/7,996 (47.861431%; displayed 47.9%) | 3,181/5,884 (54.061863%; displayed 54.1%) |
| Beyond 1,022 residues | 73/7,996 (0.912956%; displayed 0.9%) | 145/5,884 (2.464310%; displayed 2.5%) |

There are **6,043 distinct Pfam families**, including **1,660 mixed-label families** (at least one protein of each label). The mixed-family row uses **all original family assignments**, not only representatives. It is descriptive and uses both training and validation labels; it is not a held-out predictor or a minimum-support family contrast. Length is the original full-sequence length, and the threshold is strictly >1,022.

## All-hit prevalence before selecting representatives

These counts deduplicate protein–family or protein–clan pairs but allow proteins in multiple rows; they are **not** parts of a 100% bar. Top-seven lists below rank by combined prevalence independently of labels. Representative ranking differs because proteins can have multiple, sometimes overlapping alternatives; choosing a top-listed family whenever present would artificially privilege the displayed categories.

### Most prevalent family memberships

| Accession | Exact biological description | Proteins | Cohort percent |
|---|---|---:|---:|
| PF13649 | Methyltransferase domain | 151 | 1.087896% |
| PF08241 | Methyltransferase domain | 145 | 1.044669% |
| PF13489 | Methyltransferase domain | 113 | 0.814121% |
| PF00069 | Protein kinase domain | 108 | 0.778098% |
| PF00583 | Acetyltransferase (GNAT) family | 105 | 0.756484% |
| PF13847 | Methyltransferase domain | 101 | 0.727666% |
| PF07714 | Protein tyrosine and serine/threonine kinase | 100 | 0.720461% |

### Most prevalent clan memberships

| Accession | Exact biological description | Proteins | Cohort percent |
|---|---|---:|---:|
| CL0123 | Helix-turn-helix clan | 676 | 4.870317% |
| CL0063 | FAD/NAD(P)-binding Rossmann fold Superfamily | 571 | 4.113833% |
| CL0023 | P-loop containing nucleoside triphosphate hydrolase superfamily | 464 | 3.342939% |
| CL0036 | Common phosphate binding-site TIM barrel superfamily | 348 | 2.507205% |
| CL0020 | Tetratrico peptide repeat superfamily | 201 | 1.448127% |
| CL0159 | Ig-like fold superfamily (E-set) | 190 | 1.368876% |
| CL0016 | Protein kinase superfamily | 187 | 1.347262% |

## Multidomain and overlapping-assignment ambiguity

Multiple distinct Pfam assignments may describe distinct domains, repeated regions, or overlapping alternative profiles. They must not all be called separate biological domains. The existing same-clan competing-overlap flag is retained without resolving or rerunning annotations.

| Diagnostic | Soluble | Insoluble | Combined |
|---|---:|---:|---:|
| multiple_distinct_families | 2,223 (27.80%) | 1,694 (28.79%) | 3,917 (28.22%) |
| multiple_domain_hits | 2,440 (30.52%) | 1,851 (31.46%) | 4,291 (30.91%) |
| multiple_clans | 722 (9.03%) | 522 (8.87%) | 1,244 (8.96%) |
| coverage_tie | 34 (0.43%) | 28 (0.48%) | 62 (0.45%) |
| competing_clan_overlap | 1,276 (15.96%) | 1,105 (18.78%) | 2,381 (17.15%) |
| no_clan_despite_other_clan | 172 (2.15%) | 89 (1.51%) | 261 (1.88%) |

Alternative label-independent rules: highest domain bit score (ties: earliest start, accession), or earliest domain (ties: earliest end, accession). Bit scores are not uniformly calibrated across Pfam models, so that alternative is only a sensitivity diagnostic. Keep the primary top-seven categories fixed when comparing plotted distributions.

| Alternative | Level | Changed representative, combined | Changed plotted category, combined | Largest category shift, soluble / insoluble |
|---|---|---|---|
| highest_score | family | 1,303 (9.39%) | 161 (1.16%) | 0.538 / 0.850 percentage points |
| highest_score | clan | 466 (3.36%) | 306 (2.20%) | 0.513 / 0.510 percentage points |
| first_domain | family | 1,934 (13.93%) | 173 (1.25%) | 0.763 / 1.088 percentage points |
| first_domain | clan | 969 (6.98%) | 664 (4.78%) | 0.913 / 0.731 percentage points |

Exact alternative category counts, all-hit prevalence across all families/clans, and per-class ambiguity counts are saved in the statistics TSV. Per-protein representative and alternate assignments are saved separately.

## Interpretation and caption

The **clan candidate is the more interpretable of these two displays**, because it groups related domain families. However, neither top-seven display summarizes most of this biologically diverse cohort: residual categories dominate, particularly at family level. The figure should not imply that the named categories are the dominant biology of the entire dataset. A fuller family distribution would require a different visualization or many more categories.
Suggested caption: **“Soluble and insoluble proteins span diverse Pfam families; broader clan groupings reveal their domain-family composition.”**
Methods footnote: “Pfam 37.0; one coverage-based representative domain family per protein; seven categories selected by combined prevalence. Mixed-family membership uses all Pfam hits.”
This is descriptive composition, not evidence for a causal solubility mechanism or ESM representation. Broad clan names do not establish precise function for every member.

## Reference metadata and reproduction

No Pfam searches were rerun. Existing Pfam-A.clans.tsv supplies family descriptions and clan membership. Full clan names were extracted into the small pfam37_clan_names.tsv reference from pinned [Pfam 37.0 Pfam-C.gz](https://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam37.0/Pfam-C.gz), #=GF AC/ID/DE fields (Latin-1 decoding). Download MD5 verified: `5646960c822e6bada8be93afc5b4aa09`; SHA-256: `a3bd7c5fa1d07929b7a600277281cb6088a9738d6b1220a4f556ad444e4815ce`. No reference download is needed to regenerate figures.
From repository root:
```sh
python Violet969-PLM_Sol-3589a95/curation/presentation_assets/family_distribution_plot.py
```
pandas 3.1.1; matplotlib 3.11.2. PNGs have white backgrounds, 2640 × 2600 pixels, 400 dpi; SVG text remains editable. Bar proportions use exact integer ratios; only displayed labels round to one decimal. The TSV includes exact rational percentages.
Checks include cohort counts; stable ID/hash/split/label/cluster alignment; completed annotation checksum; GA thresholds; valid coordinates; domain/membership consistency; exclusive/exhaustive partitions; order-invariant assignment; and unchanged input hashes. No independent-test reads, annotation reruns or model execution. The original DeepLoc slide is unavailable locally; the panels follow the requested aligned bars and navy/purple/muted style.

Input SHA-256:

- `results/cohort_manifest.tsv`: `4368fe98cd9799f08461678859b1421153188ea7d47f744d3c6101be8b407059`
- `pfam_annotation/searched_manifest.tsv`: `a9c8d5d2d9f24b30a493da80045077f4bd464b303820aae7715e9b19c7f0a2ff`
- `pfam_annotation/annotated_cohort.tsv`: `13cf4c0fc57b0e102f7037cf86f3d8696d4f8aa0c9b8fa4b53e39c2ddef7c493`
- `pfam_annotation/domains.tsv`: `559c1a062e9eccb22a524395d150f7d4bcbcb4b429057a4a8babeeddaab5e653`
- `pfam_annotation/memberships.tsv`: `0e27586d6402ce294d60b79080dc6e87d1bf7611a7524afaf08f75b09b26b85e`
- `pfam_annotation/search_run.json`: `d6abd0b68a3326f831c3add1a0b4f3926a733843f95bd29a4d4b6f1cdeb3fdc5`
- `pfam_annotation/Pfam-A.clans.tsv`: `6f1e3f0cd744fb1ddbda4a71c83a504037d5aa99cf58136cd58a8ad5998f6574`
- `pfam_annotation/domains_flagged.tsv`: `32684bd91d18394e6dbfa059084cfee94d631ba1fb2387cac4f07090ca3c4015`
- `presentation_assets/pfam37_clan_names.tsv`: `0f7fb39958c6a002f18425e956adb9b5990e37a65d4790ce261eec913439f77a`
