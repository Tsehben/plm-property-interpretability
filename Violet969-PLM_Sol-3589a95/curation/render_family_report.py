#!/usr/bin/env python3
"""Render the family feasibility decision from saved measurements."""
from pathlib import Path
import pandas as pd
P=Path(__file__).resolve().parent
R=P/'family_results'
s=pd.read_csv(R/'feasibility_summary.tsv',sep='\t')
c=pd.read_csv(R/'coverage_by_source_label.tsv',sep='\t')
def table(d):
    return '| '+' | '.join(d.columns)+' |\n| '+' | '.join(['---']*len(d.columns))+' |\n'+'\n'.join('| '+' | '.join('NA' if pd.isna(x) else str(x) for x in row)+' |' for row in d.itertuples(index=False,name=None))+'\n'
text='''# UESolDS: heterogeneous and within-family feasibility

> Historical pilot-only snapshot. The completed [uniform Pfam annotation](UNIFORM_PFAM_REPORT.md) supersedes the annotation-coverage and zero-support conclusions below; original pilot measurements are preserved.

## Recommendation

Use **the existing cached subset of `clean_all_sources` as the primary experimental cohort**: 9,880 training + 4,000 published validation proteins (13,880 total), drawn from the same cleaned full-development parent of 69,248 training + 4,000 validation proteins. This preserves all represented source groups, labels, observed families and organisms, including single-class families and unknown annotations. It reuses the existing embeddings without additional model computation. No family or label balancing is performed. The inherited 10,000-training sample was proportional to source × label; the published validation set was already balanced by its authors, so its prevalence is not a population estimate.

This updates the previous recommendation of source-controlled UESolDS as primary. Under the clarified design, source control becomes an essential **sensitivity analysis** and matching criterion for B; it should not remove broader biological variation from A. The choice uses cohort support and scientific scope, not ESM AUROC. The family scripts never read predictions or performance tables.

**Current verdict:** the primary cohort supports heterogeneous sequence/property exploration (A), but family/organism representation-overlap claims remain severely annotation-limited. It does **not yet support an evidenced within-family contrast (B)** at either requested threshold. Therefore it cannot currently support both analyses with meaningful family controls. Zero observed qualifying families is an annotation-limited lower bound, not evidence that the underlying dataset lacks mixed-label families.

## One cohort, two linked analyses

**A — heterogeneous, primary.** Train the property probe on all eligible training records in the chosen scope, retaining unannotated proteins and families observed with only one label. Report source-stratified outcomes. Family/organism concept probes can use observed annotations only; report coverage and selection bias and never encode an unqueried family as a negative biological label. Preserve the same train/validation split. A homogeneous annotation-only subset must not silently replace A.

**B — within-family, controlled.** Derive membership strata from A using an observed Pfam accession and both labels. Report ≥3 and ≥5 proteins per class separately. For comparable-source contrasts, require those counts within the same Pfam × centre/source-proxy stratum. Source equality is only a proxy for assay comparability; record-level experimental outcome and conditions still need verification. B is a controlled evaluation/sensitivity of the primary model, not the only training set. Any future separately fitted within-family model must be clearly secondary. The empty B membership manifest is intentional.

Pfam accessions and MMseqs2 clusters are separate columns and separate concepts. A protein with several observed Pfam assignments belongs to several membership strata; union protein counts deduplicate `raw_uid`. A shared domain does not hold whole-protein identity or domain architecture constant. The unordered `pfam_set` is retained for inspection, not called an architecture. No Pfam assignments are propagated through MMseqs2 clusters, no Pfam clans are substituted for families, and no unknown membership is inferred from an identifier.

## Annotation coverage and missingness

Only the existing 100-protein annotation pilot was queried. Joining by original ID plus exact original-sequence hash and checking labels retains all 100 in the cleaned parent. Forty-nine proteins have at least one Pfam cross-reference, giving 67 protein–Pfam memberships and **66 distinct accessions**. All 66 are single-class in the observed sample. Only PF01381 appears on more than one protein: two soluble MCSG records in two detected homology groups. All remaining accessions occur once. All 49 annotated proteins are training records; validation has no attempted annotations.

| Scope | Total proteins | Attempts | Pfam-annotated proteins | Coverage | Distinct Pfam | Mixed-label Pfam |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Full cleaned parent | 73,248 | 100 | 49 | 0.0669% | 66 | 0 |
| Primary cached cohort | 13,880 | 13 | 7 | 0.0504% | 11 | 0 |

Among full-parent training records, Pfam coverage is 20/27,266 negatives (0.0734%) and 29/41,982 positives (0.0691%). Conditional on a pilot attempt it is 20/39 (51.3%) versus 29/61 (47.5%). These small pilot denominators do not establish a label-dependent missingness mechanism. The other 73,148 parent records are unqueried, not failed annotations. Of the 51 attempted records without Pfam, 28 lack an exact match and 23 have a match but no assignment through the historical retrieval procedure. Such absence is not proof that the protein lacks a Pfam domain.

Full-parent training coverage by source, shown as **Pfam assigned / attempted / all records**. A dash means there are no proteins in that class/source combination. All validation source/label strata have zero attempted annotations; their denominators are preserved in the downloadable table.

'''
rows=[]
cc=c[(c.policy=='primary_clean_all_sources')&(c.scope=='full')&(c.split=='train')]
for center,g in cc.groupby('center'):
    row={'source / centre':center}
    for label in [0,1]:
        a=g[g.label==label]
        row[f'label {label}: assigned / attempted / total']='—' if a.empty else f'{int(a.usable_pfam.iloc[0])} / {int(a.attempted.iloc[0])} / {int(a.n.iloc[0])}'
    rows.append(row)
text+=table(pd.DataFrame(rows))
text+='''
![Pfam coverage by source and label](family_results/pfam_coverage.png)

`coverage_by_source_label.tsv` distinguishes not queried, attempted without assignment, original assignments, and assignments usable under the evidence policy. Pfam annotations were collected via UniProt/RCSB cross-references in the earlier pilot; a Pfam release was not pinned and domain coordinates/scores were not recorded. They are observed accession memberships, not a new uniform profile-HMM scan. The pilot includes one tag-trimmed exact match but it has no Pfam assignment; all 49 assigned records use full-sequence exact matches. The full-sequence-only evidence sensitivity therefore changes no counts. Species annotation is also limited: only ten cached-cohort proteins have confident organism assignments. No taxonomic diversity is inferred from missing species metadata.

## Family support and independent-group counts

The following results hold separately for **both ≥3 and ≥5 proteins per class**, in the full parent and primary cached cohort:

- Qualifying Pfam families: **0**; eligible unique proteins: **0**; eligible detected homology groups: **0**.
- Qualifying same-Pfam/same-source strata: **0**; therefore their source composition and label reliability are **not estimable**, rather than an empty evidence base being called reliable.
- Families meeting each threshold in both training and validation: **0**.
- Families with at least the threshold number of class-exclusive homology groups per class: **0**.

`family_support.tsv` reports every observed family, including single-class families, separately for training, validation and combined development, and for full versus cached scope. `family_source_support.tsv` adds the centre. Counts of proteins, groups per class, groups shared across labels, class-exclusive groups, and distinct sources are separate. `source_label_counts` retains every source × label cell.

The parent's 67,207 detected groups and the cached cohort's 13,716 groups describe diversity across **all** records, not within qualifying families. MMseqs2 groups use the earlier 30%-identity/80%-bilateral-coverage analysis; distinct groups are operational independence units, not a guarantee of unrelatedness. A family with five proteins in one group does not supply five independent contrasts. Combined-development support is descriptive: it cannot substitute for threshold support in both published splits. Future group-aware evaluation should also avoid repeated proteins appearing in multiple family strata being counted as independent observations.

## Quality-filter sensitivity

All rows below use full-parent development scope. `qualifying_families`, eligible proteins/groups and same-source qualifying families are zero at **both** thresholds for every policy; the table shows the remaining support. These filters diagnose consequences and do not replace the primary cohort.

'''
v=s[(s.scope=='full')&(s.split_scope=='development')&(s.threshold_per_class==3)][['policy','n_cohort','n0','n1','annotated_proteins','distinct_pfam','mixed_label_families','qualifying_families']]
text+=table(v)
text+='''
- **Source-controlled:** the eight training-supported SG centres, useful for protocol/source sensitivity; it loses 17 annotated proteins and narrows source coverage.
- **No terminal artifact flags:** excludes any prespecified terminal His/cleavage/GS-rich flag. It disproportionately removes positive records and removes no Pfam-annotated records here. Since these are heuristic flags, this is not a justified primary exclusion or tag-stripping rule.
- **Within ESM context:** retains length ≤1,022 to assess input truncation. It loses three annotated proteins and disproportionately affects negatives; native long proteins stay in primary A.
- **Confident organism only:** collapses the parent to 69 pilot records (ten cached). This mostly selects annotation availability and is unsuitable as primary A.
- **Strict verified labels:** remains empty. No source prefix or Pfam match establishes a construct-matched solubility outcome.
- **Full-sequence Pfam evidence:** keeps every A protein but restricts usable family memberships to exact full-sequence matches. It leaves all 49 annotated proteins intact.

The full/cached, split-specific ≥3/≥5 results are in `feasibility_summary.tsv`. A label-reliability table identifies structural-genomics workflow-derived labels, structure-selected positives, DNASU outcome candidates without retrieved statuses, and unresolved-source records. All remain **record-level unverified**. No hypothetical mixed-label family is declared experimentally validated.

## What is needed before B is feasible

Annotate **all 13,880 primary cached-cohort sequences uniformly**, retaining unsuccessful searches as explicit attempted/no-hit outcomes, before choosing families. Use a pinned Pfam release and sequence-based domain search with recorded thresholds, coordinates and evidence; exact-reference mappings can supplement it. Annotate training and development without querying the locked independent test. Do not preferentially annotate only negatives or families already known to contain both labels. Retain the full cleaned parent as the expansion pool, rather than augmenting only desired family/label combinations.

Then rerun the same ≥3/≥5 and independent-group/source-support counts. Those thresholds are feasibility screens, not adequate-power guarantees. Select B family eligibility using training support; validation counts establish evaluability, not inclusion in primary A. Freeze the cohort and eligibility protocol before subsequent representation comparisons. If no adequate B exists after uniform annotation, report that limitation rather than balancing families or training exclusively on the surviving families. Original assay evidence is a separate recovery task: family annotation cannot repair workflow-derived labels.

No new annotation download, domain scan, model fit, test read or AUROC-driven cohort selection was performed in this extension. This report measures **currently evidenced feasibility** and explicitly leaves true full-cohort family support unresolved.

## Reproduction and artifacts

From the repository root with the existing curation environment:

```sh
/tmp/uesolds-venv/bin/python Violet969-PLM_Sol-3589a95/curation/family_feasibility.py
/tmp/uesolds-venv/bin/python Violet969-PLM_Sol-3589a95/curation/render_family_report.py
/tmp/uesolds-venv/bin/python -m unittest discover -s Violet969-PLM_Sol-3589a95/curation -p 'test_family_feasibility.py' -v
```

The audit accepts `--out` for an isolated reproduction; the report renderer reads the default `family_results/`. It consumes only committed audit/annotation tables, so the generated FASTAs previously moved out of the checkout are not needed.

- `primary_analysis_A_manifest.tsv` retains every cleaned-parent row, including unknown and single-class family memberships. Select `cached == True` for the operational primary experiment; retain original split roles.
- `primary_pfam_memberships.tsv` links observed Pfam accessions to original records and homology groups without collapsing them.
- `controlled_analysis_B_memberships.tsv` is an explicit empty table with schema and ≥3/≥5 threshold fields. It derives only from A, with thresholds recomputed separately for `scope=full` and `scope=cached`; use cached scope for the operational experiment. It never replaces A with mixed-label families.
- Coverage, family support, source support, filter summaries, label-reliability tables, plot and input hashes are in `family_results/`.

Checks cover unique records, exact ID/hash/label joins, accession syntax, B being a subset of A and unchanged inputs. Synthetic tests cover both threshold boundaries, source matching, single-class retention, multi-domain/duplicate memberships, and many proteins sharing a single homology group. Existing data, results and test holdout remain untouched.
'''
(P/'FAMILY_FEASIBILITY.md').write_text(text)
print('Wrote',P/'FAMILY_FEASIBILITY.md')
