from audit import *
PRIMARY=['PSI_original','PSI_NetSolP','NESG_raw_scores','NESG_unique_scores','NESG_unique_E45_S45_vs0','NESG_NetSolP','UESolDS_DNASU_proxy']
def table(d):
 def fmt(v):
  if pd.isna(v) or v=='':return '—'
  if isinstance(v,(float,np.floating)):return str(int(v)) if v.is_integer() else f'{v:.2f}'
  return str(v).replace('|','/').replace('\n',' ')
 return '| '+' | '.join(d.columns)+' |\n| '+' | '.join(['---']*len(d.columns))+' |\n'+'\n'.join('| '+' | '.join(fmt(v) for v in row)+' |' for row in d.itertuples(index=False,name=None))+'\n'
def write(name,s):(ROOT/name).write_text(s)
def main():
 c=pd.read_csv(OUT/'dataset_counts.tsv',sep='\t');f=pd.read_csv(OUT/'family_diversity.tsv',sep='\t');o=read(OUT/'dataset_overlap.tsv');supp=read(OUT/'additional_policy_counts.tsv');fs=read(OUT/'family_support.tsv.gz');dist=read(OUT/'family_clan_prevalence.tsv.gz');fold=read(OUT/'five_fold_feasibility.tsv');status=read(OUT/'pET15_status_audit.tsv')
 write('DATASET_COUNTS.md','''# Dataset sizes and outcome policies

All numbers come from the pinned downloaded files. `Soluble`/`insoluble` are deliberately blank for raw ordinal-score cohorts and the usability task. “Unclassified” means no binary outcome was assigned, not that the experimental score is missing. No class balancing was applied.

## Exact counts

'''+table(c[['policy','total','soluble','insoluble','percent_soluble','unclassified','parent_N','excluded','outcome']])+'''
`parent_N` is the input to that row's filtering policy; published releases are initially audited as delivered. NetSolP PSI excludes 990 original SoDoPe records (738 positive, 252 negative). Original Price publication N = 9,644 is not silently substituted for the recovered archive N = 9,703; its exact row membership and binary counts remain unavailable.

## Additional eligibility views

These are separate manifests, not changes to the originals. Consistent-deduplicated PSI excludes all exact-sequence binary conflicts, then counts one lexicographically selected record per agreeing sequence. This cannot resolve differences in experimental constructs removed by NetSolP tag trimming.

'''+table(supp)+'''
## Complete E × S table, as distributed

'''+table(pd.read_csv(OUT/'NESG_E_by_S.tsv',sep='\t'))+'''
E = 0 rows are never treated as experimentally demonstrated insolubility. Nine fractional E records remain unrounded. The 59 E = 0 / S > 0 records need source clarification. Their count matching the publication/mirror size difference does not establish their historical identity.

## Prespecified sensitivity rules

- `E45_S45_vs0`: integer E in {4,5}; positive S in {4,5}; negative S = 0. Primary high-expression extreme-recovery policy.
- `E345_S45_vs0`: additionally admits E = 3; tests support gained at lower expression and potentially noisier recovery assessment.
- `E45_S345_vs0`: additionally calls S = 3 positive; changes the estimand to moderate-or-high recovery versus none.
- `E345_S345_vs0`: combines those relaxations; an exploratory sensitivity, not evidence that moderate scores are unambiguous extremes.
- `E45_S45_vs01`: includes S = 1 in the low-recovery class. These are **low-recovery**, not literally no-soluble-recovery negatives.
- `expressed_ordinal`: retains integer E > 0 with all S values; no binary relabeling. E and S are ordinal, not precise continuous concentrations.

`NESG_exclusion_reasons.tsv` partitions each raw/unique policy's exclusions into noninteger score, E = 0, expression below threshold, and intermediate S. `NESG_repeated_records.tsv` preserves every repeat and conflict. A missing/ambiguous measurement is never used as a negative. The 1,323-record NetSolP release is not identical to literal score-rule filtering: three fractional-E exceptions remain in the released file. The integer-rule view contains 1,320 records (841 positive / 479 negative).

## Direct pET15 outcome cross-check

'''+table(status)+'''
Class mappings are audited against the released raw status fields. pET21 has author classes and sequence files, but not the equivalent complete raw outcome-status export.
''')
 q=f[(f['filter']=='all_GA')&f.policy.isin(PRIMARY)].copy()
 counts=q[['policy','total','annotated','annotation_percent','families','clans','groups','mixed_families','families_proteins_ge3','families_proteins_ge5','families_groups_ge3','families_groups_ge5','families_exclusive_groups_ge3','families_exclusive_groups_ge5']]
 sections=[]
 for pol in ['PSI_original','PSI_NetSolP','NESG_unique_E45_S45_vs0','NESG_NetSolP']:
  g=dist[(dist.policy==pol)&(dist.combined_rank<=10)]
  for level in ['family','clan']:
   a=g[g.level==level];base=a.drop_duplicates('accession').set_index('accession')[['name','description','combined_rank']];pivot=a.pivot(index='accession',columns='label',values=['N','percent']).fillna(0);pivot.columns=[f'{metric}_{int(lab)}' for metric,lab in pivot.columns];base=base.join(pivot).sort_values('combined_rank').reset_index().drop(columns='combined_rank');sections.append(f'### {pol}: {level}\n\n'+table(base))
 long=[]
 for pol in PRIMARY:
  a=fs[(fs.policy==pol)&(fs['filter']=='all_GA')];long.append(dict(policy=pol,families=len(a),singleton_families=int((a.proteins==1).sum()),families_at_most5=int((a.proteins<=5).sum()),largest_family=int(a.proteins.max())))
 save(pd.DataFrame(long),'family_long_tail.tsv')
 write('FAMILY_DIVERSITY.md','''# Pfam diversity and within-family feasibility

Pfam 37.0 / PyHMMER 0.11.2 annotation is complete for every unique input sequence, with existing UESolDS hits reused by exact SHA256. Full-length source sequences are preserved; NetSolP's already-trimmed constructs are annotated as delivered. All reported family statistics are descriptive, not selected by predictive performance.

## Diversity and independent support

'''+table(counts)+'''
Family and clan counts are distinct concepts. “Groups” are connected components of direct MMseqs2 hits at ≥30% identity, ≥80% coverage of each sequence and E ≤ 0.001 across the entire audit union. They are a reproducible operational homology-control criterion, not proof of evolutionary independence. Added sequences can bridge components; therefore UESolDS audit-union group counts differ from its original frozen cluster IDs, which are preserved separately and unchanged.

A protein contributes once to each distinct Pfam accession. Multiple domains may put it in multiple families; family prevalences therefore do not sum to 100%. Counts of proteins in qualifying families use the union of record IDs. ≥3/≥5 **group** thresholds require that many groups containing each class; a mixed-label group can contribute to both class counts. The stricter `exclusive_groups` columns require class-exclusive groups on each side, so the same group cannot supply both classes. Pfam sharing does not establish identical full-protein architecture.

'''+table(q[['policy','multifamily_proteins','proteins_in_groups_ge3','proteins_in_groups_ge5']])+'''
## Annotation ambiguity sensitivity

'''+table(f[f.policy.isin(['PSI_NetSolP','NESG_unique_E45_S45_vs0','NESG_NetSolP'])][['policy','filter','families','mixed_families','families_groups_ge3','families_groups_ge5']])+'''
`no_competing_overlap` removes both hits when different Pfam families from the same clan overlap in aligned coordinates. It is intentionally conservative, not official PfamScan adjudication. `HMM_coverage80` requires ≥80% of the profile aligned. Neither removes a protein from heterogeneous analysis A. The substantial reduction warns against interpreting multiple related overlapping Pfams as multiple independent biological contrasts. All qualifying-family identities and counts are in `family_support.tsv.gz`.

## Homology-aware evaluation feasibility

The prospective diagnostic split assigns each union component to one of five folds by SHA256(`dataset-audit-v1:` + lexicographically smallest sequence hash) modulo 5. Fold 0 is the prespecified evaluation fold; others are development. It is label-blind, without family or class balancing. This is a newly audited design, **not a pristine historical external holdout**. Original UESolDS roles and published NetSolP partitions remain intact in the canonical records.

'''+table(q[['policy','train_groups_ge3_eval_both','train_groups_ge5_eval_both','train_groups_ge3_eval_ge3','train_groups_ge5_eval_ge5']])+'''
The first two columns require ≥3/≥5 development groups per class and at least one evaluation group per class in the same family. The last two require ≥3/≥5 in evaluation as well. A single group per evaluation class establishes minimal evaluability, not a well-powered family-specific estimate. `five_fold_feasibility.tsv` reports every fold rather than selecting the most favorable. Published-PSI fold support is separately recorded in `published_PSI_family_feasibility.tsv`.

For strict NESG, A is biologically heterogeneous and B is feasible for a limited panel; after conservative Pfam overlap removal, only 16 families meet ≥3 groups/class and 3 meet ≥5. NetSolP PSI supports a broader B panel (75 and 21 under the same conservative rule). NetSolP NESG alone is weak for a primary controlled study (2 and 0 after that rule).

Within NetSolP PSI, pET15 has 346 mixed-label families and 48/11 with ≥3/≥5 groups per class; pET21 has 677 mixed families and 116/48 meeting those thresholds. These are same-vector strata, not guaranteed matched experimental batches.

PSI's two vector strata are recorded in `PSI_family_vector_support.tsv`; this permits same-vector contrasts without flattening natural family proportions. NESG's common study pipeline gives better study-level comparability, but complete record-level vector, condition and replicate metadata are absent. Family agreement cannot repair that gap.

## Actual named biological composition

The figures `figures/top20_family.png` and `figures/top20_clan.png` show actual biological names. `top20_named_families_clans.tsv` includes top-20 combined ranks and top-20 ranks within each class (ties included), with exact counts, denominators and percentages. Columns `_0` are insoluble and `_1` soluble. Families absent in one class have zero prevalence there.

'''+ '\n'.join(sections)+'''
The composition differs descriptively: in strict NESG, NAD(P)-binding Rossmann-like PF13450 occurs in 76/1,148 insoluble proteins (6.62%) versus 43/2,614 soluble proteins (1.64%); clan CL0063 occurs in 19.86% versus 8.49%. In NetSolP PSI, protein kinase PF00069 occurs in 116/3,726 insoluble proteins (3.11%) versus 29/7,500 soluble proteins (0.39%). These natural biological associations are retained, not normalized away. Related overlapping hits and broader architectures must be considered before attributing an effect to an individual family.

## Long tail and multidomain accounting

'''+table(pd.DataFrame(long))+'''
`representative_distributions.tsv` provides an optional mutually exclusive 100% distribution: choose the longest aligned domain, then greatest domain bit score, accession and coordinate as deterministic tie-breaks. This is label-independent but loses multidomain information; the main prevalence analysis uses all memberships. No-hit and annotated-without-clan categories are separate. `family_long_tail.png` displays the full ranked frequency distribution.

Length, all 20 amino-acid fractions, charge proxy and mean Kyte–Doolittle hydropathy are summarized by class in `biochemical_summary.tsv`. Unknown-residue sequences retain missing physicochemical summaries rather than silently imputing amino acids. These descriptions are not trained predictors or mechanistic conclusions.
''')
 focus=o[((o.A=='PSI_NetSolP')&o.B.isin(['NESG_raw_scores','NESG_unique_E45_S45_vs0','NESG_NetSolP','UESolDS_development']))|((o.A=='NESG_unique_E45_S45_vs0')&(o.B=='UESolDS_development'))]
 write('DATASET_OVERLAP.md','''# Dataset overlap, conflicting outcomes and independence

These are releases of related experiments, not six independent cohorts. Exact protein hashes, original IDs and direct homology hits are recorded separately. No sequence was modified to manufacture an exact match.

## Direct overlap

'''+table(focus[['A','B','exact_sequences','exact_conflicting_sequences','A_direct_30','B_direct_30','A_direct_90','B_direct_90']])+'''
`A_direct_*` / `B_direct_*` count records on that side with a direct hit to the other at the stated identity and ≥80% bilateral coverage. Full pairwise results, including ≥95%, are in `dataset_overlap.tsv`. Exact conflict counts concern assigned binary outcomes; raw NESG E/S rows have no binary labels and thus their zero binary-conflict count is **not** evidence of agreement. Score conflicts are retained separately.

Original tagged PSI and trimmed NetSolP PSI have zero exact full-sequence overlaps, despite 11,226 accession-matched records with identical labels and fully explained terminal offsets. All 460 UESolDS DNASU-proxy IDs map to original PSI records with agreeing labels, and every cohort sequence is contained in the original construct; none are full-sequence identical. This demonstrates why exact hashes alone understate shared biological/experimental provenance.

PSI NetSolP and strict NESG share 1,138 exact sequences; four have different binary labels. These may reflect different experimental records, constructs, thresholds or updates. Do not overwrite either outcome. At ≥90% identity, additional near-duplicate pairs disagree; their counts are in the overlap table. Common-source history means the remaining proteins are not necessarily independent assays even when sequence leakage is removed.

## Duplicates

'''+table(read(OUT/'duplicate_summary.tsv'))+'''
Binary-conflict fields are not defined as ordinal E/S disagreement. For raw NESG, 150 exact-sequence groups / 316 records have differing E or S scores. Every repeated record is retained in `NESG_repeated_records.tsv`; PSI and cross-dataset binary conflicts are retained in `duplicate_records.tsv` and `exact_conflicts.tsv`. Group support prevents repeat records from masquerading as independent examples.

## Published versus new partitions

The released PSI partitions keep all exact duplicates within one fold. Two cross-fold direct MMseqs pairs pass our local ≥30%/80%-coverage rule (maximum identity 30.7%). This is a different alignment/normalization criterion from the published ≤25% global ggsearch rule, so these hits do not by themselves refute that claim. We verified the delivered folds with our specified criterion; we did not rerun proprietary USEARCH or reproduce the full global-identity graph partitioner.

The 1,323-record NetSolP NESG release has zero exact or ≥90% direct matches to PSI, but 191 records have direct ≥30%/80% hits to curated PSI. Therefore its published external designation is not sufficient under our operational criterion. A conservative removal of components intersecting **original PSI** leaves 1,107 released NESG records (698 soluble / 409 insoluble). This view still preserves the published fractional-score labels; use the integer-score restriction as well for a strict experimental readout.

'''+table(supp[supp.policy.str.contains('no_PSI')])+'''
Component exclusion is more conservative than direct hit exclusion: intermediates from any audit dataset can bridge a component. All eligibility decisions are separate flags/manifests; originals remain intact. Because we inspected these releases for feasibility, they are not described as untouched prospective holdouts. Before future modeling, freeze the selected component-level split and isolate its evaluation labels.

## UESolDS boundary

Only the frozen 13,880 operational proteins were included in the sequence/Pfam/homology union. Exact-hash overlap with the larger original train/validation manifest is separately available in `full_training_validation_hash_overlap.tsv`; no additional UESolDS sequences were extracted. **The locked UESolDS independent test was never accessed.** Existing records, annotation files, embeddings and original cluster assignments were not edited.
''')
 write('DATASET_RECOMMENDATION.md','''# Recommendation for direct-versus-indirect encoding experiments

**Use the strict NESG score cohort as the primary experimental solubility dataset: 3,762 proteins, 2,614 soluble and 1,148 insoluble.** Preserve its full 9,272 unique-sequence score archive alongside it; use the broader expressed-protein ordinal view for a separately declared sensitivity analysis. The strongest labels are the high-expression extremes, not all S = 0 entries in the raw archive.

This is a qualified recommendation for the pinned public NESG/SoluProt archive, not a claim that the missing original 9,644-record release was recovered or that every row's assay protocol is independently verified. It is based on interpretable measured outcomes and sufficient heterogeneous family coverage, not on model performance. No ESM/SAE probes were trained.

## Candidate comparison

'''+table(q[['policy','total','families','mixed_families','groups','families_groups_ge3','families_groups_ge5']])+'''
- **Primary A:** Retain all 3,762 strict-score proteins, including single-class families, no-hit proteins and the natural 69.48% soluble prevalence. This retains 1,383 Pfam families, 243 clans and 1,820 audit homology components. Do not balance family distributions or restrict A to mixed families.
- **Controlled B, from the same A:** 390 families have both labels; 55 have ≥3 groups per class and 28 have ≥5. These cover 883 and 698 proteins respectively (union counts). Only 35/13 families meet the stricter class-exclusive-group thresholds. Under the prespecified split, 21/15 families have ≥3/≥5 development groups/class plus both labels in evaluation. These are candidate contrasts, not 55/28 independent biological mechanisms. Competing Pfam assignments reduce conservative support to 16/3 families; select a small biologically interpretable, architecture-checked panel using development data before any probes.
- **Broader biological comparison:** NetSolP PSI is the stronger option for family diversity and broad B coverage: 2,075 families, 962 mixed families, 181/71 with ≥3/≥5 groups per class. Retain its original DNASU records for auditing. If used for primary binary modeling, flag the 12 conflicting exact-sequence groups and avoid counting repeated sequences as independent. The consistent-deduplicated view has 11,110 proteins (7,402 soluble / 3,708 insoluble), 959 mixed families, and preserves 181/71 group-supported families. It is a credible secondary dataset, but its binary statuses give less quantitative control over expression/recovery than NESG scores.
- **NetSolP NESG alone:** Excellent for reproducing a published external benchmark, insufficient as the main controlled-family dataset: 13/2 group-supported families, falling to 2/0 under conservative overlap handling. The three fractional-E exceptions also need explicit treatment. Its strict integer-rule view is 1,320, not 1,323.
- **DNASU proxy in UESolDS:** 460 records / 409 positive / 51 negative; no family meets ≥3 proteins of each class. It cannot substitute for full PSI or provide the desired B analysis.

## Ordinal NESG analysis

The 7,200 unique-sequence records with integer E > 0 provide an ordinal S analysis spanning intermediate outcomes. Keep E explicit, report expression strata, and do not call S precise continuous solubility. The alternative E ≥3/extreme-S policy yields 4,436 proteins (2,992 / 1,444), 1,522 families and 66/31 group-supported families. Its larger support comes with a changed expression-conditioned population, so it remains a prespecified sensitivity rather than an automatic replacement for the primary rule.

## Training/development versus external evaluation

Do not simply train on PSI and call the full NESG archive independent: the strict cohorts share 1,138 exact sequences and extensive homology. Conversely, fitting on full NESG prevents treating all PSI as an independent external evaluation. A conservative PSI-development / NESG-external design can retain 1,176 strict NESG proteins (741 / 435) after removing components shared with original PSI, but has only seven families with ≥3 groups/class. It supports an aggregate external check better than broad within-family validation. See the explicit eligibility manifest and overlap report before choosing direction.

Do not pool PSI and NESG binary labels into one undifferentiated training target. Their related history, duplicate experiments, vector differences and binary-versus-score definitions would recreate the provenance problem we are trying to address. They can be compared as separate outcome strata in a shared annotation inventory, with shared sequence components assigned consistently.

## Remaining information and smallest useful next steps

1. Obtain the original 9,644-ID E/S export (or a dated mapping to the public 9,703 records), including the meaning of fractional E and E = 0 / S > 0 entries. The present strict integer/high-E policy avoids those anomalies but cannot establish historical membership.
2. Recover replicate/construct/condition details for the four exact PSI–strict-NESG binary disagreements and the 12 conflicting PSI sequence groups. This small targeted audit is more informative than collecting more embeddings.
3. For a development-selected B panel, inspect overlapping Pfam calls and complete architectures, and verify assay/construct comparability. Freeze a homology-aware evaluation protocol before training. Our saved five-fold feasibility split is descriptive and must not be portrayed as the original untouched external test.

The evidence supports a useful experimental-solubility representation study with a conditional target: **soluble recovery of a recombinant construct in the documented E. coli production setting**. It does not turn sequence associations into causal computations. SAE interventions on the model would test effects inside the model; experimental changes in protein solubility would still require wet-lab validation.
''')
if __name__=='__main__':main()
