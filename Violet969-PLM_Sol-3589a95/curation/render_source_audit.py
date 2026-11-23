#!/usr/bin/env python3
"""Render new audit reports from saved results; does not fit models."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
ROOT=Path(__file__).resolve().parent;OUT=ROOT/'source_audit_results'
def read(n):return pd.read_csv(OUT/n,sep='\t',keep_default_na=False)
def md(d,cols=None):
 if cols:d=d[cols]
 def cell(v):
  if isinstance(v,(float,np.floating)):return f'{v:.4f}' if np.isfinite(v) else 'undefined'
  return str(v).replace('|',' / ')
 return '\n'.join(['| '+' | '.join(map(str,d.columns))+' |','| '+' | '.join(['---']*len(d.columns))+' |']+['| '+' | '.join(cell(x) for x in row)+' |' for row in d.itertuples(index=False,name=None)])
def ci(r,metric='auroc'):
 try:return f'{float(r[metric]):.3f} [{float(r[metric+"_low"]):.3f}, {float(r[metric+"_high"]):.3f}]'
 except (ValueError,TypeError):return 'undefined'
def write(n,parts):(ROOT/n).write_text('\n\n'.join(parts)+'\n')

def diversity():
 s=read('source_diversity.tsv');c=s[s.split=='combined'];loss=read('leave_source_out_diversity.tsv');policy=read('policy_comparison.tsv');overlap=read('source_pair_overlap.tsv')
 counts=s[s.split!='combined'][['center','split','n','n0','n1','groups','families','clans']]
 cov=c[['center','n','pfam_annotated','taxonomy_assigned','families','clans','mixed_families','family_source_pairs_both_splits_ge3']].copy()
 cov['Pfam_percent']=100*cov.pfam_annotated/cov.n;cov['taxonomy_percent']=100*cov.taxonomy_assigned/cov.n
 desc=read('source_label_descriptors.tsv');geom=read('source_biochemical_neighbour_loss.tsv')
 pair=overlap[(overlap.source_a=='JCSG')&(overlap.source_b=='MCSG')]
 ref=pd.read_csv(ROOT/'pfam_annotation/Pfam-A.clans.tsv',sep='\t',header=None,names=['pfam','clan','clan_name','family_name','description'],keep_default_na=False)
 dist=read('source_family_clan_distribution.tsv');names=ref.set_index('pfam').description.to_dict();names.update(ref[ref.clan!=''].drop_duplicates('clan').set_index('clan').clan_name.to_dict())
 dist['name']=dist.accession.map(names);dist.to_csv(OUT/'source_family_clan_distribution.tsv',sep='\t',index=False,float_format='%.17g')
 tops=dist[dist.level=='family'].sort_values(['center','proteins','accession'],ascending=[True,False,True]).groupby('center').head(3)
 write('SOURCE_DIVERSITY.md',[
 '# Source contribution and biological diversity — 2026-11-12',
 '**The operational cohort remains 13,880 proteins, 6,043 Pfam families, 557 clans and 13,716 detected homology groups.** Removing a source can remove real biological diversity even when its solubility labels are unsuitable for a primary property readout. All measurements below refer to this cached cohort, not the full cleaned pool.',
 '## Scope and counting',
 'Pfam 37.0 GA-passing memberships and existing MMseqs2 groups are reused without changes. Each protein counts once in cohort totals; protein–family and protein–clan pairs are deduplicated. Multidomain proteins can occur in several membership rows, so family/clan percentages need not sum to 100%. Family membership is not full architecture. Groups use the existing 30% identity, 80% bilateral coverage, E≤0.001 criterion; these are operational diversity units, not guarantees of evolutionary or statistical independence. Training and validation groups are disjoint.',
 '## Counts by source and original split',md(counts),
 '## Biological coverage and mixed-label support',md(cov),
 'The exact species-consensus taxonomy annotation is independent of the centre proxy. Low coverage is substantial and source-dependent: only 1/460 DNASU-proxy proteins has a consensus assignment in the current exact full-sequence procedure. This is an annotation limitation, not evidence of low organism diversity. Ambiguous taxonomic matches, engineered taxa and unmatched proteins remain distinct in [taxonomy_coverage.tsv](source_audit_results/taxonomy_coverage.tsv).',
 '## Marginal contribution: remove one source, retain everything else',
 'Losses below are observed set differences on the combined cohort. A family can remain present but cease to have both labels, hence loss of mixed-family support differs from loss of a family. Family × source opportunities disappear with the removed source even when other sources still contain the family.',
 md(loss,['center','removed_n','unique_families_lost','unique_clans_lost','groups_lost','lost_mixed_families','lost_mixed_family_source_pairs','lost_family_source_pairs_both_splits_ge3']),
 '![Marginal family contributions](source_audit_results/marginal_family_diversity.png)',
 'PDB_serial contributes **1,494 exclusive families, 63 exclusive clans and 2,845 groups absent from all other sources**; it also supplies one side of 293 otherwise-lost mixed-label families. This biological value supports retaining its sequences for concept analysis; it does not turn structure-selected positives into standardized solubility measurements. RefSeq adds 56 exclusive families despite being positive-only. These observations rule out class imbalance alone as an exclusion criterion.',
 '## JCSG and MCSG: complementary, not interchangeable',md(pair),
 'Across both splits, their family union is **2,503**, with **449 shared families**, **1,098 JCSG-only** and **956 MCSG-only** relative to each other. Only **10 detected homology groups** span the two sources, and no exact sequence is shared. The previously reported 322-family intersection used training only; the population here includes validation and is explicitly different.',
 'Together they retain 4,774 proteins, 4,743 groups and 582 mixed-label families, but only ten families have ≥3/class in both splits; requiring the same centre leaves four family/centre pairs. Those pairs involve related methyltransferase Pfams and are not four independent biological mechanisms. Under the existing no-clan-competition annotation sensitivity, previous analysis found no such robust pairs. These centre choices do not fix label validity.',
 '## Biochemical coverage, composition and construct observations',
 'The complete [source × split × label descriptor table](source_audit_results/source_label_descriptors.tsv) reports length, charge proxy, KD hydropathy, entropy, all 20 amino-acid fractions, terminal His/cleavage/linker flags, non-M starts, repeated IDs, annotation coverage and existing APR/disorder/TM predictions. Missing predictions remain missing; a motif flag is an observation, not proof of an expression tag or truncated construct. Means and 5th/50th/95th percentiles are supplied without selecting features by label association.',
 'A prespecified 3D grid uses training-fitted quintiles of log full length, charge proxy and KD hydropathy. No single source removal eliminates an occupied coarse cell: other sources span that coarse range. This is a limitation of a low-resolution coverage measure, not evidence that distributions or biological diversity are unchanged.',
 'A complementary 21-dimensional comparison uses full-sequence AA fractions and log length, standardized on training only. For every protein in the removed source, measure the nearest retained protein outside its MMseqs group, compared with its nearest different-group neighbour before removal. The table describes geometry in this chosen descriptor space, not a universal biological distance.',
 md(geom,['center','n','baseline_nearest_other_group_median','retained_nearest_other_group_median','increase_p95','fraction_distance_increased']),
 'SGPP shows the largest affected fraction (about 55%), consistent with a distinctive sampled biochemical distribution. Zero median distance increase in several sources does not imply no tail losses: source-specific families and upper-tail neighbour distances still change. The feature axes and scale are fixed without labels.',
 '## Named family distribution by source',
 'The three most prevalent memberships per source are shown below. All family/clan counts and percentages, including long tails, are in [source_family_clan_distribution.tsv](source_audit_results/source_family_clan_distribution.tsv). Selection here uses source-wise frequency for description, not solubility association.',
 md(tops,['center','accession','name','proteins','denominator']),
 '## Reproduction and limitations',
 'Run `source_validity_audit.py --stage diversity`, then the record-recovery scripts and `source_audit_supplement.py`, then `render_source_audit.py`. Fixed protocol and input hashes are in source_audit_results/. Coarse cells and nearest-neighbour measurements are descriptive; source-removal estimates mix source size with unique contribution and are not a ranking of source quality. Counts are observed in this development sample, not the entire originating programme. No annotation, clustering or embedding extraction was rerun.'
 ])

def shortcuts():
 s=read('source_prediction_metrics.tsv');b=read('solubility_metrics.tsv');d=read('solubility_paired_differences.tsv')
 fam=read('source_within_family.tsv');family=read('solubility_family_strata.tsv');cross=read('cross_centre_uncertainty.tsv')
 source_table=[]
 for model,g in s.groupby('model',sort=False):
  row={'model':model}
  for key in ['accuracy','balanced_accuracy','macro_f1','macro_ovr_auroc','macro_ap','logloss']:
   r=g[g.metric==key].iloc[0];row[key]=f'{r.value:.3f} [{r.low:.3f}, {r.high:.3f}]'
  source_table.append(row)
 pooled=b[(b.protocol=='all_train')&(b.evaluation=='ALL')].copy()
 base=[]
 for center,g in b[b.protocol=='all_train'].groupby('evaluation',sort=False):
  r=g.iloc[0];row=dict(evaluation=center,n=r.n,n0=r.n0,n1=r.n1,prevalence=r.prevalence)
  for model in ['composition','ESM_layer8']:
   r=g[g.model==model].iloc[0];row[model+'_AUROC_CI']=ci(r);row[model+'_AP']=r.ap
  base.append(row)
 conditional=fam.groupby('model',sort=False).agg(families=('pfam','nunique'),accuracy=('accuracy','mean'),balanced_accuracy=('balanced_accuracy','mean'),macro_ovr_auroc=('macro_ovr_auroc','mean'),logloss=('logloss','mean')).reset_index()
 sourceinc=read('source_paired_differences.tsv')
 transfer=cross[cross.evaluation.isin(['JCSG','MCSG'])].copy();transfer['AUROC_95CI']=transfer.apply(ci,axis=1)
 loo=b[b.protocol.str.startswith('without_') & b.evaluation.eq('ALL')].copy();loo['AUROC_95CI']=loo.apply(ci,axis=1)
 sf=family.groupby(['kind','model']).agg(strata=('stratum','nunique'),mean_auroc=('auroc','mean'),mean_ap=('ap','mean')).reset_index()
 delta=d[(d.protocol=='all_train')&(d.metric=='auroc')]
 write('SOURCE_SHORTCUTS.md',[
 '# Sequence-accessible source information and solubility sensitivity — 2026-11-12',
 '**Source-associated information is accessible from sequence and frozen ESM activations.** It does not follow that the solubility readout uses a particular source cue or that a centre effect is technical rather than biological. The original-record audit independently finds weak evidence for many negative labels.',
 '## Fixed evaluation protocol',
 'The models use exactly 9,880 training and 4,000 validation records, with original roles, labels and detected homology groups preserved. No detected group, exact sequence or reused identifier spans these splits. This protects against detected leakage at the existing MMseqs threshold, not all remote homology or domain overlap. Validation has already been inspected repeatedly; all results are exploratory development evidence. The independent test remains locked.',
 'Source prediction uses all 17 frozen header-derived centre/repository classes. Fifteen occur in validation (BSGC and CESG do not). Composition is 20 AA fractions plus log length on the first 1,022 residues; frozen ESM-2 layer-8 activations are cached 1,280-dimensional means on the same visible residues. Pfam uses full-sequence multi-hot memberships, a training-only family vocabulary and an explicit no-hit indicator; it has a modest full-length information advantage over ESM on long proteins. Validation-only families are ignored as unseen features, not counted as absent biological annotation.',
 'L2 logistic regression, C=1, lbfgs, tolerance 1e-4, maximum 5,000 iterations, no class weights or hyperparameter search. Dense sequence inputs are standardized using only the relevant fitting training subset. Encoders see training only. All fits converged (source: 74–425 iterations; solubility transfer: 7–209). No ESM weights were loaded or updated. The existing all-source solubility predictions are reused, including their saved precision; new fits are limited to training-source sensitivities.',
 'Approximate 95% intervals use 200 resamples of validation homology groups, keeping group members together. Paired differences use identical resamples. These are conditional on fitted models and do not include training-fit or label uncertainty; there is no multiplicity adjustment. AUROC is undefined for single-class strata. AP is interpreted with observed positive prevalence, and is trivially 1 for all-positive strata. For tiny minorities—especially SSGCID with one negative—bootstrap intervals cannot support reliable population inference; undefined bootstrap draws are counted and excluded.',
 '## Source prediction',md(pd.DataFrame(source_table)),
 'The training-prior predictor always selects the most frequent training source (PDB_serial): 18.7% validation accuracy. Balanced accuracy is mean recall over the 15 represented validation classes; macro F1 uses those classes too. Macro one-vs-rest AUROC/AP omit absent validation classes, and log loss retains all 17 probability columns. The explicit prior model accounts for unequal class sizes; chance accuracy is not assumed to be 1/17.',
 'ESM improves classification accuracy over composition by **7.80 percentage points** (paired 95% CI 6.25–9.33), and balanced accuracy by about 8.43 points. However, ESM has much worse log loss (4.269 versus 2.117), indicating poor confidence calibration under this fixed high-dimensional probe and source shift. The model has not been calibrated or tuned after seeing validation. AUROC/accuracy cannot conceal that failure.',
 'Adding ESM to composition improves accuracy over composition by **7.72 points** (6.02–9.46). This demonstrates incremental linear predictive information under the specified models; it is not a causal adjustment or proof of independence from all nonlinear composition effects. Source-target identity itself mixes laboratory, repository and naming conventions.',
 '[Per-centre precision, recall, support, one-vs-rest AUROC/AP and intervals](source_audit_results/source_prediction_per_centre.tsv), [confusion matrices](source_audit_results/source_confusion.tsv), [paired increments](source_audit_results/source_paired_differences.tsv) and all probability tables are retained.',
 '## Source distinguishability within families',md(conditional),
 'There are **27 qualifying families** with at least two sources having ≥3 records in both training and validation. Restrict candidate sources to those supported in that family and renormalize every model equally. A family-specific training-frequency prior provides a stronger baseline than a global prior. These are equal-family descriptive means; memberships overlap, and they are not independent replicates or a validation-selected optimal subset. Exact per-family counts and metrics are in source_within_family.tsv and source_within_family_support.tsv.',
 'Per-family intervals are saved in [source_within_family_uncertainty.tsv](source_audit_results/source_within_family_uncertainty.tsv). They resample validation homology groups within each fixed stratum; draws missing a candidate source are excluded and counted. There is no independence claim across overlapping families and no macro-mean confidence interval.',
 'Source discrimination survives family conditioning: composition and ESM both show conditional ranking signal. ESM does not dominate every metric: mean conditional one-vs-rest AUROC is about 0.662 for ESM and 0.681 for composition. Sharing one Pfam does not hold all domains, architecture, organism or chemistry fixed. These results do not identify an experimental artefact.',
 '## Solubility: pooled and within-centre, same fitted predictions',md(pd.DataFrame(base)),
 'Pooled ESM AUROC is 0.6769 versus 0.6565 for matched composition; the paired difference is 0.0204 [0.0030, 0.0360]. The metadata-only source model is 0.7791 [0.7664, 0.7915]. These compare the same 4,000 records. They do **not** share an estimand with a five-centre mean, and subtracting a pooled score from a within-centre average is not a causal source-reliance estimate.',
 'Per-centre paired ESM minus composition AUROC:',md(delta,['evaluation','n','difference','low','high','valid_bootstraps']),
 'ESM has higher ranking performance than composition within JCSG and MCSG, but is weaker in SGPP. SECSG improves relative to an inverted composition signal while remaining near/below chance. Centre-wise variation, tiny strata and possible label-generation differences prevent interpreting these as universal biochemical laws. A constant metadata source model scores 0.5 within a centre by construction; that does not exclude finer source-related sequence signals.',
 '## Cross-centre training',md(transfer,['protocol','model','evaluation','n0','n1','AUROC_95CI','ap','prevalence']),
 'Each row trains only on the named centre’s original training proteins and evaluates the named, unchanged validation stratum. JCSG→MCSG and MCSG→JCSG are separate transfer tasks, not paired population comparisons. Their results do not upgrade either centre’s label validity. The complete matrix includes all validation sources, with intervals in cross_centre_uncertainty.tsv.',
 '## Leave-one-source-out training: fixed 4,000-protein validation population',md(loo,['protocol','model','n','AUROC_95CI','ap']),
 'All training-source omissions above evaluate the **same** validation proteins; no removed source is also silently removed from evaluation. Separate omitted-centre-only results and paired AUROC/AP changes are saved. Removing PDB training lowers pooled ESM AUROC to about 0.634; this changes training size, class composition and biology simultaneously. It is a sensitivity result, not a causal estimate of a PDB shortcut or a reason to discard PDB proteins.',
 '## Family-overlapping versus family × source solubility evaluation',md(sf),
 'These use global all-source fitted predictions, not separate family-tuned fits. Training support is ≥3/class, and the reported strata also have ≥3/class in validation: 69 Pfam families versus only four family/source cells. The source model still ranks labels within families (mean AUROC 0.776), showing that family matching alone leaves source differences. Within family/source cells the source model is constant, while sequence probes can vary. Only three related methyltransferase accessions supply the four cells; this tiny, overlapping set cannot establish a general within-family mechanism. [Individual metrics and group-bootstrap intervals](source_audit_results/solubility_family_strata.tsv) accompany the means. No family associations were erased or rebalanced.',
 '## Interpretation and feasibility',
 'High pooled AUROC does not validate labels. Near-chance within-centre performance does not establish that ESM lacks solubility information. Predictable source information does not show that the solubility probe uses it. Predictive readouts or overlapping probe weights do not establish causal computation; ablation must distinguish a change in the probe from a biologically meaningful output change.',
 'The cached cohort suffices for this audit. Expanding embeddings is not the next priority: previous analysis identified 59,368 extra clean training proteins, but their uniform Pfam annotations are unknown and unchanged validation limits controlled contrasts. First recover verified construct-specific outcomes. Mean-pooled caches suffice for these logistic probes, but future SAE interventions require appropriate token-level forward passes and reconstruction controls; no SAE analysis was performed here.',
 'Reproduce with `source_validity_audit.py --stage probes`, `source_audit_uncertainty.py`, then `render_source_audit.py`. Protocol, fit diagnostics, predictions and input hashes are saved in source_audit_results/. The procedure does not rerun annotations, clustering or ESM extraction.'
 ])

def selection():
 p=read('policy_comparison.tsv');combined=p[p.split=='combined'];splits=p[p.split!='combined']
 write('DATASET_SELECTION.md',[
 '# Primary dataset decision — 2026-11-12',
 '**Recommendation: policy C, one canonical heterogeneous cohort with separate evidence and analysis-eligibility flags. Retain all 13,880 proteins for biological-concept work, but do not yet designate any complete source as a validated primary binary solubility cohort.** Explicit DNASU outcomes are the first recovery priority; currently no reliable, comparable binary property core has been established.',
 'This differs from treating all-source A as a scientifically validated solubility benchmark. A remains useful for exploratory published-label prediction and broad biological concepts. A good pooled score, source matching or restriction to a single centre cannot repair unsupported negative labels.',
 '## Candidate policies: measured coverage',md(combined,['policy','n','n0','n1','sources','families','clans','groups','mixed_families','families_both_splits_ge3','family_source_pairs_both_splits_ge3','families_lost','groups_lost']),
 '“Both splits ≥3” requires at least three proteins per label independently in training and validation. Pfam families, clans and homology groups are different units. A “mixed” family in the pooled cohort is weaker than a supported held-out within-family comparison. C_verified_property_core = 0 means that no comparable binary assay cohort was established; it does not deny the two recovered MCSG soluble statuses or the explicit DNASU repository statuses.',
 'Original split/class counts for each policy:',md(splits,['policy','split','n','n0','n1']),
 'Exact source composition for every policy is in [policy_source_counts.tsv](source_audit_results/policy_source_counts.tsv); the [canonical manifest](source_audit_results/policy_manifest.tsv) retains original IDs, sequence hashes, source, split, label and homology group and adds eligibility/evidence columns only. Existing datasets and prior results are not overwritten.',
 '## A — all-source heterogeneous property cohort',
 'Retains all biological diversity and observed family–label correlations. Suitable for descriptive concept associations, source diagnostics and exploratory prediction of **published UESolDS labels**. Unsuitable as a verified solubility benchmark: PDB positives are structure-selected, TargetTrack negatives can reflect selection/cloning, and eSOL linkage is unresolved. Interventions changing this label readout would not establish a direct biochemical solubility mechanism.',
 '## B — selected-source heterogeneous property cohort',
 'DNASU is the only currently identified repository stratum with directly recoverable explicit positive and negative outcome fields. Its 460 proteins (409 positive, 51 negative), 405 families and 456 groups offer a candidate, not certified, property subset. Twelve audited statuses agree with released labels, but construct reconciliation is incomplete for all twelve and assay comparability is unresolved. It loses 5,638 families relative to A and has no ≥3/class family contrast in both splits. Thus it is not yet a substitute for the heterogeneous direct-versus-indirect study, even if its published labels are more plausible.',
 'No alternative selected-source list has demonstrated better negative-outcome evidence. Excluding PDB/RefSeq solely because they are positive-only sacrifices diversity without repairing retained negatives. Selection should follow recoverable measurement definitions, never an AUROC ranking.',
 '## C — canonical cohort plus evidence-qualified views (preferred)',
 'Keep every original protein and family association. `concept_eligible=True` for all records; annotation-specific availability determines which concepts can be analyzed without imputing biology. `published_label_exploratory_only=True` preserves the released outcome as a descriptive target. `verified_property_eligible=False` until comparable construct-linked binary evidence is adjudicated. Separate flags record exact TargetTrack trial matches, explicit soluble stages and matching DNASU statuses. These flags distinguish recovery from certification.',
 'For primary property labels, prioritize reconciled DNASU Protein_Soluble/Tested_Not_Soluble outcomes under a documented assay definition. Add JCSG, MCSG or other centres only when individual assayed constructs have explicit outcomes—including actual negatives—and compatible conditions. Do not use selected/cloned states, absent scores, centre names or successful structure deposition alone as primary binary truth. PDB and other extension proteins remain useful for family, structural and biochemical concept probing.',
 'This is one dataset with transparent analysis eligibility, not a balanced replacement dataset. The primary label subset may ultimately have a narrower scope than the concept cohort; report that scope rather than claiming direct generalization to all proteins. If family × source support stays sparse after evidence recovery, controlled analysis B must remain exploratory or require new measured examples.',
 '## D — single or limited centres',
 'JCSG and MCSG provide complementary families and ten shared detected groups. Their union retains 2,503/6,043 families, losing 3,540, with four supported family/source cells. Neither has recovered explicit negative solubility measurements in the audited sample. The eight training-supported centres retain 10,058 proteins but still lose 1,792 families and do not increase the four controlled cells. Source comparability improves only at the proxy level; assay and negative-label comparability remain unproven. These are diagnostic views, not preferred primary cohorts.',
 '## Suitability for the planned research',
 '| Use | Current defensible scope | What remains necessary |\n|---|---|---|\n| Biological-concept probes | Heterogeneous canonical cohort with sequence/Pfam and documented predictor annotations | Respect annotation missingness, domain overlap and model-derived evidence |\n| Correlational published-label probes | Exploratory all-source results with matched centre/family strata | Explicitly name the outcome as published UESolDS label, not validated solubility |\n| Primary solubility readout | Not yet established | Construct-linked positive **and negative** outcomes with a consistent assay definition |\n| SAE intervention feasibility | Can test technical reconstruction/ablation machinery separately | Valid property readout, token-level activations, matched controls, model-output effects; probe change alone is not causal biochemical evidence |\n| Controlled within-family conclusions | Limited exploratory views | More independent families/groups with reliable two-class outcomes and fixed evaluation support |',
 '## Final decision and smallest next step',
 '**We cannot yet certify a reliable, biologically heterogeneous binary solubility cohort from the available UESolDS evidence.** We can preserve a useful heterogeneous concept cohort and now have a concrete route to stronger labels. Start with the fixed twelve DNASU records already audited: reconcile the original fusion sequences and obtain their expression/solubility score definitions and protocol versions. Then extend that exact rule to all 51 DNASU negatives and corresponding positives. Ask for the pre-cleaning UESolDS sequence-to-assay crosswalk or the original NESG/DNASU assay export; the current terminal differences are directly measurable, not hypothetical.',
 'For structural-genomics negatives, request centre LIMS results for an actual soluble-expression test on the matching construct; a TargetTrack selection/cloning record cannot resolve the question. If those results cannot be recovered, a small matched experimental panel with measured soluble fractions is the next reliable option. No filtering, class balancing, label rewriting or new embeddings can substitute for that missing evidence.',
 'Supporting evidence: [source validity](SOURCE_AUDIT.md), [diversity](SOURCE_DIVERSITY.md), [sequence-source and property probes](SOURCE_SHORTCUTS.md). The independent test remains unopened, all original records and previous artifacts are preserved, and cohort policy was not selected to maximize AUROC.'
 ])

if __name__=='__main__':diversity();shortcuts();selection()
