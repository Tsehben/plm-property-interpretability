#!/usr/bin/env python3
"""Render uniform Pfam results; no model scores consumed."""
from pathlib import Path
import json
import pandas as pd
P=Path(__file__).resolve().parent;R=P/'pfam_annotation'
def table(d):
    def fmt(x):
        if pd.isna(x):return 'NA'
        if isinstance(x,float):return f'{x:.4g}'
        return str(x)
    return '| '+' | '.join(d.columns)+' |\n| '+' | '.join(['---']*len(d.columns))+' |\n'+'\n'.join('| '+' | '.join(fmt(x) for x in row)+' |' for row in d.itertuples(index=False,name=None))+'\n'
a=pd.read_csv(R/'annotated_cohort.tsv',sep='\t',keep_default_na=False)
s=pd.read_csv(R/'feasibility_summary.tsv',sep='\t');main=s[s.policy=='GA_all']
dom=pd.read_csv(R/'domains_flagged.tsv',sep='\t',keep_default_na=False)
meta=json.loads((R/'search_run.json').read_text())
family=pd.read_csv(R/'family_support.tsv',sep='\t');sources=pd.read_csv(R/'family_source_support.tsv',sep='\t')
assigned=a.pfam_status.eq('assigned');n=int(assigned.sum());f=int(main.iloc[0].families);mixed=int(main.iloc[0].mixed_families)
text=f'''# Uniform Pfam annotation and controlled-analysis feasibility

All **13,880** frozen cached-cohort proteins were searched against **21,979 Pfam 37.0 profiles**, using pinned **PyHMMER 0.11.2** and the per-profile sequence/domain gathering thresholds. All input sequences, labels, original split roles and homology-cluster assignments were preserved. No protein was selected by label, source, prior annotation or model performance. The locked test was never read.

**Coverage:** {n:,}/13,880 proteins ({100*n/13880:.2f}%) have at least one GA-passing domain, spanning **{f:,} observed Pfam accessions**. **{mixed:,} families contain both labels** across development. This replaces the pilot-only feasibility limitation; the old pilot report remains a historical snapshot.

## Primary cohort and interpretation

Keep **A: all 9,880 training + 4,000 validation proteins** as the heterogeneous primary experiment, including single-class families and searched/no-hit records. Use B as a controlled evaluation on family strata derived from that same A cohort. Do not rebalance family distributions or train the primary property probe only on B. The cohort was frozen before this annotation run; no AUROC was computed or consulted.

'''
for _,r in main.iterrows():
    text+=f"At **≥{int(r.k)} proteins per class in training**, {int(r.training_eligible):,} families qualify; {int(r.training_eligible_val_both):,} also have both labels in validation. Of these, {int(r.both_splits_threshold):,} meet ≥{int(r.k)} per class separately in **both splits**, and {int(r.both_splits_k_exclusive_groups):,} also meet that threshold using class-exclusive detected homology groups in each split. Same-source matching leaves {int(r.same_source_val_both_pairs):,} training-supported family × source pairs with both validation classes, including {int(r.same_source_both_splits_threshold_pairs):,} pairs meeting the threshold in both splits.\n\n"
if main.both_splits_threshold.max()>0:
    text+='**B is now feasible as an exploratory within-Pfam contrast**, with the strongest candidates requiring support in both splits and within the same source. Feasibility is not a power calculation or proof of reliable biochemical solubility labels. A shared Pfam domain does not establish identical whole-protein architecture, organism, construct or assay.\n\n'
else:
    text+='**B remains limited:** although training-supported families may be minimally evaluable, none meets the requested class-support threshold independently in both published splits. Do not treat pooled support as sufficient held-out support.\n\n'
conservative=s[s.policy=='GA_without_clan_competition'].set_index('k')
text+=f"**The source-matched version is much less secure:** after removing competing clan-overlap assignments, only {int(conservative.loc[3,'both_splits_threshold'])} families retain ≥3 per class in both splits and {int(conservative.loc[5,'both_splits_threshold'])} retain ≥5; no same-source pair reaches either threshold in both splits. This sensitivity is conservative rather than official clan resolution, but it shows that the nominal same-source support is annotation-dependent.\n\n"
text+='''## Annotation coverage and missingness

All records were attempted successfully. `searched_no_GA_hit` means no domain met the pinned profile thresholds; it does not prove absence of a biological family. There are no unqueried records and no failed searches hidden as no-hit outcomes. Full sequences were annotated even when the cached ESM representation sees only the first 1,022 residues.

'''
rows=[]
for (split,y),g in a.groupby(['split','label']):
    nh=g.pfam_status.eq('assigned').sum();rows.append(dict(split=split,label=y,n=len(g),assigned=int(nh),no_hit=int(len(g)-nh),coverage_percent=round(100*nh/len(g),2)))
text+=table(pd.DataFrame(rows))
text+='\nSource × label coverage (training and validation pooled here; split-specific counts are saved):\n\n'
rows=[]
for center,g in a.groupby('center'):
    row={'source':center}
    for y in [0,1]:
        v=g[g.label==y];nh=v.pfam_status.eq('assigned').sum()
        row[f'label{y}: assigned/total']='—' if not len(v) else f'{nh}/{len(v)} ({100*nh/len(v):.1f}%)'
    rows.append(row)
text+=table(pd.DataFrame(rows))
text+='\n![Uniform annotation coverage](pfam_annotation/coverage.png)\n\n'
text+=f"There are {len(dom):,} retained domain instances; repeat copies count once toward each protein–family membership. {dom[dom.clan_competing_overlap].raw_uid.nunique():,} proteins have competing overlapping different-family hits within the same clan. Primary annotation retains GA-passing evidence; a conservative sensitivity removes both competing alternatives. Clans never substitute for Pfam families.\n\n"
text+='''## Family and homology-group support

Pooled-development counts below describe available family diversity; they do **not** select B families. B eligibility uses training counts, and validation support is reported separately. Protein totals are unions across qualifying families, avoiding double-counting multidomain proteins.

'''
text+=table(main[['k','pooled_qualifying','pooled_eligible_proteins','pooled_eligible_groups','pooled_eligible_groups0','pooled_eligible_groups1']].rename(columns={'k':'threshold','pooled_qualifying':'families','pooled_eligible_proteins':'proteins','pooled_eligible_groups':'all groups','pooled_eligible_groups0':'label 0 groups','pooled_eligible_groups1':'label 1 groups'}))
text+='\nTraining-selected families with both validation classes:\n\n'
text+=table(main[['k','training_eligible','training_eligible_val_both','both_splits_threshold','eval_eligible_unique_proteins','eval_eligible_val_proteins','eval_eligible_groups0','eval_eligible_groups1']])
text+='''
`family_support.tsv` gives each family's per-class protein counts, distinct groups per class, groups shared across labels, class-exclusive groups, source counts and split support. Group counts are not summed across families. A group appearing under both labels is shared, not two independent observations. Detected MMseqs2 groups are operational independence units, not guarantees against remote homology; their existing 30%-identity/80%-bilateral-coverage definition is unchanged. No group spans training and validation in this cohort. Qualifying family counts based on class-exclusive groups are a conservative secondary support check.

Examples below have ≥5 training proteins per label and ≥3 validation proteins per label, ordered by validation minority-class count, **not model performance**. These are illustrations, not a selected replacement for A or the complete B membership list.

'''
t=family[(family.policy=='GA_all')&(family.split_scope=='train')];v=family[(family.policy=='GA_all')&(family.split_scope=='validation')]
x=t.merge(v,on='pfam',suffixes=('_train','_val'));x=x[(x.n0_train>=5)&(x.n1_train>=5)&(x.n0_val>=3)&(x.n1_val>=3)].copy();x['minority_val']=x[['n0_val','n1_val']].min(axis=1)
x=x.sort_values(['minority_val','pfam'],ascending=[False,True]).head(12)
text+=table(x[['pfam','n0_train','n1_train','groups0_train','groups1_train','n0_val','n1_val','groups0_val','groups1_val']])
text+='''
## Same-source contrasts and label reliability

Same-source eligibility requires the threshold in the **same Pfam × source stratum in training**, then evaluates validation support in that stratum. Source counts below refer to unique proteins within the union of such minimally evaluable strata; a protein may be in multiple family strata, but contributes once per source/split/label cell.

'''
bs=pd.read_csv(R/'B_source_distribution.tsv',sep='\t');bs=bs[(bs.threshold==5)&(bs.contrast=='same_source_validation_both')]
if len(bs):
    b=bs.groupby(['center','label'])[['n_proteins','n_groups']].sum().reset_index();text+=table(b)
text+='\nSame-source pairs meeting ≥3 proteins per class in both original splits:\n\n'
st=sources[(sources.policy=='GA_all')&(sources.split_scope=='train')]
sv=sources[(sources.policy=='GA_all')&(sources.split_scope=='validation')]
sx=st.merge(sv,on=['pfam','center'],suffixes=('_train','_val'))
sx=sx[(sx.n0_train>=3)&(sx.n1_train>=3)&(sx.n0_val>=3)&(sx.n1_val>=3)]
text+=table(sx[['pfam','center','n0_train','n1_train','groups0_train','groups1_train','n0_val','n1_val','groups0_val','groups1_val']])
text+='\nThe four nominal ≥3-per-class pairs involve PF08241, PF13489 and PF13649: related methyltransferase models in clan CL0063 (NADP_Rossmann). Their overlapping memberships are correlated annotation hypotheses, not four independent replications. None survives the conservative competing-domain exclusion at the same support threshold. Clan information is used to diagnose ambiguity, never to merge these into a single substitute family.\n'
text+='\n'+table(main[['k','same_source_training_pairs','same_source_val_both_pairs','same_source_both_splits_threshold_pairs','same_source_eval_unique_proteins','same_source_eval_val_proteins','same_source_eval_unique_groups']])
text+='''
Record-level experimental reliability is **unchanged and unverified for every record**. TargetTrack/SG statuses can encode workflow progress rather than confirmed insolubility; PDB records are structure-selected positives; DNASU-like records are candidates for explicit-status recovery but those statuses have not been recovered. Same-centre comparisons reduce a provenance shortcut without certifying comparable assay conditions. `family_label_reliability.tsv` gives the source × label distribution and reliability category inside each family; `family_source_support.tsv` contains the per-class group counts for every source stratum. Pure-positive sources can contribute to A and pooled family diversity but cannot independently supply a two-class B contrast.

## Fixed quality sensitivities

Domain-quality policies retain all A proteins and change usable family memberships only. Cohort-quality policies are diagnostic subsets, not alternate primary cohorts selected for performance. A stricter filter can lose power or preferentially remove one label; counts must be read alongside its source/label coverage table.

'''
text+=table(s[['policy','k','n_cohort','annotated_proteins','mixed_families','pooled_qualifying','training_eligible_val_both','both_splits_threshold','same_source_both_splits_threshold_pairs']])
text+='''
Under the conservative no-clan-competition policy, the three families with ≥5 proteins per label in both splits are PF00076, PF00400 and PF04055. These remain cross-source candidates, not verified same-assay comparisons.

The 80% HMM-span filter assesses partial-domain evidence, not homology coverage. The clan-competition sensitivity excludes both different-family same-clan hits with overlapping aligned residues; it is conservative and can remove legitimate nested domains. It is not presented as official PfamScan post-processing. The ESM-visible-domain policy requires alignment end ≤1,022; it avoids assigning a family based solely on a domain outside the cached input. Tags remain intact in the primary sequences.

## Reproducibility and next experiment

Follow [PFAM_METHODS.md](PFAM_METHODS.md) for the pinned reference URLs, checksums, package versions and exact search thresholds. `search_run.json` verifies all 21,979 profiles completed against all 13,880 sequences. `profile_thresholds.tsv` preserves each versioned accession and its sequence/domain cutoffs; `domains.tsv` preserves coordinates, scores and E-values. `annotated_cohort.tsv` includes no-hit rows; `memberships.tsv` records individual Pfam memberships. `analysis_B_memberships.tsv` links training-eligible strata back to A and flags validation and same-source evaluability. `training_selected_families.tsv` records support in both splits without changing split roles.

Next, fit property/concept models on the heterogeneous A training set and evaluate property predictions within the predefined B strata. Report all per-family/source counts, family-macro results, and uncertainty using detected homology groups as resampling units; account for proteins belonging to several Pfam strata. Use the ≥3 and ≥5 requirements as declared sensitivity levels, not a search over model outcomes. No representation probe was run for this annotation task. Solubility conclusions remain exploratory until original experimental outcomes are verified.

The old pilot-based report correctly described the evidence then available; this uniform run supersedes its zero-support conclusion. Original datasets, old analyses, family/split assignments and the locked test remain preserved.
'''
(P/'UNIFORM_PFAM_REPORT.md').write_text(text)
print('Wrote',P/'UNIFORM_PFAM_REPORT.md')
