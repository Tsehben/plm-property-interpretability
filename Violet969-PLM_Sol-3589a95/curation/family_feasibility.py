#!/usr/bin/env python3
"""Audit A/B feasibility from the same curated cohort; no models or test access.
Pfam memberships come only from the existing pilot's evidence, matched by ID+hash.
Unannotated proteins remain in analysis A. See FAMILY_FEASIBILITY.md.
"""
import argparse
import hashlib
import json
import os
import re
from pathlib import Path
import numpy as np
import pandas as pd
os.environ.setdefault('MPLCONFIGDIR', '/tmp/uesolds-matplotlib')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

BASE=Path(__file__).resolve().parents[1]
ROOT=Path(__file__).resolve().parent

def read(p): return pd.read_csv(p,sep='\t',keep_default_na=False)
def save(d,p): d.to_csv(p,sep='\t',index=False,float_format='%.8g')
def digest(p): return hashlib.sha256(p.read_bytes()).hexdigest()

def family_counts(m,by):
    """Protein membership counts and distinct detected clusters, never cluster=family."""
    cols=by+['n0','n1','n_proteins','groups0','groups1','groups_total','groups_shared_labels','groups_exclusive0','groups_exclusive1','n_sources','n_train','n_validation','verified_outcomes','source_label_counts']
    rows=[]
    for key,g in m.groupby(by,sort=True):
        key=key if isinstance(key,tuple) else (key,)
        g=g.drop_duplicates('raw_uid')
        zero=g[g.label==0]; one=g[g.label==1]
        c0=set(zero.cluster);c1=set(one.cluster)
        rows.append(dict(zip(by,key),n0=len(zero),n1=len(one),n_proteins=len(g),groups0=len(c0),groups1=len(c1),groups_total=len(c0|c1),groups_shared_labels=len(c0&c1),groups_exclusive0=len(c0-c1),groups_exclusive1=len(c1-c0),n_sources=g.center.nunique(),n_train=int(g.split.eq('train').sum()),n_validation=int(g.split.eq('validation').sum()),verified_outcomes=0,source_label_counts=json.dumps({f'{c}|{y}':int(n) for (c,y),n in g.groupby(['center','label']).size().items()},sort_keys=True)))
    return pd.DataFrame(rows,columns=cols)

def qualify(t,k): return t[(t.n0>=k)&(t.n1>=k)]

def main(out):
    out.mkdir(parents=True,exist_ok=True)
    inputs=[ROOT/'results/records.tsv',ROOT/'results/cohort_manifest.tsv',BASE/'pilot_annotation/uesolds_pilot_annotations.tsv',BASE/'pilot_annotation/sample_100.json']
    hashes={str(p.relative_to(BASE)):digest(p) for p in inputs}
    d=read(inputs[0]); manifest=read(inputs[1]); ann=read(inputs[2])
    sample=json.loads(inputs[3].read_text())
    h={r['id']:hashlib.sha256(r['seq'].encode()).hexdigest() for r in sample}
    assert len(h)==len(sample)==100
    ann['sha256']=ann.original_id.map(h)
    assert ann.sha256.notna().all()
    assert not ann.duplicated(['original_id','sha256']).any()
    # Preserve all curated proteins even when evidence is absent.
    d=d.merge(manifest[['raw_uid','cluster','clean','source_controlled','strict_label']],on='raw_uid',validate='one_to_one')
    d=d.merge(ann[['original_id','sha256','match_method','evidence_source','raw_label']],on=['original_id','sha256'],how='left',validate='many_to_one')
    for col in ['match_method','evidence_source']:
        d[col]=d[col].fillna('')
    assert (d.loc[d.annotation_attempted,'raw_label'].astype(int)==d.loc[d.annotation_attempted,'label']).all()
    assert set(d.split)<=set(['train','validation'])
    assert d.raw_uid.is_unique
    d['pfam_status']=np.select([~d.annotation_attempted,d.pfam_assigned],['not_queried','assigned'],default='attempted_no_assignment')
    d['label_reliability']=d.id_class.map({'SG_target':'workflow-derived; record-level outcome unverified','PDB_serial':'structure-selected positive; no standardized assay','DNASU_clone':'explicit-outcome candidate; record-level status unverified'}).fillna('original source/measurement unresolved')
    d['verified_outcome']=False
    d['pfam_ids']=d.pfam_accessions.map(lambda s:sorted(set(s.split(';'))-set([''])))
    assert all(re.fullmatch(r'PF\d{5}',pf) for ids in d.pfam_ids for pf in ids)
    d['pfam_set']=d.pfam_ids.map(';'.join) # unordered membership set, NOT domain architecture
    primary=d[d.clean].copy()
    save(primary[['raw_uid','legacy_uid','split','label','center','id_class','cluster','cached','sha256','pfam_status','pfam_set','organism_status','source_organism','match_method','evidence_source','label_reliability','verified_outcome']],out/'primary_analysis_A_manifest.tsv')
    filters={
        'primary_clean_all_sources':d.clean,
        'source_controlled_sensitivity':d.source_controlled,
        'no_terminal_artifact_flags':d.clean&~(d.terminal_his|d.terminal_cleavage|d.terminal_linker),
        'within_ESM_context':d.clean&~d.esm_truncated,
        'confident_organism_only':d.clean&d.organism_assigned,
        'strict_verified_labels':d.clean&d.strict_label,
        # Evidence policy only: retains the SAME A proteins, filters family assignments.
        'full_sequence_Pfam_evidence':d.clean,
    }
    coverage=[]; support=[]; source_support=[]; summaries=[]; b_rows=[]; quality=[]
    for policy,mask in filters.items():
        for scope in ['full','cached']:
            g=d[mask&(d.cached if scope=='cached' else True)].copy()
            if policy=='full_sequence_Pfam_evidence':
                g.loc[g.match_method!='exact_full_sequence','pfam_ids']=g.loc[g.match_method!='exact_full_sequence','pfam_ids'].map(lambda _:[])
            g['usable_pfam']=g.pfam_ids.map(bool)
            for (split,center,label),v in g.groupby(['split','center','label']):
                coverage.append(dict(policy=policy,scope=scope,split=split,center=center,label=label,n=len(v),attempted=int(v.annotation_attempted.sum()),not_queried=int((~v.annotation_attempted).sum()),attempted_no_assignment=int((v.annotation_attempted&~v.pfam_assigned).sum()),assigned_original=int(v.pfam_assigned.sum()),usable_pfam=int(v.usable_pfam.sum()),coverage=v.usable_pfam.mean(),success_among_attempted=v.usable_pfam.sum()/v.annotation_attempted.sum() if v.annotation_attempted.sum() else np.nan))
            m=g.explode('pfam_ids').rename(columns={'pfam_ids':'pfam'})
            m=m[m.pfam.notna() & m.pfam.ne('')].drop_duplicates(['raw_uid','pfam'])
            for split in ['development','train','validation']:
                v=g if split=='development' else g[g.split==split]
                z=m if split=='development' else m[m.split==split]
                t=family_counts(z,['pfam']); ts=family_counts(z,['pfam','center'])
                for table,store in [(t,support),(ts,source_support)]:
                    if len(table): store.append(table.assign(policy=policy,scope=scope,split_scope=split))
                for k in [3,5]:
                    eligible=qualify(t,k); same=qualify(ts,k)
                    proteins=z[z.pfam.isin(eligible.pfam)]
                    independent=eligible[(eligible.groups_exclusive0>=k)&(eligible.groups_exclusive1>=k)]
                    same_ids=set(same.pfam)
                    same_proteins=z.merge(same[['pfam','center']],on=['pfam','center'],how='inner')
                    # An evaluated contrast needs threshold support in BOTH published splits.
                    ttr=family_counts(m[m.split=='train'],['pfam']); tva=family_counts(m[m.split=='validation'],['pfam'])
                    eval_ids=set(qualify(ttr,k).pfam)&set(qualify(tva,k).pfam)
                    summaries.append(dict(policy=policy,scope=scope,split_scope=split,threshold_per_class=k,n_cohort=len(v),n0=int(v.label.eq(0).sum()),n1=int(v.label.eq(1).sum()),n_groups=v.cluster.nunique(),annotated_proteins=int(v.usable_pfam.sum()),attempted=int(v.annotation_attempted.sum()),distinct_pfam=len(t),single_class_families=int(((t.n0==0)|(t.n1==0)).sum()),mixed_label_families=int(((t.n0>0)&(t.n1>0)).sum()),qualifying_families=len(eligible),eligible_unique_proteins=proteins.raw_uid.nunique(),eligible_unique_groups=proteins.cluster.nunique(),families_with_k_exclusive_groups_each=len(independent),same_source_qualifying_families=len(same_ids),same_source_family_pairs=len(same),same_source_unique_proteins=same_proteins.raw_uid.nunique(),families_meeting_threshold_in_both_splits=len(eval_ids)))
                    if policy=='primary_clean_all_sources' and split=='development':
                        b=proteins[['raw_uid','legacy_uid','pfam','split','label','center','cluster','cached','label_reliability']].copy()
                        b['scope']=scope;b['threshold_per_class']=k;b['same_source_qualifies']=[(p,c) in set(zip(same.pfam,same.center)) for p,c in zip(b.pfam,b.center)]
                        b_rows.append(b)
            quality.append(dict(policy=policy,scope=scope,n=len(g),n0=int(g.label.eq(0).sum()),n1=int(g.label.eq(1).sum()),n_centers=g.center.nunique(),confident_organisms=g.loc[g.organism_assigned,'source_organism'].nunique(),verified_outcomes=0))
    save(pd.DataFrame(coverage),out/'coverage_by_source_label.tsv')
    save(pd.concat(support,ignore_index=True),out/'family_support.tsv')
    save(pd.concat(source_support,ignore_index=True),out/'family_source_support.tsv')
    summary=pd.DataFrame(summaries);save(summary,out/'feasibility_summary.tsv');save(pd.DataFrame(quality),out/'quality_filter_counts.tsv')
    b=pd.concat(b_rows,ignore_index=True);save(b,out/'controlled_analysis_B_memberships.tsv')
    assert set(b.raw_uid)<=set(primary.raw_uid)
    memberships=primary.explode('pfam_ids').rename(columns={'pfam_ids':'pfam'})
    memberships=memberships[memberships.pfam.notna()&memberships.pfam.ne('')].drop_duplicates(['raw_uid','pfam'])
    save(memberships[['raw_uid','pfam','pfam_set','split','label','center','cluster','cached','match_method','evidence_source','label_reliability','verified_outcome']],out/'primary_pfam_memberships.tsv')
    save(primary.groupby(['center','split','label','label_reliability']).size().rename('n').reset_index(),out/'label_reliability_by_source.tsv')
    # Explicit zero-row schemas preserve the absence of verified qualifying families.
    fig,ax=plt.subplots(figsize=(9,4))
    cov=pd.DataFrame(coverage);cov=cov[(cov.policy=='primary_clean_all_sources')&(cov.scope=='full')&(cov.split=='train')]
    plotted=cov.pivot(index='center',columns='label',values='coverage')*100
    plotted.plot.bar(ax=ax,color=['#d98268','#4e8ca8']);ax.set_ylabel('Training proteins with usable Pfam (%)');ax.set_xlabel('Source proxy / centre');ax.legend(['label 0','label 1']);fig.tight_layout();fig.savefig(out/'pfam_coverage.png',dpi=160);plt.close(fig)
    # All original audit inputs stay intact and no test path is accepted by the CLI.
    assert all(digest(BASE/p)==h for p,h in hashes.items())
    (out/'provenance.json').write_text(json.dumps(dict(inputs=hashes,script_sha256=digest(Path(__file__)),pfam_method='existing pilot UniProt/RCSB cross-references; no new domain scan or cluster propagation',pfam_release='not pinned in historical pilot',thresholds=[3,5],primary='clean_all_sources',test_access=False,assertions='unique records; ID/hash/label joins; Pfam identifiers; B subset of A; input integrity'),indent=2)+'\n')
    print(summary[(summary.scope=='full')&(summary.split_scope=='development')].to_string(index=False))

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--out',type=Path,default=ROOT/'family_results');a=ap.parse_args();main(a.out)
