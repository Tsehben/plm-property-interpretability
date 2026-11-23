#!/usr/bin/env python3
"""Uniform Pfam coverage and within-family feasibility; no model performance input."""
import argparse, hashlib, json, os
from pathlib import Path
import pandas as pd
import numpy as np
from family_feasibility import family_counts,qualify
os.environ.setdefault('MPLCONFIGDIR','/tmp/uesolds-matplotlib')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parent

def read(p):return pd.read_csv(p,sep='\t',keep_default_na=False)
def save(d,p):d.to_csv(p,sep='\t',index=False,float_format='%.9g')
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def run(out):
    runmeta=json.loads((out/'search_run.json').read_text());assert runmeta['completed']
    assert sha(out/'domains.tsv')==runmeta['domains_sha256']
    d=read(out/'searched_manifest.tsv');assert len(d)==13880
    rec=read(ROOT/'results/records.tsv')
    d=d.merge(rec[['raw_uid','terminal_his','terminal_cleavage','terminal_linker','esm_truncated','organism_assigned']],on='raw_uid',validate='one_to_one')
    dom=read(out/'domains.tsv');assert dom.raw_uid.isin(d.raw_uid).all()
    assert ((dom.sequence_score>=dom.sequence_ga)&(dom.domain_score>=dom.domain_ga)).all()
    assert ((dom.ali_from>=1)&(dom.ali_to<=dom.sequence_length)&(dom.ali_to>=dom.ali_from)).all()
    clans=pd.read_csv(out/'Pfam-A.clans.tsv',sep='\t',names=['pfam','clan','clan_name','family_name','description'],keep_default_na=False)
    dom=dom.merge(clans[['pfam','clan']],on='pfam',how='left',validate='many_to_one').fillna({'clan':''})
    assert dom.pfam.isin(clans.pfam).all()
    dom['clan_competing_overlap']=False;overlaps=[]
    # Flag BOTH alternatives rather than choosing one using labels/source.
    for (uid,clan),g in dom[dom.clan.ne('')].groupby(['raw_uid','clan']):
        items=list(g.itertuples())
        for i,a in enumerate(items):
            for b in items[i+1:]:
                overlap=min(a.ali_to,b.ali_to)-max(a.ali_from,b.ali_from)+1
                if a.pfam!=b.pfam and overlap>0:
                    dom.loc[[a.Index,b.Index],'clan_competing_overlap']=True
                    overlaps.append(dict(raw_uid=uid,clan=clan,pfam_a=a.pfam,pfam_b=b.pfam,overlap_residues=overlap))
    save(pd.DataFrame(overlaps,columns=['raw_uid','clan','pfam_a','pfam_b','overlap_residues']),out/'clan_competition.tsv')
    save(dom,out/'domains_flagged.tsv')
    d['label_reliability']=d.id_class.map({'SG_target':'workflow-derived; construct outcome unverified','DNASU_clone':'explicit-outcome candidate; status not recovered','PDB_serial':'structure-selected positive; standardized assay absent'}).fillna('original source/measurement unresolved')
    membership=dom[['raw_uid','pfam','pfam_versioned','name','clan']].drop_duplicates(['raw_uid','pfam']).merge(d,on='raw_uid',validate='many_to_one')
    save(membership,out/'memberships.tsv')
    byid=membership.groupby('raw_uid').pfam.agg(lambda a:';'.join(sorted(set(a))))
    d['pfam_accessions']=d.raw_uid.map(byid).fillna('');d['pfam_status']=np.where(d.pfam_accessions.ne(''),'assigned','searched_no_GA_hit')
    d['annotation_attempted']=True;d['pfam_release']='37.0';d['verified_solubility_outcome']=False
    save(d,out/'annotated_cohort.tsv')
    policies={
        'GA_all':(np.ones(len(d),bool),np.ones(len(dom),bool)),
        'GA_hmm_coverage_80pct':(np.ones(len(d),bool),dom.hmm_coverage>=0.8),
        'GA_without_clan_competition':(np.ones(len(d),bool),~dom.clan_competing_overlap),
        'GA_domains_visible_to_ESM':(np.ones(len(d),bool),dom.ali_to<=1022),
        'source_controlled':(d.source_controlled,np.ones(len(dom),bool)),
        'no_terminal_artifact_flags':(~(d.terminal_his|d.terminal_cleavage|d.terminal_linker),np.ones(len(dom),bool)),
        'within_ESM_context':(~d.esm_truncated,np.ones(len(dom),bool)),
    }
    coverage=[];summary=[];support=[];ssupport=[];candidates=[];eligible=[]
    for policy,(rm,dm) in policies.items():
        a=d[rm];m=dom[dm][['raw_uid','pfam']].drop_duplicates().merge(a,on='raw_uid',validate='many_to_one')
        a=a.copy();a['assigned']=a.raw_uid.isin(m.raw_uid)
        for (split,center,label),g in a.groupby(['split','center','label']):
            coverage.append(dict(policy=policy,split=split,center=center,label=label,n=len(g),attempted=len(g),assigned=int(g.assigned.sum()),no_usable_hit=int((~g.assigned).sum()),coverage=g.assigned.mean()))
        tabs={}
        for scope in ['development','train','validation']:
            mm=m if scope=='development' else m[m.split==scope]
            t=family_counts(mm,['pfam']);st=family_counts(mm,['pfam','center']);tabs[scope]=(t,st)
            support.append(t.assign(policy=policy,split_scope=scope));ssupport.append(st.assign(policy=policy,split_scope=scope))
        for k in [3,5]:
            t,st=tabs['development'];tr,str_=tabs['train'];va,sva=tabs['validation']
            combined=qualify(t,k);train=qualify(tr,k);valid=qualify(va,k)
            train_ids=set(train.pfam);val_both=set(va[(va.n0>0)&(va.n1>0)].pfam)
            eval_ids=train_ids&val_both;strict_ids=train_ids&set(valid.pfam)
            sg=qualify(str_,k);sgkeys=set(zip(sg.pfam,sg.center))
            sv=sva[(sva.n0>0)&(sva.n1>0)];evalsg=sgkeys&set(zip(sv.pfam,sv.center))
            strictsg=sgkeys&set(zip(qualify(sva,k).pfam,qualify(sva,k).center))
            em=m[m.pfam.isin(combined.pfam)]
            evalm=m[m.pfam.isin(eval_ids)]
            sg_m=m.merge(pd.DataFrame(sorted(evalsg),columns=['pfam','center']),on=['pfam','center'])
            iv=qualify(tr,k);iv=iv[(iv.groups_exclusive0>=k)&(iv.groups_exclusive1>=k)]
            ivv=va[(va.groups_exclusive0>=k)&(va.groups_exclusive1>=k)]
            summary.append(dict(policy=policy,k=k,n_cohort=len(a),annotated_proteins=m.raw_uid.nunique(),families=len(t),single_class_families=int(((t.n0==0)|(t.n1==0)).sum()),mixed_families=int(((t.n0>0)&(t.n1>0)).sum()),pooled_qualifying=len(combined),pooled_eligible_proteins=em.raw_uid.nunique(),pooled_eligible_groups=em.cluster.nunique(),pooled_eligible_groups0=em[em.label==0].cluster.nunique(),pooled_eligible_groups1=em[em.label==1].cluster.nunique(),training_eligible=len(train_ids),training_eligible_val_both=len(eval_ids),both_splits_threshold=len(strict_ids),both_splits_k_exclusive_groups=len(set(iv.pfam)&set(ivv.pfam)),eval_eligible_unique_proteins=evalm.raw_uid.nunique(),eval_eligible_unique_groups=evalm.cluster.nunique(),eval_eligible_groups0=evalm[evalm.label==0].cluster.nunique(),eval_eligible_groups1=evalm[evalm.label==1].cluster.nunique(),eval_eligible_val_proteins=evalm[evalm.split=='validation'].raw_uid.nunique(),same_source_training_pairs=len(sgkeys),same_source_val_both_pairs=len(evalsg),same_source_both_splits_threshold_pairs=len(strictsg),same_source_val_both_families=len({f for f,c in evalsg}),same_source_eval_unique_proteins=sg_m.raw_uid.nunique(),same_source_eval_unique_groups=sg_m.cluster.nunique(),same_source_eval_val_proteins=sg_m[sg_m.split=='validation'].raw_uid.nunique()))
            if policy=='GA_all':
                cand=tr.merge(va,on='pfam',how='left',suffixes=('_train','_validation')).fillna(0)
                cand=cand[cand.pfam.isin(train_ids)].copy();cand['threshold']=k;candidates.append(cand)
                b=m[m.pfam.isin(train_ids)].copy();b['threshold']=k
                b['validation_both_classes']=b.pfam.isin(eval_ids)
                b['same_source_training_eligible']=[(f,c) in sgkeys for f,c in zip(b.pfam,b.center)]
                b['same_source_validation_both']=[(f,c) in evalsg for f,c in zip(b.pfam,b.center)]
                eligible.append(b[['raw_uid','pfam','split','label','center','cluster','threshold','validation_both_classes','same_source_training_eligible','same_source_validation_both','label_reliability']])
    save(pd.DataFrame(coverage),out/'coverage_by_source_label.tsv')
    save(pd.DataFrame(summary),out/'feasibility_summary.tsv')
    save(pd.concat(support),out/'family_support.tsv');save(pd.concat(ssupport),out/'family_source_support.tsv')
    save(pd.concat(candidates),out/'training_selected_families.tsv');save(pd.concat(eligible),out/'analysis_B_memberships.tsv')
    b_all=pd.concat(eligible)
    source_rows=[]
    for k,g in b_all.groupby('threshold'):
        for kind,mask in [('training_eligible',np.ones(len(g),bool)),('validation_both',g.validation_both_classes),('same_source_validation_both',g.same_source_validation_both)]:
            for (split,center,label),v in g[mask].groupby(['split','center','label']):
                source_rows.append(dict(threshold=k,contrast=kind,split=split,center=center,label=label,n_proteins=v.raw_uid.nunique(),n_groups=v.cluster.nunique(),verified_outcomes=0))
    save(pd.DataFrame(source_rows),out/'B_source_distribution.tsv')
    save(membership.groupby(['pfam','center','split','label','label_reliability']).size().rename('n').reset_index(),out/'family_label_reliability.tsv')
    # B is always a view of A. Preserve all 13,880 A rows and their roles.
    old=read(ROOT/'results/cohort_manifest.tsv');old=old[old.clean&old.cached]
    assert set(d.raw_uid)==set(old.raw_uid)
    for col in ['split','label','sha256','legacy_uid','cluster']:
        assert d.set_index('raw_uid')[col].sort_index().equals(old.set_index('raw_uid')[col].sort_index())
    assert set(pd.concat(eligible).raw_uid)<=set(d.raw_uid)
    assert not set(d[d.split=='train'].cluster)&set(d[d.split=='validation'].cluster)
    cov=pd.DataFrame(coverage);v=cov[cov.policy=='GA_all']
    plot=v.groupby(['center','label'])[['n','assigned']].sum();plot['percent']=100*plot.assigned/plot.n
    fig,ax=plt.subplots(figsize=(11,5));plot.percent.unstack().plot.bar(ax=ax,color=['#d98268','#4e8ca8']);ax.set_ylim(0,100);ax.set_ylabel('Proteins assigned at least one Pfam (%)');ax.set_xlabel('Source / centre');ax.legend(['label 0','label 1']);fig.tight_layout();fig.savefig(out/'coverage.png',dpi=160);plt.close(fig)
    (out/'analysis_checks.json').write_text(json.dumps(dict(n_A=len(d),all_attempted=True,split_labels_hashes_preserved=True,clusters_disjoint_across_splits=True,B_subset_A=True,domain_thresholds_and_coordinates_valid=True,analysis_script_sha256=sha(Path(__file__))),indent=2)+'\n')
    print(pd.DataFrame(summary).to_string(index=False))
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--out',type=Path,default=ROOT/'pfam_annotation');a=ap.parse_args();run(a.out)
