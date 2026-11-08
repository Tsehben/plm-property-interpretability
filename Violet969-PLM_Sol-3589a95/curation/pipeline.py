#!/usr/bin/env python3
"""Train/development-only UESolDS audit. Never opens the independent test file.
Run prepare, homology, finish in order; see README.md. All outputs are additive.
"""
import argparse, collections, hashlib, importlib.util, json, os, platform, re, subprocess, sys, tempfile
from pathlib import Path
os.environ.setdefault('MPLCONFIGDIR', '/tmp/uesolds-matplotlib')
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from sklearn.preprocessing import OneHotEncoder, StandardScaler
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score, average_precision_score
import sklearn

BASE = Path(__file__).resolve().parents[1]
AA = 'ACDEFGHIKLMNPQRSTVWY'
KD = dict(zip('ACDEFGHIKLMNPQRSTVWY', [1.8,2.5,-3.5,-3.5,2.8,-0.4,-3.2,4.5,-3.9,3.8,1.9,-3.5,-1.6,-3.5,-4.5,-0.8,-0.7,4.2,-0.9,-1.3]))
# Reuse source parsing without running the legacy baseline main or changing its outputs.
oldargv = sys.argv[:]
sys.argv = ['baseline', str(BASE), 'unused']
spec = importlib.util.spec_from_file_location('legacy', BASE/'baseline_source_vs_sequence/baseline.py')
legacy = importlib.util.module_from_spec(spec)
spec.loader.exec_module(legacy)
sys.argv = oldargv

def sha(s): return hashlib.sha256(s.encode()).hexdigest()
def file_sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def save(df, p): df.to_csv(p, sep='\t', index=False, float_format='%.8g')
def read(p): return pd.read_csv(p, sep='\t', keep_default_na=False)
def counts(df, keys):
    return df.groupby(keys, dropna=False).agg(n=('label','size'), soluble=('label','sum'), unique_sequences=('sha256','nunique')).reset_index().assign(insoluble=lambda d:d.n-d.soluble, prevalence=lambda d:d.soluble/d.n)

def prepare(out):
    rows=[]; inputs={}
    for split, name in [('train','Train_dataset.fasta'), ('validation','validation_dataset.fasta')]:
        p=BASE/'embedding_dataset'/name; inputs[str(p.relative_to(BASE))]=file_sha(p)
        eligible_index=0
        for i,(header,s) in enumerate(legacy.parse_fasta(p)):
            pid=header.split()[0]; y=int(header.split()[-1].split('-')[-1]); assert y in (0,1)
            excluded=split=='train' and pid in legacy.EXCLUDE_FROM_TRAIN
            cls=legacy.id_class(pid); center=legacy.source_center(pid,cls)
            n=len(s); c=collections.Counter(s)
            rows.append(dict(raw_uid=f'{split}:{i}', raw_index=i, legacy_uid='' if excluded else f'{split}:{eligible_index}', original_id=pid, header=header, split=split,label=y,seq=s,sha256=sha(s),id_class=cls,center=center,length=n,
                prior_holdout_exclusion=excluded, nonstandard=''.join(sorted(set(s)-set(AA))), invalid=not s or bool(set(s)-set(AA+'BXZJUO')), ambiguous=bool(set(s)&set('BXZJ')), esm_truncated=n>1022,
                terminal_his=bool(re.search('H{6,}',s[:30]) or re.search('H{6,}',s[-30:])),
                terminal_cleavage=any(m in s[:40] or m in s[-40:] for m in ['ENLYFQ','LEVLFQ','LVPRGS']),
                terminal_linker=bool(re.search(r'[GS]{8,}',s[:40]) or re.search(r'[GS]{8,}',s[-40:])),
                non_m_start=not s.startswith('M'), charge_proxy=(c['K']+c['R']-c['D']-c['E'])/max(n,1),
                hydropathy=sum(KD.get(a,0)*v for a,v in c.items())/max(n,1), entropy=-sum(v/max(n,1)*np.log2(v/max(n,1)) for v in c.values()),
                **{'frac_'+a:c[a]/max(n,1) for a in AA}))
            if not excluded: eligible_index+=1
    df=pd.DataFrame(rows)
    df['conflicting_sequence']=df.groupby('sha256').label.transform('nunique')>1
    df['duplicate_sequence']=df.duplicated('sha256',False)
    df['repeated_identifier']=df.duplicated('original_id',False)
    # Cache keys refer to post-exclusion positions, not original FASTA positions.
    man=read(BASE/'dev_cohort/dev_cohort_manifest.tsv')
    mapping=df[df.legacy_uid!=''].merge(man, left_on='legacy_uid',right_on='record_uid',suffixes=('','_cache'),validate='one_to_one')
    assert len(mapping)==14000
    for col in ['sha256','original_id','label']:
        assert (mapping[col].astype(str)==mapping[col+'_cache'].astype(str)).all(), col
    seqpairs=read(BASE/'esm2_probe/uesolds_pairs.tsv')
    # The existing pairs TSV may name its sequence column differently.
    pairmap=dict(seqpairs.iloc[:,:2].itertuples(index=False,name=None))
    assert all(pairmap[r.legacy_uid]==r.seq for r in mapping.itertuples())
    df['cached']=df.legacy_uid.isin(man.record_uid)
    save(mapping[['raw_uid','raw_index','legacy_uid','original_id','sha256']],out/'cache_mapping.tsv')
    # Pilot annotations are attached by ID AND exact sequence, never by ID alone.
    sample=json.loads((BASE/'pilot_annotation/sample_100.json').read_text())
    sample_hash={r['id']:sha(r['seq']) for r in sample}
    ann=read(BASE/'pilot_annotation/uesolds_pilot_annotations.tsv')
    ann['sha256']=ann.original_id.map(sample_hash)
    ann=ann[['original_id','sha256','organism_status','source_organism','pfam_accessions']]
    df=df.merge(ann,on=['original_id','sha256'],how='left',validate='many_to_one')
    df['annotation_attempted']=df.organism_status.notna()
    df['organism_assigned']=df.source_organism.fillna('').ne('') & df.organism_status.fillna('').eq('confident')
    df['pfam_assigned']=df.pfam_accessions.fillna('').ne('')
    save(df.drop(columns='seq'),out/'records.tsv')
    save(df[df.duplicate_sequence|df.repeated_identifier].drop(columns='seq'),out/'redundancy.tsv')
    save(counts(df,['split','id_class','center','label']),out/'source_counts.tsv')
    features=['length','charge_proxy','hydropathy','entropy']+['frac_'+a for a in AA]
    save(df.groupby(['split','center','label'])[features].agg(['mean','median','std']).reset_index(),out/'descriptor_summary.tsv')
    flags=['terminal_his','terminal_cleavage','terminal_linker','non_m_start','esm_truncated','ambiguous','invalid']
    save(df.groupby(['split','center','label'])[flags].agg(['sum','mean']).reset_index(),out/'artifact_summary.tsv')
    save(df.groupby(['split','center','label'])[['annotation_attempted','organism_assigned','pfam_assigned']].agg(['sum','mean']).reset_index(),out/'annotation_coverage.tsv')
    save(df[df.annotation_attempted].groupby(['center','label'])[['organism_assigned','pfam_assigned']].agg(['size','sum','mean']).reset_index(),out/'annotation_attempted_only.tsv')
    # Direction and magnitude of differences, including sources with both labels.
    effects=[]
    for split,d in df.groupby('split'):
        for center,g in [('POOLED',d),*list(d.groupby('center'))]:
            if g.label.nunique()<2: continue
            for f in features+flags:
                a=g[g.label==1][f].astype(float); b=g[g.label==0][f].astype(float)
                sd=np.sqrt((a.var()+b.var())/2)
                effects.append(dict(split=split,center=center,feature=f,n0=len(b),n1=len(a),mean0=b.mean(),mean1=a.mean(),difference=a.mean()-b.mean(),standardized_difference=(a.mean()-b.mean())/sd if sd>0 else np.nan))
    save(pd.DataFrame(effects),out/'within_source_effects.tsv')
    # Choose sources using training only, before fitting or examining new scores.
    tr=df[(df.split=='train')&~df.prior_holdout_exclusion&~df.invalid&~df.ambiguous&~df.conflicting_sequence]
    ct=pd.crosstab(tr.center,tr.label).reindex(columns=[0,1],fill_value=0)
    centers=sorted(ct.index[(ct[0]>=100)&(ct[1]>=100)&ct.index.isin(tr[tr.id_class=='SG_target'].center)])
    config=dict(min_training_per_class=100,centers=centers,identity=0.30,coverage=0.80,evalue=0.001,sensitivity=7.5,seed=42,C=1.0,strict_label_status='unavailable: no record-level experimental outcomes in release')
    (out/'policy.json').write_text(json.dumps(config,indent=2)+'\n')
    # All valid train/dev sequences for MMseqs2, no label-dependent sequence filtering.
    for split in ['train','validation']:
        with (out/f'{split}.fasta').open('w') as f:
            for r in df[(df.split==split)&~df.invalid].itertuples(): f.write(f'>{r.raw_uid}\n{r.seq}\n')
    with (out/'all_development.fasta').open('w') as f:
        for r in df[~df.invalid].itertuples(): f.write(f'>{r.raw_uid}\n{r.seq}\n')
    for p in [BASE/'dev_cohort/dev_cohort_manifest.tsv',BASE/'esm2_probe/uesolds_pairs.tsv',BASE/'pilot_annotation/uesolds_pilot_annotations.tsv',BASE/'pilot_annotation/sample_100.json',*sorted((BASE/'esm2_probe/emb').glob('*'))]:
        if p.is_file(): inputs[str(p.relative_to(BASE))]=file_sha(p)
    (out/'input_hashes.json').write_text(json.dumps(inputs,indent=2)+'\n')
    (out/'environment.json').write_text(json.dumps(dict(python=sys.version,platform=platform.platform(),numpy=np.__version__,pandas=pd.__version__,sklearn=sklearn.__version__,matplotlib=matplotlib.__version__),indent=2)+'\n')
    # Source prevalence plot and conditional biochemical effect heatmap.
    fig,axs=plt.subplots(1,2,figsize=(13,6))
    for ax,(split,d) in zip(axs,df.groupby('split')):
        tab=pd.crosstab(d.center,d.label).reindex(columns=[0,1],fill_value=0)
        tab.plot.barh(stacked=True,ax=ax,color=['#d98268','#4e8ca8']); ax.set_title(split); ax.set_xlabel('Records'); ax.legend(['label 0','label 1'])
    fig.tight_layout();fig.savefig(out/'source_distributions.png',dpi=160);plt.close(fig)
    e=pd.DataFrame(effects); cols=['length','charge_proxy','hydropathy','entropy','terminal_his','non_m_start']
    mat=e[(e.split=='train')&e.feature.isin(cols)].pivot(index='center',columns='feature',values='standardized_difference').reindex(columns=cols)
    fig,ax=plt.subplots(figsize=(9,6)); im=ax.imshow(mat,aspect='auto',cmap='RdBu_r',vmin=-1,vmax=1);ax.set_yticks(range(len(mat)),mat.index);ax.set_xticks(range(len(cols)),cols,rotation=35,ha='right');fig.colorbar(im,ax=ax,label='Standardized difference: label 1 minus label 0');fig.tight_layout();fig.savefig(out/'conditional_effects.png',dpi=160);plt.close(fig)
    print('Prepared',len(df),'records. Training-selected centers:',centers,flush=True)

def homology(out, mmseqs):
    commands=[]
    with tempfile.TemporaryDirectory(prefix='uesolds-mmseqs-') as tmp:
        # Direct cross-split search avoids relying on cluster representatives for leakage.
        cmd=[mmseqs,'easy-search',str(out/'validation.fasta'),str(out/'train.fasta'),str(out/'cross_split_hits.tsv'),tmp+'/search','--min-seq-id','0.30','-c','0.80','--cov-mode','0','-s','7.5','-e','0.001','--max-seqs','10000','--alignment-mode','3','--threads','4','--format-output','query,target,fident,qcov,tcov,evalue,bits']
        commands.append(cmd)
        with (out/'mmseqs_search.log').open('w') as f: subprocess.run(cmd,stdout=f,stderr=subprocess.STDOUT,check=True)
        cmd=[mmseqs,'easy-cluster',str(out/'all_development.fasta'),str(out/'diversity'),tmp+'/cluster','--min-seq-id','0.30','-c','0.80','--cov-mode','0','-s','7.5','-e','0.001','--alignment-mode','3','--cluster-mode','1','--linclust-version','1','--threads','4']
        commands.append(cmd)
        with (out/'mmseqs_cluster.log').open('w') as f: subprocess.run(cmd,stdout=f,stderr=subprocess.STDOUT,check=True)
    (out/'homology_run.json').write_text(json.dumps(dict(version=subprocess.check_output([mmseqs,'version'],text=True).strip(),binary_sha256=file_sha(mmseqs),commands=commands,input_hashes={p.name:file_sha(p) for p in [out/'train.fasta',out/'validation.fasta',out/'all_development.fasta']}),indent=2)+'\n')

def finish(out):
    df=read(out/'records.tsv'); policy=json.loads((out/'policy.json').read_text())
    run=json.loads((out/'homology_run.json').read_text())
    assert all(file_sha(out/p)==v for p,v in run['input_hashes'].items())
    hits=pd.read_csv(out/'cross_split_hits.tsv',sep='\t',names=['query','target','identity','qcov','tcov','evalue','bits'])
    assert (hits.identity>=0.30).all() and (hits.qcov>=0.80).all() and (hits.tcov>=0.80).all()
    assert hits['query'].isin(df.loc[df.split=='validation','raw_uid']).all()
    assert hits.target.isin(df.loc[df.split=='train','raw_uid']).all()
    clusters=pd.read_csv(out/'diversity_cluster.tsv',sep='\t',names=['cluster','raw_uid'])
    assert clusters.raw_uid.is_unique
    assert set(clusters.raw_uid)==set(df.loc[~df.invalid,'raw_uid'])
    df=df.merge(clusters,on='raw_uid',how='left',validate='one_to_one')
    df['homology_to_validation']=df.raw_uid.isin(hits.target)
    # Conservatively also remove training members of any cluster spanning validation.
    valclusters=set(df.loc[df.split=='validation','cluster'].dropna())
    df['cluster_spans_validation']=df.cluster.isin(valclusters)&df.split.eq('train')
    validation_ids=set(df.loc[df.split=='validation','original_id'])
    df['identifier_in_validation']=df.split.eq('train')&df.original_id.isin(validation_ids)
    valhash=set(df.loc[df.split=='validation','sha256'])
    df['exact_in_validation']=df.split.eq('train')&df.sha256.isin(valhash)
    df['duplicate_within_split']=df.duplicated(['split','sha256'])
    reasons=[]
    for r in df.itertuples():
        reasons.append(';'.join(k for k in ['prior_holdout_exclusion','invalid','ambiguous','conflicting_sequence','duplicate_within_split','exact_in_validation','homology_to_validation','cluster_spans_validation','identifier_in_validation'] if getattr(r,k)))
    df['exclusion_reasons']=reasons
    df['clean']=df.exclusion_reasons.eq('')
    df['source_controlled']=df.clean&df.center.isin(policy['centers'])&df.id_class.eq('SG_target')
    df['strict_label']=False
    df['strict_label_reason']='no construct-matched experimental outcome recovered'
    manifest_cols=['raw_uid','legacy_uid','original_id','split','label','sha256','center','id_class','length','cached','cluster','exclusion_reasons','clean','source_controlled','strict_label','strict_label_reason']
    save(df[manifest_cols],out/'cohort_manifest.tsv')
    cohorts={'original':~df.prior_holdout_exclusion,'clean_all_sources':df.clean,'source_controlled':df.source_controlled,'strict_label':df.strict_label}
    stats=[]
    for name,mask in cohorts.items():
        save(df.loc[mask,manifest_cols],out/f'{name}_manifest.tsv')
        for scope,scope_mask in [('full',np.ones(len(df),bool)),('cached',df.cached)]:
            for split in ['train','validation']:
                g=df[mask&scope_mask&df.split.eq(split)]
                stats.append(dict(cohort=name,scope=scope,split=split,n=len(g),soluble=int(g.label.sum()),insoluble=int(len(g)-g.label.sum()),unique_sequences=g.sha256.nunique(),clusters=g.cluster.nunique(),centers=g.center.nunique(),annotation_attempted=int(g.annotation_attempted.sum()),organism_assigned=int(g.organism_assigned.sum()),pfam_assigned=int(g.pfam_assigned.sum())))
        save(counts(df[mask],['split','center']),out/f'{name}_source_counts.tsv')
    save(pd.DataFrame(stats),out/'cohort_statistics.tsv')
    save(hits.merge(df[['raw_uid','center','label']],left_on='query',right_on='raw_uid').drop(columns='raw_uid').merge(df[['raw_uid','center','label']],left_on='target',right_on='raw_uid',suffixes=('_validation','_train')).drop(columns='raw_uid'),out/'homology_pairs_annotated.tsv')
    # Cached-only model comparison uses precisely the same rows for every feature set.
    vectors={}
    for p in sorted((BASE/'esm2_probe/emb').glob('uesolds_shard*.npz')):
        with np.load(p,allow_pickle=False) as z:
            for k,v in zip(z['keys'],z['whole']):
                assert str(k) not in vectors
                vectors[str(k)]=v
    pairs=read(BASE/'esm2_probe/uesolds_pairs.tsv');seqmap=dict(pairs.iloc[:,:2].itertuples(index=False,name=None))
    results=[];predictions=[]
    for name,mask in cohorts.items():
        d=df[mask&df.cached]; tr=d[d.split=='train'];va=d[d.split=='validation']
        if tr.label.nunique()<2 or va.empty: continue
        y=tr.label.to_numpy(); yv=va.label.to_numpy()
        def comp(g):
            ss=[seqmap[k][:1022] for k in g.legacy_uid]
            return np.array([[s.count(a)/len(s) for a in AA]+[np.log(len(s))] for s in ss])
        enc=OneHotEncoder(handle_unknown='ignore',sparse_output=False).fit(tr[['id_class','center']])
        features={'source':(enc.transform(tr[['id_class','center']]),enc.transform(va[['id_class','center']])), 'composition':(comp(tr),comp(va)), 'ESM_layer8':(np.stack([vectors[k] for k in tr.legacy_uid]),np.stack([vectors[k] for k in va.legacy_uid]))}
        for model,(x,xv) in features.items():
            if model!='source':
                scaler=StandardScaler().fit(x); x=scaler.transform(x);xv=scaler.transform(xv)
            assert np.isfinite(x).all() and np.isfinite(xv).all()
            fit=LogisticRegression(C=1,solver='lbfgs',max_iter=5000,random_state=42).fit(x,y)
            assert fit.n_iter_[0]<5000
            p=fit.predict_proba(xv)[:,1]
            for uid,yy,pp in zip(va.raw_uid,yv,p): predictions.append(dict(cohort=name,model=model,raw_uid=uid,label=yy,prediction=pp))
            groups={'POOLED':np.ones(len(va),bool),**{c:va.center.eq(c).to_numpy() for c in sorted(va.center.unique())}}
            center_scores=[]
            for group,idx in groups.items():
                yy=yv[idx]; defined=len(set(yy))==2
                auc=roc_auc_score(yy,p[idx]) if defined else np.nan
                results.append(dict(cohort=name,model=model,group=group,n_train=len(tr),n_eval=int(idx.sum()),n0=int((yy==0).sum()),n1=int(yy.sum()),auroc=auc,average_precision=average_precision_score(yy,p[idx]) if defined else np.nan))
                if group in policy['centers'] and defined: center_scores.append(auc)
            results.append(dict(cohort=name,model=model,group='MACRO_training_selected_centers',n_train=len(tr),n_eval=len(center_scores),auroc=np.mean(center_scores) if center_scores else np.nan))
            print(name,model,'pooled',results[-len(groups)-1]['auroc'],'macro',results[-1]['auroc'],flush=True)
    save(pd.DataFrame(results),out/'baseline_results.tsv');save(pd.DataFrame(predictions),out/'predictions.tsv')
    clean=df[df.clean];tr=clean[clean.split=='train'];va=clean[clean.split=='validation']
    assert not set(tr.sha256)&set(va.sha256)
    assert not set(tr.cluster.dropna())&set(va.cluster.dropna())
    assert not set(tr.raw_uid)&set(hits[hits['query'].isin(va.raw_uid)].target)
    assert not set(tr.original_id)&set(va.original_id)
    assert len(vectors)==14000 and df.cached.sum()==14000
    # Verify input integrity after the entire workflow, including cached artifacts.
    assert all(file_sha(BASE/p)==h for p,h in json.loads((out/'input_hashes.json').read_text()).items())
    summary=dict(records=len(df),cache_index_displacements=int((read(out/'cache_mapping.tsv').raw_uid!=read(out/'cache_mapping.tsv').legacy_uid).sum()),exact_duplicate_records=int(df.duplicate_sequence.sum()),conflicting_records=int(df.conflicting_sequence.sum()),repeated_identifier_records=int(df.repeated_identifier.sum()),cross_split_pairs=len(hits),train_with_hits=hits.target.nunique(),validation_with_hits=hits['query'].nunique(),clusters=df.cluster.nunique(),selected_centers=policy['centers'],strict_label_n=0,checks='passed: cache mapping, finite models, convergence, exact/cluster/hit/identifier separation, input hashes; independent test never opened')
    (out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps(summary,indent=2))

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('stage',choices=['prepare','homology','finish']);ap.add_argument('--out',type=Path,default=Path(__file__).resolve().parent/'results');ap.add_argument('--mmseqs',default='mmseqs');a=ap.parse_args();a.out.mkdir(parents=True,exist_ok=True)
    if a.stage=='prepare': prepare(a.out)
    elif a.stage=='homology':homology(a.out,a.mmseqs)
    else:finish(a.out)
