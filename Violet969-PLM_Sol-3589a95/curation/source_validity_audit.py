#!/usr/bin/env python3
"""Additive source audit. Existing development artifacts only; no test or model extraction.
Usage: python curation/source_validity_audit.py [--stage diversity|probes|all]
"""
import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','VECLIB_MAXIMUM_THREADS']:
    os.environ.setdefault(k,'2')
os.environ.setdefault('MPLCONFIGDIR','/tmp/uesolds-matplotlib')
import argparse,hashlib,importlib.util,json,platform,time,warnings
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import sparse
from sklearn.preprocessing import StandardScaler,MultiLabelBinarizer
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score,average_precision_score,log_loss,confusion_matrix
from sklearn.exceptions import ConvergenceWarning
from threadpoolctl import threadpool_limits
import sklearn
ROOT=Path(__file__).resolve().parent;BASE=ROOT.parent;OUT=ROOT/'source_audit_results'
AA='ACDEFGHIKLMNPQRSTVWY';BOOT=200

def read(p):return pd.read_csv(p,sep='\t',keep_default_na=False,float_precision='round_trip')
def save(d,name):d.to_csv(OUT/name,sep='\t',index=False,float_format='%.17g')
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def emit(s):print(s,flush=True)

def load():
    spec=importlib.util.spec_from_file_location('frozen_validation',ROOT/'presentation_assets/family_distribution_plot.py')
    mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
    d,m,families,h=mod.load_validated()
    inputs=[ROOT/'results/records.tsv',ROOT/'biological_annotation/protein_annotations.tsv',
            ROOT/'pfam_annotation/Pfam-A.clans.tsv',BASE/'esm2_probe/uesolds_pairs.tsv',
            BASE/'dev_cohort/dev_cohort_manifest.tsv',ROOT/'results/predictions.tsv',
            BASE/'esm2_probe/emb/uesolds_extraction_metadata.json',ROOT/'source_audit_results/protocol.json']
    h={str(ROOT/k):v for k,v in h.items()};h.update({str(p):sha(p) for p in inputs})
    identity=['raw_uid','legacy_uid','sha256','original_id','label','split','center','length','cluster']
    bio=read(inputs[1]).sort_values('raw_uid').reset_index(drop=True)
    pd.testing.assert_frame_equal(d[identity],bio[identity])
    records=read(inputs[0]).set_index('raw_uid').loc[d.raw_uid].reset_index()
    for c in ['sha256','label','split','center','length']:assert np.array_equal(records[c],d[c])
    for c in ['non_m_start','repeated_identifier','terminal_his','terminal_cleavage','terminal_linker','charge_proxy','hydropathy','entropy']:
        d[c]=records[c].values
    for c in ['taxonomy_status','organism','taxonomy_id','full_disorder_fraction','full_aggregation_fraction','full_tm_fraction']:
        d[c]=bio[c].values
    for a in AA:d['frac_'+a]=records['frac_'+a].values
    d['pfam_annotated']=d.raw_uid.map(lambda u:bool(families.get(u,set())))
    d['taxonomy_assigned']=d.taxonomy_status.eq('assigned_exact_species_consensus')
    ref=pd.read_csv(inputs[2],sep='\t',header=None,names=['pfam','clan','clan_name','family_name','description'],keep_default_na=False)
    m=m[['raw_uid','pfam']].merge(d[['raw_uid','split','label','center','cluster']],on='raw_uid',validate='many_to_one').merge(ref,on='pfam',validate='many_to_one')
    pairs=read(inputs[3]);seq=dict(pairs.iloc[:,:2].itertuples(index=False,name=None))
    old=read(inputs[4]).set_index('record_uid')
    for r in d.itertuples():
        s=seq[r.legacy_uid];assert len(s)==r.length and hashlib.sha256(s.encode()).hexdigest()==r.sha256
        assert old.loc[r.legacy_uid,'sha256']==r.sha256 and old.loc[r.legacy_uid,'label']==r.label
    t=d[d.split=='train'];v=d[d.split=='validation']
    for c in ['sha256','cluster','original_id']:assert not set(t[c])&set(v[c])
    assert d.cluster.ne('').all()
    return d,m,families,seq,h

def support(m,keys):
    if m.empty:return pd.DataFrame(columns=keys+['n0','n1','groups0','groups1'])
    grouped=m.groupby(keys+['label']).agg(n=('raw_uid','nunique'),groups=('cluster','nunique')).unstack('label',fill_value=0)
    for metric in ['n','groups']:
        for label in [0,1]:
            if (metric,label) not in grouped.columns:grouped[(metric,label)]=0
    grouped=grouped[[('n',0),('n',1),('groups',0),('groups',1)]]
    grouped.columns=['n0','n1','groups0','groups1']
    return grouped.reset_index()


def supported_sets(m,keys,k):
    z=support(m,keys);return set(map(tuple,z.loc[(z.n0>=k)&(z.n1>=k),keys].values))

def summary(g,m):
    q=m[m.raw_uid.isin(g.raw_uid)];f=support(q,['pfam']);c=support(q,['pfam','center'])
    out=dict(n=len(g),n0=int(g.label.eq(0).sum()),n1=int(g.label.eq(1).sum()),groups=g.cluster.nunique(),
             groups0=g[g.label==0].cluster.nunique(),groups1=g[g.label==1].cluster.nunique(),sources=g.center.nunique(),
             pfam_annotated=int(g.pfam_annotated.sum()),families=q.pfam.nunique(),clans=q.loc[q.clan!='','clan'].nunique(),
             mixed_families=int(((f.n0>0)&(f.n1>0)).sum()),mixed_family_source_pairs=int(((c.n0>0)&(c.n1>0)).sum()),
             taxa=g.loc[g.taxonomy_assigned,'taxonomy_id'].nunique(),taxonomy_assigned=int(g.taxonomy_assigned.sum()),
             biochemical_cells=g.bio_cell.nunique())
    for k in [3,5]:
        out[f'families_ge{k}']=int(((f.n0>=k)&(f.n1>=k)).sum())
        out[f'family_source_pairs_ge{k}']=int(((c.n0>=k)&(c.n1>=k)).sum())
        for keys,name in [(['pfam'],'families'),(['pfam','center'],'family_source_pairs')]:
            tr=supported_sets(q[q.split=='train'],keys,k);va=supported_sets(q[q.split=='validation'],keys,k)
            out[f'{name}_both_splits_ge{k}']=len(tr&va)
    return out

def diversity(d,m,seq):
    features=np.c_[np.log(d.length),d.charge_proxy,d.hydropathy]
    bins=[np.unique(np.quantile(features[d.split=='train',j],[.2,.4,.6,.8])) for j in range(3)]
    d['bio_cell']=[':'.join(map(str,r)) for r in np.column_stack([np.digitize(features[:,j],bins[j]) for j in range(3)])]
    (OUT/'biochemical_bins.json').write_text(json.dumps(dict(features=['log_full_length','charge_proxy','KD_hydropathy'],training_boundaries=[b.tolist() for b in bins]),indent=2)+'\n')
    rows=[];marg=[];overlap=[];dist=[];desc=[];familyrows=[]
    overall=summary(d,m)
    for center,g in d.groupby('center'):
        rest=d[d.center!=center];q=m[m.center==center];other=m[m.center!=center]
        for split in ['train','validation','combined']:
            gg=g if split=='combined' else g[g.split==split]
            rows.append(dict(center=center,split=split,**summary(gg,m)))
        fs=set(q.pfam);cs=set(q.loc[q.clan!='','clan']);gs=set(g.cluster)
        kept=summary(rest,m)
        marg.append(dict(center=center,removed_n=len(g),unique_families_lost=len(fs-set(other.pfam)),
                         unique_clans_lost=len(cs-set(other.clan)),groups_lost=len(gs-set(rest.cluster)),
                         shared_families=len(fs&set(other.pfam)),exclusive_biochemical_cells_lost=len(set(g.bio_cell)-set(rest.bio_cell)),
                         **{'lost_'+k:overall[k]-kept[k] for k in ['mixed_families','mixed_family_source_pairs','families_both_splits_ge3','family_source_pairs_both_splits_ge3','families_both_splits_ge5','family_source_pairs_both_splits_ge5','taxa']}))
        for level,col in [('family','pfam'),('clan','clan')]:
            qq=q if col=='pfam' else q[q.clan!='']
            counts=qq.groupby(col).raw_uid.nunique()
            for cat,count in counts.items():dist.append(dict(center=center,level=level,accession=cat,proteins=count,denominator=len(g),fraction=count/len(g)))
        for (split,label),gg in g.groupby(['split','label']):
            for f in ['length','charge_proxy','hydropathy','entropy']+['frac_'+a for a in AA]+['terminal_his','terminal_cleavage','terminal_linker','non_m_start','repeated_identifier','pfam_annotated','taxonomy_assigned','full_disorder_fraction','full_aggregation_fraction','full_tm_fraction']:
                x=pd.to_numeric(gg[f],errors='coerce').astype(float)
                desc.append(dict(center=center,split=split,label=label,feature=f,n=len(gg),observed=x.notna().sum(),missing=x.isna().sum(),mean=x.mean(),q05=x.quantile(.05),median=x.median(),q95=x.quantile(.95)))
    centers=sorted(d.center.unique())
    for i,c in enumerate(centers):
        for b in centers[i+1:]:
            x=d[d.center==c];y=d[d.center==b];fx=set(m[m.center==c].pfam);fy=set(m[m.center==b].pfam)
            gx=set(x.cluster);gy=set(y.cluster)
            overlap.append(dict(source_a=c,source_b=b,shared_families=len(fx&fy),family_jaccard=len(fx&fy)/len(fx|fy) if fx|fy else np.nan,
                                a_only_families=len(fx-fy),b_only_families=len(fy-fx),shared_groups=len(gx&gy),
                                exact_sequence_overlap=len(set(x.sha256)&set(y.sha256)),shared_bio_cells=len(set(x.bio_cell)&set(y.bio_cell))))
    save(pd.DataFrame(rows),'source_diversity.tsv');save(pd.DataFrame(marg),'leave_source_out_diversity.tsv')
    save(pd.DataFrame(overlap),'source_pair_overlap.tsv');save(pd.DataFrame(dist),'source_family_clan_distribution.tsv')
    save(pd.DataFrame(desc),'source_label_descriptors.tsv')
    save(support(m,['split','center','pfam']),'family_source_support.tsv')
    save(d.groupby(['center','split','label','taxonomy_status']).size().reset_index(name='n'),'taxonomy_coverage.tsv')
    policies={'A_all':np.ones(len(d),bool),'B_DNASU_candidate':d.center.eq('DNASU_clone'),
              'C_canonical':np.ones(len(d),bool),'C_verified_property_core':np.zeros(len(d),bool),
              'D_JCSG':d.center.eq('JCSG'),'D_MCSG':d.center.eq('MCSG'),
              'D_JCSG_MCSG':d.center.isin(['JCSG','MCSG']),'D_training_supported_centres':d.source_controlled}
    ps=[];manifest=d[['raw_uid','legacy_uid','sha256','split','label','center','cluster']].copy()
    for policy,mask in policies.items():
        manifest[policy]=mask
        for split in ['train','validation','combined']:
            g=d[mask & (True if split=='combined' else d.split.eq(split))]
            r=summary(g,m);r.update(policy=policy,split=split,families_lost=overall['families']-r['families'] if split=='combined' else '',
                                       groups_lost=overall['groups']-r['groups'] if split=='combined' else '')
            ps.append(r)
    manifest['concept_eligible']=True;manifest['published_label_exploratory_only']=True
    manifest['verified_property_eligible']=False
    manifest['property_eligibility_reason']='No recovered construct-matched comparable binary assay evidence; candidate flags do not certify labels'
    save(manifest,'policy_manifest.tsv');save(pd.DataFrame(ps),'policy_comparison.tsv')
    save(pd.concat([d[mask].groupby(['split','center','label']).size().reset_index(name='n').assign(policy=policy) for policy,mask in policies.items()]),'policy_source_counts.tsv')
    sample=[]
    for _,g in d[d.center.isin(['JCSG','MCSG'])].groupby(['center','split','label']):
        g=g.sort_values(['length','raw_uid']);sample.append(g.iloc[[0,len(g)//2,len(g)-1]])
    s=pd.concat(sample).sort_values(['center','split','label','length'])
    save(s[['raw_uid','legacy_uid','original_id','split','label','center','length','sha256','terminal_his','terminal_cleavage','non_m_start']],'recovery_sample.tsv')
    # Descriptive figures only. All family memberships retained, no category balancing.
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    tab=pd.DataFrame(marg).set_index('center').sort_values('unique_families_lost')
    fig,ax=plt.subplots(figsize=(8,5));ax.barh(tab.index,tab.unique_families_lost,color='#829FAF');ax.set_xlabel('Observed Pfam families lost if source removed');fig.tight_layout();fig.savefig(OUT/'marginal_family_diversity.png',dpi=180);plt.close(fig)
    return d

def bootstrap_indices(groups,B=BOOT):
    _,inv=np.unique(np.asarray(groups),return_inverse=True);members=[np.flatnonzero(inv==k) for k in range(inv.max()+1)]
    rng=np.random.default_rng(42)
    return [np.concatenate([members[k] for k in rng.integers(0,len(members),len(members))]) for _ in range(B)] if members else []

def binary(y,p):
    y=np.asarray(y);p=np.asarray(p)
    if not len(y):return dict(auroc=np.nan,ap=np.nan,brier=np.nan,logloss=np.nan)
    return dict(auroc=roc_auc_score(y,p) if len(set(y))==2 else np.nan,
                ap=average_precision_score(y,p) if y.sum()>0 else np.nan,
                brier=np.mean((p-y)**2),logloss=log_loss(y,np.c_[1-p,p],labels=[0,1]))

def binary_summary(g,p,ci=True):
    y=g.label.to_numpy();out=dict(n=len(g),n0=int((y==0).sum()),n1=int(y.sum()),groups=g.cluster.nunique(),prevalence=float(y.mean()),**binary(y,p))
    if ci:
        draws=[binary(y[ix],p[ix]) for ix in bootstrap_indices(g.cluster)]
        for key in ['auroc','ap','brier','logloss']:
            z=[r[key] for r in draws if np.isfinite(r[key])]
            out[key+'_valid_bootstraps']=len(z)
            out[key+'_low'],out[key+'_high']=np.quantile(z,[.025,.975]) if len(z)>=20 else [np.nan,np.nan]
    return out

def fit_predict(x,xv,y,scale=True):
    if scale:
        s=StandardScaler().fit(x);x=s.transform(x);xv=s.transform(xv)
    assert np.isfinite(x.data if sparse.issparse(x) else x).all()
    with warnings.catch_warnings():
        warnings.simplefilter('error',ConvergenceWarning)
        fit=LogisticRegression(C=1,solver='lbfgs',max_iter=5000,tol=1e-4,random_state=42).fit(x,y)
    return fit.predict_proba(xv),fit.classes_,int(max(fit.n_iter_))

def multi(y,p,nclass):
    guess=p.argmax(axis=1);cm=confusion_matrix(y,guess,labels=np.arange(nclass));tp=cm.diagonal();actual=cm.sum(1);pred=cm.sum(0);present=actual>0
    recall=np.divide(tp,actual,out=np.zeros(nclass),where=actual>0);f1=np.divide(2*tp,actual+pred,out=np.zeros(nclass),where=(actual+pred)>0)
    auc=[roc_auc_score(y==k,p[:,k]) for k in range(nclass) if 0<(y==k).sum()<len(y)]
    ap=[average_precision_score(y==k,p[:,k]) for k in range(nclass) if 0<(y==k).sum()<len(y)]
    return dict(accuracy=float((guess==y).mean()),balanced_accuracy=float(recall[present].mean()),macro_f1=float(f1[present].mean()),
                macro_ovr_auroc=float(np.mean(auc)),macro_ap=float(np.mean(ap)),logloss=log_loss(y,p,labels=np.arange(nclass)),represented_classes=int(present.sum()))

def probes(d,m,families,seq,hashes):
    vectors={}
    for path in sorted((BASE/'esm2_probe/emb').glob('uesolds_shard*.npz')):
        hashes[str(path)]=sha(path)
        with np.load(path,allow_pickle=False) as z:
            for key,vec in zip(z['keys'],z['whole']):
                assert str(key) not in vectors;vectors[str(key)]=vec
    assert len(vectors)==14000
    meta=json.loads((BASE/'esm2_probe/emb/uesolds_extraction_metadata.json').read_text())
    assert meta['layer_index']==8 and meta['max_residues']==1022 and meta['frozen']
    esm=np.stack([vectors[k] for k in d.legacy_uid]).astype(float);assert esm.shape==(13880,1280) and np.isfinite(esm).all()
    c=np.array([[s.count(a)/len(s) for a in AA]+[np.log(len(s))] for s in [seq[k][:1022] for k in d.legacy_uid]])
    train=np.flatnonzero(d.split.eq('train'));val=np.flatnonzero(d.split.eq('validation'));v=d.iloc[val].reset_index(drop=True);t=d.iloc[train]
    vocabulary=sorted(set(m[m.split=='train'].pfam));vb=set(vocabulary)
    enc=MultiLabelBinarizer(classes=vocabulary,sparse_output=True).fit([[]])
    fam=enc.transform([sorted(families.get(u,set())&vb) for u in d.raw_uid])
    fam=sparse.hstack([fam,sparse.csr_matrix((~d.pfam_annotated).astype(int).to_numpy()[:,None])],format='csr')
    source_names=sorted(t.center.unique());source_index={s:i for i,s in enumerate(source_names)}
    ys=d.center.map(source_index).to_numpy();yv=ys[val];nclass=len(source_names)
    features={'composition':c,'Pfam':fam,'ESM_layer8':esm,'composition_plus_ESM':np.c_[c,esm]}
    source_predictions={};source_results=[];source_per=[];fitrows=[];cmrows=[]
    prior=np.bincount(ys[train],minlength=nclass)/len(train)
    source_predictions['training_prior']=np.tile(prior,(len(val),1))
    draws=bootstrap_indices(v.cluster)
    for name,x in features.items():
        emit('Source probe '+name)
        pp,classes,it=fit_predict(x[train],x[val],ys[train],name!='Pfam')
        assert np.array_equal(classes,np.arange(nclass));source_predictions[name]=pp
        fitrows.append(dict(task='source',protocol='all_train',model=name,n_train=len(train),iterations=it,n_features=x.shape[1]))
    for name,p in source_predictions.items():
        metrics=multi(yv,p,nclass);boot=[multi(yv[ix],p[ix],nclass) for ix in draws]
        for metric,value in metrics.items():
            if metric=='represented_classes':continue
            lo,hi=np.quantile([b[metric] for b in boot],[.025,.975]);source_results.append(dict(model=name,metric=metric,value=value,low=lo,high=hi,n=len(v),classes_in_training=nclass,classes_in_validation=len(set(yv))))
        pred=p.argmax(axis=1)
        for k,s in enumerate(source_names):
            hit=yv==k;guessed=pred==k;n=int(hit.sum());pr=hit[pred==k].mean() if guessed.sum() else np.nan;rec=(pred[hit]==k).mean() if n else np.nan
            b=binary_summary(v.assign(label=hit.astype(int)),p[:,k])
            source_per.append(dict(model=name,center=s,support=n,predicted=int(guessed.sum()),precision=pr,recall=rec,**b))
        cm=confusion_matrix(yv,pred,labels=np.arange(nclass))
        for i,actual in enumerate(source_names):
            for j,predicted in enumerate(source_names):cmrows.append(dict(model=name,actual=actual,predicted=predicted,n=int(cm[i,j])))
        save(pd.DataFrame(p,columns=source_names).assign(raw_uid=v.raw_uid.to_numpy()),f'source_predictions_{name}.tsv')
    increments=[]
    for alt,base in [('ESM_layer8','composition'),('composition_plus_ESM','composition'),('ESM_layer8','Pfam')]:
        for metric in ['accuracy','balanced_accuracy','macro_f1','logloss']:
            differences=np.array([multi(yv[ix],source_predictions[alt][ix],nclass)[metric]-multi(yv[ix],source_predictions[base][ix],nclass)[metric] for ix in draws])
            point=multi(yv,source_predictions[alt],nclass)[metric]-multi(yv,source_predictions[base],nclass)[metric]
            lo,hi=np.quantile(differences,[.025,.975]);increments.append(dict(model=alt,reference=base,metric=metric,difference=point,low=lo,high=hi))
    # Same-family multiclass predictions; candidate source support is explicit.
    conditional=[];source_family_support=[]
    for f,g in m.groupby('pfam'):
        tc=g[g.split=='train'].groupby('center').raw_uid.nunique();vc=g[g.split=='validation'].groupby('center').raw_uid.nunique()
        eligible=sorted(s for s in tc.index if tc[s]>=3 and vc.get(s,0)>=3)
        if len(eligible)<2:continue
        ids=set(g[(g.split=='validation')&g.center.isin(eligible)].raw_uid);ix=np.flatnonzero(v.raw_uid.isin(ids));y=yv[ix]
        for s in eligible:source_family_support.append(dict(pfam=f,center=s,train_n=int(tc[s]),validation_n=int(vc[s])))
        pri=np.zeros(nclass)
        for s in eligible:pri[source_index[s]]=tc[s]
        pri/=pri.sum()
        for name,p in {**source_predictions,'family_training_prior':np.tile(pri,(len(val),1))}.items():
            # Restrict the candidate identities equally for all models, then renormalize.
            candidates=np.array([source_index[s] for s in eligible]);pp=p[ix][:,candidates];pp=pp/pp.sum(1,keepdims=True)
            yy=np.array([eligible.index(v.center.iloc[i]) for i in ix])
            conditional.append(dict(pfam=f,model=name,n=len(ix),centres=len(eligible),**multi(yy,pp,len(eligible))))
    save(pd.DataFrame(source_results),'source_prediction_metrics.tsv');save(pd.DataFrame(source_per),'source_prediction_per_centre.tsv')
    save(pd.DataFrame(cmrows),'source_confusion.tsv');save(pd.DataFrame(increments),'source_paired_differences.tsv')
    save(pd.DataFrame(conditional),'source_within_family.tsv');save(pd.DataFrame(source_family_support),'source_within_family_support.tsv')
    emit('Source probes and uncertainty complete')
    # Existing baseline probabilities are reused; new fits change training source composition only.
    old=read(ROOT/'results/predictions.tsv');old=old[old.cohort=='clean_all_sources']
    predictions={}
    for model in ['composition','ESM_layer8','source']:
        p=old[old.model==model].set_index('raw_uid').loc[v.raw_uid]
        assert np.array_equal(p.label,v.label);predictions[('all_train',model)]=p.prediction.to_numpy()
    for omitted in source_names:
        ii=train[d.iloc[train].center.ne(omitted).to_numpy()]
        for model,x in [('composition',c),('ESM_layer8',esm)]:
            emit('Solubility leave out '+omitted+' '+model)
            pp,cls,it=fit_predict(x[ii],x[val],d.iloc[ii].label.to_numpy());assert list(cls)==[0,1]
            predictions[('without_'+omitted,model)]=pp[:,1]
            fitrows.append(dict(task='solubility',protocol='without_'+omitted,model=model,n_train=len(ii),iterations=it,n_features=x.shape[1]))
    for center in ['JCSG','MCSG']:
        ii=train[d.iloc[train].center.eq(center).to_numpy()]
        for model,x in [('composition',c),('ESM_layer8',esm)]:
            emit('Solubility train only '+center+' '+model)
            pp,cls,it=fit_predict(x[ii],x[val],d.iloc[ii].label.to_numpy());predictions[('only_'+center,model)]=pp[:,1]
            fitrows.append(dict(task='solubility',protocol='only_'+center,model=model,n_train=len(ii),iterations=it,n_features=x.shape[1]))
    metrics=[];predrows=[];deltas=[]
    strata=[('ALL',np.arange(len(v)))]+[(s,np.flatnonzero(v.center.eq(s))) for s in sorted(v.center.unique())]
    for (protocol,model),p in predictions.items():
        predrows.append(pd.DataFrame(dict(raw_uid=v.raw_uid,protocol=protocol,model=model,probability=p)))
        for source,ix in strata:
            # CIs for pooled and the removed/sole training source; all-source comparison has every source CI.
            ci=protocol=='all_train' or source=='ALL' or protocol in ['without_'+source,'only_'+source]
            metrics.append(dict(protocol=protocol,model=model,evaluation=source,**binary_summary(v.iloc[ix],p[ix],ci)))
    for source,ix in strata:
        g=v.iloc[ix];bootix=bootstrap_indices(g.cluster)
        for protocol in sorted({key[0] for key in predictions}):
            pairs=[('ESM_layer8','composition',protocol,'ESM_minus_composition')] if protocol=='all_train' or source=='ALL' or protocol in ['without_'+source,'only_'+source] else []
            if protocol!='all_train':pairs += [(model,model,'all_train','training_change') for model in ['composition','ESM_layer8']]
            for alt,base,baseproto,contrast in pairs:
                # Source omission deltas reported on fixed pooled and omitted-source populations; cross-centre matrix also saved above.
                if contrast=='training_change' and source!='ALL' and protocol not in ['without_'+source,'only_'+source]:continue
                pa=predictions[(protocol,alt)][ix];pb=predictions[(baseproto,base)][ix];y=g.label.to_numpy()
                for metric in ['auroc','ap']:
                    av=binary(y,pa)[metric];bv=binary(y,pb)[metric]
                    if not np.isfinite(av+bv):continue
                    diffs=[binary(y[j],pa[j])[metric]-binary(y[j],pb[j])[metric] for j in bootix]
                    z=np.array(diffs);z=z[np.isfinite(z)];lo,hi=np.quantile(z,[.025,.975]) if len(z)>=20 else [np.nan,np.nan]
                    deltas.append(dict(protocol=protocol,model=alt,reference_protocol=baseproto,reference_model=base,contrast=contrast,evaluation=source,metric=metric,difference=av-bv,low=lo,high=hi,valid_bootstraps=len(z),n=len(ix)))
    familymetrics=[];feas=[]
    for keys,kind in [(['pfam'],'family'),(['pfam','center'],'family_source')]:
        trset=supported_sets(m[m.split=='train'],keys,3)
        for key,g in m[m.split=='validation'].groupby(keys):
            key=key if isinstance(key,tuple) else (key,)
            if key not in trset:continue
            ix=np.flatnonzero(v.raw_uid.isin(set(g.raw_uid)));vg=v.iloc[ix]
            counts=vg.label.value_counts();ok=min(counts.get(0,0),counts.get(1,0))>=3
            feas.append(dict(kind=kind,stratum='|'.join(key),n=len(vg),n0=counts.get(0,0),n1=counts.get(1,0),validation_ge3_each=ok))
            if not ok:continue
            for model in ['composition','ESM_layer8','source']:
                familymetrics.append(dict(kind=kind,stratum='|'.join(key),model=model,**binary_summary(vg,predictions[('all_train',model)][ix])))
    save(pd.DataFrame(metrics),'solubility_metrics.tsv');save(pd.concat(predrows),'solubility_predictions.tsv')
    save(pd.DataFrame(deltas),'solubility_paired_differences.tsv');save(pd.DataFrame(familymetrics),'solubility_family_strata.tsv')
    save(pd.DataFrame(feas),'solubility_family_feasibility.tsv');save(pd.DataFrame(fitrows),'model_fits.tsv')
    save(v[['raw_uid','sha256','split','label','center','cluster']],'evaluation_manifest.tsv')
    emit('All probes complete')

def main(stage):
    start=time.time();d,m,f,seq,h=load();emit('Frozen cohort and identities verified')
    if stage in ['all','diversity']:d=diversity(d,m,seq);emit('Source diversity and policy tables complete')
    if stage in ['all','probes']:
        with threadpool_limits(limits=2):probes(d,m,f,seq,h)
    assert all(sha(Path(p))==digest for p,digest in h.items())
    record=dict(stage=stage,python=platform.python_version(),numpy=np.__version__,pandas=pd.__version__,sklearn=sklearn.__version__,
                script_sha256=sha(Path(__file__)),seconds=time.time()-start,input_sha256=h,checks='passed; no independent-test access',
                cohort_n=13880,training_n=9880,validation_n=4000)
    (OUT/f'run_{stage}.json').write_text(json.dumps(record,indent=2)+'\n');emit('Input hashes unchanged')

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--stage',choices=['all','diversity','probes'],default='all');main(p.parse_args().stage)
