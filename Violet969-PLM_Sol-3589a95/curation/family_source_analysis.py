#!/usr/bin/env python3
"""Focused family/source audit and metadata probes. No test or embedding reads."""
import argparse,hashlib,json,os,platform,warnings
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import sparse
from sklearn.preprocessing import MultiLabelBinarizer,OneHotEncoder
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score,average_precision_score,log_loss,brier_score_loss
import sklearn
os.environ.setdefault('MPLCONFIGDIR','/tmp/uesolds-matplotlib')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parent

def read(p):return pd.read_csv(p,sep='\t',keep_default_na=False)
def save(d,p):d.to_csv(p,sep='\t',index=False,float_format='%.17g')
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def count_groups(m,keys):
    """Count distinct proteins/groups, keeping Pfam and clusters separate."""
    cols=keys+['n0','n1','n','groups0','groups1','groups_shared','exclusive_groups0','exclusive_groups1']
    rows=[]
    for key,g in m.groupby(keys,sort=True):
        key=key if isinstance(key,tuple) else (key,)
        z=g[g.label==0];o=g[g.label==1];c0=set(z.cluster);c1=set(o.cluster)
        rows.append(dict(zip(keys,key),n0=z.raw_uid.nunique(),n1=o.raw_uid.nunique(),n=g.raw_uid.nunique(),groups0=len(c0),groups1=len(c1),groups_shared=len(c0&c1),exclusive_groups0=len(c0-c1),exclusive_groups1=len(c1-c0)))
    return pd.DataFrame(rows,columns=cols)

def decompose(m,first,second):
    g=m.copy();g['w']=1/g.groupby('raw_uid').raw_uid.transform('size');g['wy']=g.w*g.label
    total=g.w.sum();p=g.wy.sum()/total;variance=p*(1-p)
    a=g.groupby(first)[['w','wy']].sum();a['p']=a.wy/a.w
    cells=g.groupby([first,second])[['w','wy']].sum().reset_index();cells['p']=cells.wy/cells.w
    cells['parent_p']=cells[first].map(a.p)
    between=((a.w/total)*(a.p-p)**2).sum()
    conditional=((cells.w/total)*(cells.p-cells.parent_p)**2).sum()
    residual=((cells.w/total)*cells.p*(1-cells.p)).sum()
    assert np.isclose(between+conditional+residual,variance)
    return dict(n_unique=g.raw_uid.nunique(),positive_prevalence=p,total_variance=variance,first=first,second=second,between_first=between,within_first_between_second=conditional,within_both=residual,between_first_fraction=between/variance,within_first_between_second_fraction=conditional/variance,within_both_fraction=residual/variance)

def supported(t,k):return t[(t.n0>=k)&(t.n1>=k)]

def main(out):
    out.mkdir(parents=True,exist_ok=True)
    inputs=[ROOT/'pfam_annotation/annotated_cohort.tsv',ROOT/'pfam_annotation/memberships.tsv',ROOT/'pfam_annotation/domains_flagged.tsv',ROOT/'results/cohort_manifest.tsv']
    hashes={str(p.relative_to(ROOT)):sha(p) for p in inputs}
    d=read(inputs[0]);m=read(inputs[1]);dom=read(inputs[2]);pool=read(inputs[3]);pool=pool[pool.clean]
    assert len(d)==13880 and d.raw_uid.is_unique
    assert set(d.raw_uid)==set(pool[pool.cached].raw_uid)
    assert not set(d[d.split=='train'].cluster)&set(d[d.split=='validation'].cluster)
    assert not set(pool[pool.split=='train'].cluster)&set(d[d.split=='validation'].cluster)
    meta=d[['raw_uid','split','label','center','id_class','cluster']]
    m=m[['raw_uid','pfam']].drop_duplicates().merge(meta,on='raw_uid',validate='many_to_one')
    clean_m=dom[~dom.clan_competing_overlap][['raw_uid','pfam']].drop_duplicates().merge(meta,on='raw_uid',validate='many_to_one')
    fam=count_groups(m,['split','pfam']);tr=fam[fam.split=='train'];va=fam[fam.split=='validation']
    focus=sorted(set(supported(tr,3).pfam)&set(supported(va,3).pfam));assert len(focus)==69
    save(pd.DataFrame({'pfam':focus}),out/'focus_69_families.tsv')
    fm=m[m.pfam.isin(focus)]
    counts=count_groups(fm,['split','pfam','center'])
    grid=pd.MultiIndex.from_product([['train','validation'],focus,sorted(d.center.unique())],names=['split','pfam','center']).to_frame(index=False)
    counts=grid.merge(counts,how='left',on=['split','pfam','center']).fillna(0)
    counts['prevalence']=np.where(counts.n>0,counts.n1/counts.n,np.nan)
    save(counts,out/'focus_69_family_source_counts.tsv')
    save(count_groups(fm,['split','pfam']),out/'focus_69_family_counts.tsv')
    variation=[];decomp=[]
    for split,g in fm.groupby('split'):
        for first,second in [('pfam','center'),('center','pfam')]:decomp.append(dict(split=split,**decompose(g,first,second)))
        for f,v in g.groupby('pfam'):
            c=counts[(counts.split==split)&(counts.pfam==f)&(counts.n>0)];mix=c[(c.n0>0)&(c.n1>0)]
            variation.append(dict(split=split,pfam=f,n=v.raw_uid.nunique(),n_sources=len(c),mixed_source_cells=len(mix),cells_ge3_each=len(supported(c,3)),cells_ge5_each=len(supported(c,5)),proteins_in_mixed_source_cells=int(mix.n.sum()),largest_source_fraction=c.n.max()/c.n.sum(),min_source_prevalence=c.prevalence.min(),max_source_prevalence=c.prevalence.max()))
    save(pd.DataFrame(variation),out/'focus_69_variation.tsv');save(pd.DataFrame(decomp),out/'variance_decomposition.tsv')
    # Centre comparisons: entire annotated A, not just the 69 focus families.
    centers=[];center_support=[]
    for policy,mm in [('GA_all',m),('no_clan_competition',clean_m)]:
        for center in ['JCSG','MCSG']:
            a=d[d.center==center];b=mm[mm.center==center]
            for split in ['train','validation']:
                z=a[a.split==split];bm=b[b.split==split]
                t=count_groups(bm,['pfam']);center_support.append(t.assign(policy=policy,center=center,split=split))
                centers.append(dict(policy=policy,center=center,split=split,n=len(z),n0=int(z.label.eq(0).sum()),n1=int(z.label.eq(1).sum()),prevalence=z.label.mean(),annotated=bm.raw_uid.nunique(),coverage=bm.raw_uid.nunique()/len(z),groups=z.cluster.nunique(),groups0=z[z.label==0].cluster.nunique(),groups1=z[z.label==1].cluster.nunique(),families=len(t),mixed_families=int(((t.n0>0)&(t.n1>0)).sum()),families_ge3=len(supported(t,3)),families_ge5=len(supported(t,5)),verified_outcomes=0))
    save(pd.DataFrame(centers),out/'center_comparison.tsv');save(pd.concat(center_support),out/'center_family_support.tsv')
    # Larger training pool: direct evidence and explicit validation-only ceiling.
    expansion=[];ceilings=[];proxy=[]
    for center in ['ALL',*sorted(d.center.unique())]:
        full=pool[(pool.split=='train')&(True if center=='ALL' else pool.center.eq(center))]
        cached=full[full.cached];extra=full[~full.cached];known=set(cached.cluster)
        for scope,a in [('cached',cached),('additional',extra),('full',full)]:
            expansion.append(dict(center=center,scope=scope,n=len(a),n0=int(a.label.eq(0).sum()),n1=int(a.label.eq(1).sum()),groups=a.cluster.nunique(),groups0=a[a.label==0].cluster.nunique(),groups1=a[a.label==1].cluster.nunique(),groups_not_in_cached_center=len(set(a.cluster)-known),uniform_Pfam_annotated=int(a.cached.sum())))
        for policy,mm in [('GA_all',m),('no_clan_competition',clean_m)]:
            v=mm[(mm.split=='validation')&(True if center=='ALL' else mm.center.eq(center))]
            t=mm[(mm.split=='train')&(True if center=='ALL' else mm.center.eq(center))]
            vc=count_groups(v,['pfam']);tc=count_groups(t,['pfam'])
            for k in [3,5]:
                vi=set(supported(vc,k).pfam);ti=set(supported(tc,k).pfam)
                vg=vc[(vc.exclusive_groups0>=k)&(vc.exclusive_groups1>=k)]
                ceilings.append(dict(center=center,policy=policy,k=k,current_both_splits=len(vi&ti),validation_only_ceiling=len(vi),max_additional_families=len(vi-ti),validation_exclusive_groups_ceiling=len(vg),ceiling_note='potential only; additional training Pfam unknown'))
    # SAME known homology group is a follow-up clue, never a transferred Pfam call.
    links=m[(m.split=='train')&m.pfam.isin(focus)][['pfam','cluster']].drop_duplicates()
    extras=pool[(pool.split=='train')&~pool.cached]
    candidates=extras.merge(links,on='cluster',how='inner')
    save(candidates[['raw_uid','center','label','cluster','pfam']].assign(evidence='same detected group as annotated cached training protein; family unverified'),out/'unannotated_cluster_link_candidates.tsv')
    for (f,c),g in candidates.groupby(['pfam','center']):
        proxy.append(dict(pfam=f,center=c,candidate_n0=g[g.label==0].raw_uid.nunique(),candidate_n1=g[g.label==1].raw_uid.nunique(),candidate_groups=g.cluster.nunique(),new_global_groups=0))
    save(pd.DataFrame(proxy),out/'cluster_link_support_proxy.tsv')
    save(pd.DataFrame(expansion),out/'training_pool_expansion.tsv');save(pd.DataFrame(ceilings),out/'fixed_validation_ceilings.tsv')
    # Fixed-C metadata probes; all A training records, no family balancing.
    predictions=[];scores=[];fitmeta=[];coef=[];bootstrap=[]
    v=d[d.split=='validation'].copy().reset_index(drop=True);yv=v.label.to_numpy()
    vpos={u:i for i,u in enumerate(v.raw_uid)}
    groups=[('pooled','ALL',np.arange(len(v)))]
    for c,g in v.groupby('center'):groups.append(('source',c,g.index.to_numpy()))
    vf=fm[fm.split=='validation']
    for f,g in vf.groupby('pfam'):groups.append(('focus_family',f,np.array([vpos[u] for u in g.raw_uid.unique()])))
    for (f,c),g in vf.groupby(['pfam','center']):
        tg=fm[(fm.split=='train')&(fm.pfam==f)&(fm.center==c)]
        if tg.label.eq(0).sum()>=3 and tg.label.eq(1).sum()>=3:
            groups.append(('family_source_train_ge3',f+'|'+c,np.array([vpos[u] for u in g.raw_uid.unique()])))
    for c in ['JCSG','MCSG']:
        groups.append(('focus_union_center',c,np.array(sorted({vpos[u] for u in vf[vf.center==c].raw_uid}))))
    groups.append(('focus_union','69_families',np.array(sorted({vpos[u] for u in vf.raw_uid}))))
    # Multi-hot family models are constant only for the complete assignment set.
    sets=m.groupby('raw_uid').pfam.agg(lambda x:';'.join(sorted(set(x))))
    v['pfam_set']=v.raw_uid.map(sets).fillna('NO_GA_HIT')
    for (ps,c),g in v.groupby(['pfam_set','center']):
        if g.label.nunique()==2:groups.append(('exact_family_set_source',ps+'|'+c,g.index.to_numpy()))
    for protocol,held in [('published_split',None),('leave_JCSG_out','JCSG'),('leave_MCSG_out','MCSG')]:
        train=d[(d.split=='train')&(True if held is None else d.center.ne(held))].copy()
        evalmask=np.ones(len(v),bool) if held is None else v.center.eq(held).to_numpy()
        for policy,mm in [('GA_all',m),('no_clan_competition',clean_m)]:
            if held is not None and policy!='GA_all':continue
            mapping=mm.groupby('raw_uid').pfam.agg(lambda x:sorted(set(x))).to_dict()
            ts=[mapping.get(u,[]) for u in train.raw_uid];vs=[mapping.get(u,[]) for u in v.raw_uid]
            mlb=MultiLabelBinarizer(sparse_output=True).fit(ts);known=set(mlb.classes_)
            ft=mlb.transform(ts);fv=mlb.transform([[f for f in fs if f in known] for fs in vs])
            ft=sparse.hstack([ft,np.array([not fs for fs in ts])[:,None]],format='csr');fv=sparse.hstack([fv,np.array([not fs for fs in vs])[:,None]],format='csr')
            enc=OneHotEncoder(handle_unknown='ignore').fit(train[['id_class','center']])
            st=enc.transform(train[['id_class','center']]);sv=enc.transform(v[['id_class','center']])
            models={'family':(ft,fv),'source':(st,sv),'family_plus_source':(sparse.hstack([ft,st],format='csr'),sparse.hstack([fv,sv],format='csr'))}
            for name,(xt,xv) in models.items():
                fit=LogisticRegression(C=1,solver='lbfgs',max_iter=5000,random_state=42).fit(xt,train.label)
                assert fit.n_iter_[0]<5000
                p=fit.predict_proba(xv)[:,1];assert np.isfinite(p).all()
                fitmeta.append(dict(protocol=protocol,policy=policy,model=name,n_train=len(train),n_features=xt.shape[1],iterations=int(fit.n_iter_[0]),validation_unseen_family_proteins=sum(bool(set(fs)-known) for fs in vs)))
                for u,y,pp,c,cl in zip(v.raw_uid,yv,p,v.center,v.cluster):
                    if held is None or c==held:predictions.append(dict(protocol=protocol,policy=policy,model=name,raw_uid=u,label=y,center=c,cluster=cl,prediction=pp))
                modelrows=[]
                for typ,key,idx in groups:
                    idx=idx[evalmask[idx]]
                    if not len(idx):continue
                    y=yv[idx];defined=len(set(y))==2
                    row=dict(protocol=protocol,policy=policy,model=name,group_type=typ,group=key,n=len(idx),n0=int((y==0).sum()),n1=int(y.sum()),auroc=roc_auc_score(y,p[idx]) if defined else np.nan,AP=average_precision_score(y,p[idx]) if defined else np.nan,log_loss=log_loss(y,p[idx],labels=[0,1]),brier=brier_score_loss(y,p[idx]))
                    scores.append(row);modelrows.append(row)
                for typ in ['source','focus_family','family_source_train_ge3','exact_family_set_source']:
                    sub=[r for r in modelrows if r['group_type']==typ and r['n0']>=3 and r['n1']>=3]
                    scores.append(dict(protocol=protocol,policy=policy,model=name,group_type='macro',group=typ+'_val_ge3',n=len(sub),auroc=float(np.mean([r['auroc'] for r in sub])) if sub else np.nan,AP=float(np.mean([r['AP'] for r in sub])) if sub else np.nan))
                if held is None:
                    # Pooled CI: resample validation homology groups, preserve all members.
                    clusters=sorted(v.cluster.unique());indices=[np.flatnonzero(v.cluster.to_numpy()==cl) for cl in clusters]
                    rng=np.random.default_rng(42);vals=[]
                    for _ in range(200):
                        idx=np.concatenate([indices[i] for i in rng.integers(0,len(indices),len(indices))])
                        vals.append(roc_auc_score(yv[idx],p[idx]))
                    bootstrap.append(dict(policy=policy,model=name,replicates=200,unit='detected validation homology cluster',auroc=roc_auc_score(yv,p),low=np.quantile(vals,.025),high=np.quantile(vals,.975)))
                print(protocol,policy,name,'AUROC',roc_auc_score(yv[evalmask],p[evalmask]),flush=True)
    save(pd.DataFrame(scores),out/'baseline_group_results.tsv');save(pd.DataFrame(predictions),out/'baseline_predictions.tsv');save(pd.DataFrame(fitmeta),out/'baseline_fits.tsv');save(pd.DataFrame(bootstrap),out/'pooled_cluster_bootstrap.tsv')
    # Source-by-family prevalence heatmaps; gray cells have no observed proteins.
    fig,axes=plt.subplots(1,2,figsize=(13,15))
    source_order=sorted(d.center.unique())
    for ax,split in zip(axes,['train','validation']):
        mat=counts[counts.split==split].pivot(index='pfam',columns='center',values='prevalence').reindex(index=focus,columns=source_order)
        cmap=plt.get_cmap('RdBu_r').copy();cmap.set_bad('#dddddd')
        im=ax.imshow(mat,aspect='auto',vmin=0,vmax=1,cmap=cmap);ax.set_title(split);ax.set_yticks(range(len(focus)),focus,fontsize=6);ax.set_xticks(range(len(source_order)),source_order,rotation=90,fontsize=8)
    fig.colorbar(im,ax=axes.ravel().tolist(),label='Fraction label 1',shrink=.7);fig.savefig(out/'family_source_prevalence.png',dpi=170,bbox_inches='tight');plt.close(fig)
    assert all(sha(ROOT/p)==h for p,h in hashes.items())
    (out/'provenance.json').write_text(json.dumps(dict(inputs=hashes,script_sha256=sha(Path(__file__)),python=platform.python_version(),sklearn=sklearn.__version__,numpy=np.__version__,pandas=pd.__version__,C=1,seed=42,train_n=9880,validation_n=4000,focus_families=69,notes=['69-family focus is development-selected descriptive analysis, not pristine evaluation','All models train on heterogeneous A or explicitly declared source-transfer sensitivity','No new annotation, embedding extraction or test access'],checks='input hashes preserved; fixed cohort; detected homology groups disjoint; models converged'),indent=2)+'\n')
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--out',type=Path,default=ROOT/'family_source_results');a=ap.parse_args();main(a.out)
