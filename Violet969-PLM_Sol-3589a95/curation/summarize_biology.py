#!/usr/bin/env python3
"""Descriptive annotation QC/EDA. No predictive model, tuning or feature selection."""
import argparse,json,os
from pathlib import Path
os.environ.setdefault('MPLCONFIGDIR','/tmp/uesolds-matplotlib')
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from annotate_biology import OUT,ROOT,IDENTITY,read,save,sha,write_json,load,load_disorder,check_inputs
from biological_features import AA,regions,clip_regions,regional_summary,biochemistry,aggregation

# Prespecified concept representatives, independent of association strength.
FEATURES=['full_length','full_entropy_bits','full_charge_density_ph7p0','full_pi_estimate',
          'full_frac_charged','full_frac_polar_uncharged','full_frac_aromatic','full_hydropathy_mean',
          'full_hydropathy_w9_max','full_hydropathy_w19_max','full_aggregation_fraction',
          'full_aggregation_mean','full_disorder_fraction','full_disorder_mean']

def group_summary(d,keys,features):
    long=d.melt(id_vars=keys,value_vars=features,var_name='feature',value_name='value')
    return long.groupby(keys+['feature'],dropna=False).value.agg(n='count',mean='mean',sd='std',median='median',minimum='min',maximum='max').reset_index()

def effects(d,keys,features,min_class=3):
    rows=[]
    groups=d.groupby(keys,dropna=False) if keys else [((),d)]
    for name,g in groups:
        if not isinstance(name,tuple):name=(name,)
        a=g[g.label==0];b=g[g.label==1]
        if min(len(a),len(b))<min_class:continue
        for f in features:
            x=pd.to_numeric(a[f],errors='coerce').dropna();y=pd.to_numeric(b[f],errors='coerce').dropna()
            sd=np.sqrt((x.var()+y.var())/2)
            rows.append(dict(zip(keys,name))|dict(feature=f,n0=len(x),n1=len(y),mean0=x.mean(),mean1=y.mean(),
                difference=y.mean()-x.mean(),standardized_difference=(y.mean()-x.mean())/sd if sd>0 else np.nan,
                groups0=a.cluster.nunique(),groups1=b.cluster.nunique()))
    return pd.DataFrame(rows)

def validate(out):
    original,seqs=load();d=pd.read_csv(out/'protein_annotations.tsv',sep='\t',keep_default_na=False)
    pd.testing.assert_frame_equal(d[IDENTITY],original[IDENTITY])
    for c in ['pfam_accessions','pfam_status','pfam_release','label_reliability']:
        assert d[c].astype(str).tolist()==original[c].astype(str).tolist(),c
    checks={'record_count':len(d),'identity_columns_unchanged':True,'input_hashes_unchanged':True}
    check_inputs(out)
    legacy=read(ROOT/'results/records.tsv').set_index('raw_uid').loc[d.raw_uid]
    for new,old in [('full_charge_proxy','charge_proxy'),('full_hydropathy_mean','hydropathy'),('full_entropy_bits','entropy')]+[('full_frac_'+a,'frac_'+a) for a in AA]:
        np.testing.assert_allclose(d[new].astype(float).values,legacy[old].astype(float).values,rtol=2e-7,atol=1e-8)
    checks['existing_sequence_descriptors_reproduced']=True
    statuses=read(out/'annotation_status.tsv')
    assert not statuses.duplicated(['raw_uid','concept']).any()
    assert statuses.groupby('concept').size().eq(13880).all()
    reg=read(out/'regions.tsv');lengths=d.set_index('raw_uid').length
    assert ((reg.start>=1)&(reg.end>=reg.start)&(reg.end<=reg.raw_uid.map(lengths))).all()
    pdom=read(out/'pfam_domains.tsv')
    for start,end in [('ali_from','ali_to'),('env_from','env_to')]:
        assert ((pdom[start]>=1)&(pdom[end]>=pdom[start])&(pdom[end]<=pdom.raw_uid.map(lengths))).all()
    original_domains=pd.read_csv(ROOT/'pfam_annotation/domains.tsv',sep='\t',float_precision='round_trip')
    preserved_domains=pd.read_csv(out/'pfam_domains.tsv',sep='\t',float_precision='round_trip')
    pd.testing.assert_frame_equal(original_domains,preserved_domains[original_domains.columns])
    assert d.groupby('cluster').split.nunique().max()==1
    dis=load_disorder(out)
    with np.load(out/'aggregation_profiles.npz',allow_pickle=False) as z:ag={k:z[k] for k in z.files}
    with np.load(out/'sequence_profiles.npz',allow_pickle=False) as z:sq={k:z[k] for k in z.files}
    assert ag['raw_uid'].tolist()==d.raw_uid.tolist()==sq['raw_uid'].tolist()
    assert np.array_equal(ag['offsets'],sq['offsets'])
    assert len(dis)==len(d)
    mismatch=0;max_error=0.
    for i,r in enumerate(d.itertuples()):
        start,end=ag['offsets'][i:i+2];s=seqs[r.raw_uid]
        assert bytes(sq['aa_ascii'][start:end]).decode()==s
        av=ag['scores'][start:end];dv,status,h=dis[r.raw_uid]
        assert len(av)==len(dv)==len(s) and h==r.sha256
        # Recompute scalar summaries from retained positional arrays.
        for concept,v,threshold,minsize in [('disorder',dv,.5,1),('aggregation',av,-.02,5)]:
            mask=v>=threshold if concept=='disorder' else (v>threshold)&np.array([a!='P' for a in s])
            rr=regions(mask,minsize)
            for scope,limit in [('full',len(s)),('esm',min(1022,len(s)))]:
                summary=regional_summary(v[:limit],clip_regions(rr,limit))
                for key in ['mean','maximum','fraction','n_regions','longest']:
                    saved=float(getattr(r,f'{scope}_{concept}_{key}'));expected=summary[key]
                    error=abs(saved-expected);max_error=max(max_error,error)
                    if not np.isclose(saved,expected,atol=1e-6):mismatch+=1
    assert mismatch==0,mismatch
    checks.update(positional_records_checked=len(d),domain_coordinates_checked=len(pdom),regional_intervals_checked=len(reg),
                  domains_preserved_without_score_loss=True,detected_homology_groups_disjoint=True,
                  scalar_position_max_error=max_error,all_records_accounted_for=True,
                  test_access='No test file is referenced or read by these scripts.')
    write_json(out/'validation_checks.json',checks)

def run(out):
    validate(out)
    d=pd.read_csv(out/'protein_annotations.tsv',sep='\t',keep_default_na=False,na_values=['NA'])
    for f in ['full_tm_fraction','full_membrane_associated_fraction']:
        if f in d and d[f].notna().all() and f not in FEATURES:FEATURES.append(f)
    st=read(out/'annotation_status.tsv').merge(d[['raw_uid','split','label','center','length','terminal_his','terminal_cleavage','terminal_linker']],on='raw_uid',validate='many_to_one')
    st['length_bin']=pd.cut(st.length,[0,100,300,600,1022,np.inf],labels=['1-100','101-300','301-600','601-1022','>1022']).astype(str)
    for key,groups in [('overall',['concept','status']),('source_label_split',['split','center','label','concept','status']),
                       ('length',['length_bin','concept','status']),('construct',['terminal_his','terminal_cleavage','terminal_linker','concept','status'])]:
        save(st.groupby(groups,dropna=False).size().rename('n').reset_index(),out/f'coverage_{key}.tsv')
    member=read(ROOT/'pfam_annotation/memberships.tsv')[['raw_uid','pfam']].drop_duplicates()
    fam=d.merge(member,on='raw_uid',how='inner',validate='one_to_many')
    fs=st.merge(member,on='raw_uid',how='inner')
    save(fs.groupby(['pfam','split','label','concept','status']).size().rename('n').reset_index(),out/'coverage_family.tsv')
    features=FEATURES+['full_frac_'+a for a in AA]
    save(group_summary(d,['split','label'],features),out/'distributions_split_label.tsv')
    save(group_summary(d,['split','center','label'],features),out/'distributions_source_label.tsv')
    save(effects(d,['split'],features),out/'label_associations_pooled.tsv')
    save(effects(d,['split','center'],features),out/'label_associations_source.tsv')
    save(effects(fam,['split','pfam'],FEATURES),out/'label_associations_family.tsv')
    save(effects(fam,['split','pfam','center'],FEATURES),out/'label_associations_family_source.tsv')
    n=fam.groupby(['split','pfam']).raw_uid.transform('size')
    stable=fam[n>=5]
    within=group_summary(stable,['split','pfam'],FEATURES)
    global_sd=d.groupby('split')[FEATURES].std().stack().rename('global_sd').reset_index();global_sd.columns=['split','feature','global_sd']
    within=within.merge(global_sd,on=['split','feature']);within['within_sd_over_global_sd']=within.sd/within.global_sd
    save(within,out/'family_biochemical_variation.tsv')
    support=fam.groupby(['split','pfam']).agg(n=('raw_uid','nunique'),groups=('cluster','nunique'),sources=('center','nunique'),label_classes=('label','nunique')).reset_index()
    save(support,out/'family_biochemical_support.tsv')
    # Between/within-source descriptive variance, no adjustment or residualization.
    partitions=[]
    for split,g in d.groupby('split'):
        for f in FEATURES:
            overall=g[f].mean();ss=((g[f]-overall)**2).sum()
            between=sum(len(h)*(h[f].mean()-overall)**2 for _,h in g.groupby('center'))
            partitions.append(dict(split=split,feature=f,between_source_fraction=between/ss if ss else np.nan))
    save(pd.DataFrame(partitions),out/'source_variance.tsv')
    corr=d[FEATURES].corr(method='spearman');corr.to_csv(out/'concept_spearman.tsv',sep='\t',float_format='%.8g')
    fig,axes=plt.subplots(2,3,figsize=(12,7))
    plot_features=['full_length','full_charge_density_ph7p0','full_pi_estimate','full_hydropathy_mean','full_aggregation_fraction','full_disorder_fraction']
    for ax,f in zip(axes.flat,plot_features):
        for label,color in [(0,'#a84c45'),(1,'#356995')]:
            v=d.loc[d.label==label,f].dropna()
            ax.hist(v,bins=40,density=True,alpha=.5,color=color,label=f'label {label} (n={len(v)})')
        ax.set_title(f.removeprefix('full_'));ax.set_ylabel('Density')
        if f=='full_length':ax.set_yscale('log')
    axes.flat[0].legend(fontsize=8);fig.suptitle('Descriptive development distributions; labels include workflow-derived outcomes')
    fig.tight_layout();fig.savefig(out/'concept_distributions.png',dpi=160);plt.close(fig)
    means=d.groupby('center')[FEATURES].mean();z=(means-d[FEATURES].mean())/d[FEATURES].std()
    fig,ax=plt.subplots(figsize=(13,8));im=ax.imshow(z,vmin=-2,vmax=2,cmap='RdBu_r',aspect='auto')
    ax.set_yticks(range(len(z)),[f'{c} (n={(d.center==c).sum()})' for c in z.index]);ax.set_xticks(range(len(FEATURES)),[f.removeprefix('full_') for f in FEATURES],rotation=60,ha='right')
    ax.set_title('Source means in units of pooled protein SD; no source filtering');fig.colorbar(im,ax=ax,label='Standardized source mean')
    fig.tight_layout();fig.savefig(out/'source_concepts.png',dpi=160);plt.close(fig)
    fig,ax=plt.subplots(figsize=(11,9));im=ax.imshow(corr,vmin=-1,vmax=1,cmap='RdBu_r')
    labels=[f.removeprefix('full_') for f in FEATURES];ax.set_xticks(range(len(labels)),labels,rotation=60,ha='right');ax.set_yticks(range(len(labels)),labels)
    ax.set_title('Spearman relationships: sequence, aggregation and learned disorder');fig.colorbar(im,ax=ax)
    fig.tight_layout();fig.savefig(out/'concept_relationships.png',dpi=160);plt.close(fig)
    fig,ax=plt.subplots(figsize=(12,5))
    values=[within.loc[within.feature==f,'within_sd_over_global_sd'].dropna().values for f in FEATURES]
    ax.boxplot(values,showfliers=False,tick_labels=labels);ax.axhline(1,color='gray',linestyle='--');ax.tick_params(axis='x',rotation=60)
    ax.set_ylabel('Within-family SD / pooled split SD');ax.set_title('Families with ≥5 proteins per split; domain memberships overlap')
    fig.tight_layout();fig.savefig(out/'within_family_variation.png',dpi=160);plt.close(fig)
    # Machine-readable dictionary includes all emitted columns, units and scope.
    dictionary=[]
    for c in d:
        evidence='metadata'
        if c.startswith(('full_','esm_')):evidence='direct_sequence_derived'
        if 'aggregation' in c:evidence='biophysical_sequence_propensity'
        if 'disorder' in c:evidence='pretrained_model_prediction'
        if any(x in c for x in ['domain','pfam','taxonomy','organism','cluster','reference_accession','match_']):evidence='external_biological_annotation'
        if any(x in c for x in ['helix','sheet','coil','rsa_','exposed','tm_','signal_peptide']):evidence='pretrained_model_prediction'
        unit='dimensionless'
        if c.endswith(('_start','_end','_cleavage_site')):unit='1-based inclusive residue coordinate'
        elif 'net_charge' in c or 'original_context_charge' in c:unit='elementary charge'
        elif 'charge_density' in c or 'sidechain_charge' in c:unit='elementary charge per residue'
        elif 'pi_estimate' in c:unit='pH'
        elif 'entropy' in c:unit='bits'
        elif 'hydropathy' in c:unit='Kyte-Doolittle scale'
        dictionary.append(dict(column=c,dtype=str(d[c].dtype),evidence=evidence,unit=unit,scope='full' if c.startswith('full_') else 'first_1022' if c.startswith('esm_') else 'record',missing=int(d[c].isna().sum())))
    save(pd.DataFrame(dictionary),out/'feature_dictionary.tsv')
    summary=dict(records=len(d),splits=d.split.value_counts().to_dict(),labels=d.groupby(['split','label']).size().to_dict().__str__(),
        truncation_records=int(d.esm_truncated.sum()),regions_by_concept=read(out/'regions.tsv').concept.value_counts().to_dict(),
        aggregation_terminal_sensitive_records=int(d.aggregation_terminal_sensitivity_regions_changed.sum()),
        aggregation_length_extrapolation_records=int(d.aggregation_length_extrapolation.sum()),
        pfam_assigned=int(d.pfam_status.eq('assigned').sum()),
        taxonomy_statuses=d.taxonomy_status.value_counts().to_dict(),
        qualified_family_split_strata=int(len(support[support.n>=5])),
        same_source_family_label_contrasts=int(len(read(out/'label_associations_family_source.tsv'))/len(FEATURES)),
        features_prespecified=FEATURES)
    write_json(out/'eda_summary.json',summary)
    write_json(out/'output_hashes.json',{p.name:sha(p) for p in sorted(out.iterdir()) if p.is_file() and p.name!='output_hashes.json'})
    print(json.dumps(summary,indent=2))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--out',type=Path,default=OUT);a=p.parse_args();run(a.out)
