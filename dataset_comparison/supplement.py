from audit import *
def main():
 r=read(OUT/'records.tsv.gz');m=read(OUT/'policy_groups.tsv.gz');m.label=pd.to_numeric(m.label,errors='coerce');r.label=pd.to_numeric(r.label,errors='coerce');r.exp=pd.to_numeric(r.exp,errors='coerce');r.sol=pd.to_numeric(r.sol,errors='coerce');mem=read(OUT/'family_memberships.tsv.gz')
 duplicates=[];dup_records=[]
 for name,q in r.groupby('dataset'):
  g=q.groupby('sha256').agg(records=('record_id','size'),labels=('label','nunique'))
  bad=set(g.index[g.labels>1]);dup_records.append(q[q.sha256.isin(g.index[g.records>1])].drop(columns='sequence').assign(conflicting_binary_labels=lambda x:x.sha256.isin(bad)))
  duplicates.append(dict(dataset=name,records=len(q),unique_sequences=q.sha256.nunique(),repeated_sequences=int((g.records>1).sum()),repeated_records=int(g.loc[g.records>1,'records'].sum()),conflicting_binary_sequences=len(bad),conflicting_binary_records=int(q.sha256.isin(bad).sum())))
 save(pd.DataFrame(duplicates),'duplicate_summary.tsv');save(pd.concat(dup_records),'duplicate_records.tsv')
 raw=r[r.dataset=='NESG_raw'];q=raw[~raw.sha256.duplicated(False)].copy();reasons=[]
 for name,pool in [('raw',raw),('unique',q)]:
  for tag,emin,smin,sneg in [('E45_S45_vs0',4,4,0),('E345_S45_vs0',3,4,0),('E45_S345_vs0',4,3,0),('E345_S345_vs0',3,3,0),('E45_S45_vs01',4,4,1)]:
   reason=pd.Series('included',index=pool.index);reason.loc[~((pool.sol>=smin)|(pool.sol<=sneg))]='intermediate_S';reason.loc[pool.exp<emin]='insufficient_expression_for_policy';reason.loc[pool.exp==0]='E_zero_solubility_uninterpretable';reason.loc[~pool.exp.isin(range(6))|~pool.sol.isin(range(6))]='noninteger_score_unresolved'
   for why,n in reason.value_counts().items():reasons.append(dict(pool=name,policy=tag,reason=why,N=n))
 save(pd.DataFrame(reasons),'NESG_exclusion_reasons.tsv');save(raw[(raw.exp%1!=0)|((raw.exp==0)&(raw.sol!=0))].drop(columns='sequence'),'NESG_score_anomalies.tsv')
 # Additional eligibility views, retaining all published records in canonical files.
 extras=[]
 for pol in ['PSI_original','PSI_NetSolP']:
  x=m[m.policy==pol];bad=x.groupby('sha256').label.nunique();x=x[~x.sha256.isin(bad.index[bad>1])].sort_values('record_id').drop_duplicates('sha256').copy();x.policy=pol+'_consistent_deduplicated';extras.append(x)
 x=m[m.policy=='NESG_NetSolP'].merge(r[['record_id','exp','sol']],on='record_id');x=x[x.exp.isin([4,5])&x.sol.isin([0,4,5])].drop(columns=['exp','sol']);x.policy='NESG_NetSolP_integer_rule';extras.append(x)
 psi_groups=set(m.loc[m.policy=='PSI_original','audit_cluster'])
 for pol in ['NESG_NetSolP','NESG_unique_E45_S45_vs0']:
  x=m[(m.policy==pol)&~m.audit_cluster.isin(psi_groups)].copy();x.policy=pol+'_no_PSI_component';extras.append(x)
 extras=pd.concat(extras,ignore_index=True);save(extras,'additional_eligibility.tsv')
 summary=[]
 for pol,q in extras.groupby('policy'):
  x=q.merge(mem,on='sha256');t=x.groupby(['pfam','label']).record_id.nunique().unstack(fill_value=0).reindex(columns=[0,1],fill_value=0);g=x.groupby(['pfam','label']).audit_cluster.nunique().unstack(fill_value=0).reindex(columns=[0,1],fill_value=0)
  row=dict(policy=pol,total=len(q),soluble=int((q.label==1).sum()),insoluble=int((q.label==0).sum()),families=x.pfam.nunique(),mixed_families=int(((t[0]>0)&(t[1]>0)).sum()),groups=q.audit_cluster.nunique())
  for k in [3,5]:row[f'families_groups_ge{k}']=int(((g[0]>=k)&(g[1]>=k)).sum())
  summary.append(row)
 save(pd.DataFrame(summary),'additional_policy_counts.tsv')
 # Family contrasts within the two known construct-vector strata.
 vectors=read(OUT/'PSI_construct_mapping.tsv').set_index('curated_record').vector
 original=r[r.dataset=='PSI_original'].set_index('record_id').vector
 vector_support=[];vector_summary=[]
 for pol in ['PSI_original','PSI_NetSolP']:
  q=m[m.policy==pol].copy();q['vector']=q.record_id.map(original if pol=='PSI_original' else vectors);x=q.merge(mem,on='sha256')
  for (vector,fam),g in x.groupby(['vector','pfam']):
   z=dict(policy=pol,vector=vector,pfam=fam,name=g.family_name.iloc[0])
   for lab in [0,1]:z[f'n{lab}']=int((g.label==lab).sum());z[f'g{lab}']=g.loc[g.label==lab,'audit_cluster'].nunique()
   vector_support.append(z)
  for (vector,lab),g in q.groupby(['vector','label']):vector_summary.append(dict(policy=pol,vector=vector,label=lab,N=len(g)))
 save(pd.DataFrame(vector_support),'PSI_family_vector_support.tsv');save(pd.DataFrame(vector_summary),'PSI_vector_counts.tsv')
 # Evaluate every deterministic fold without choosing the best one.
 partition_rows=[];family_rows=[]
 for pol in ['PSI_original','PSI_NetSolP','NESG_unique_E45_S45_vs0','NESG_NetSolP']:
  q=m[m.policy==pol];x=q.merge(mem,on='sha256')
  for fold in range(5):
   t=x[x.audit_fold!=fold].groupby(['pfam','label']).audit_cluster.nunique().unstack(fill_value=0).reindex(columns=[0,1],fill_value=0);v=x[x.audit_fold==fold].groupby(['pfam','label']).audit_cluster.nunique().unstack(fill_value=0).reindex(index=t.index,columns=[0,1],fill_value=0).fillna(0)
   for k in [3,5]:
    supported=(t[0]>=k)&(t[1]>=k);both=(v[0]>=1)&(v[1]>=1);strong=(v[0]>=k)&(v[1]>=k)
    partition_rows.append(dict(policy=pol,holdout_fold=fold,train_groups_per_class=k,train_qualifying_families=int(supported.sum()),holdout_both_labels=int((supported&both).sum()),holdout_groups_ge_threshold=int((supported&strong).sum())))
 save(pd.DataFrame(partition_rows),'five_fold_feasibility.tsv')
 # Preserve source-file names, plus explicit published fold-based family feasibility.
 q=m[m.policy=='PSI_NetSolP'].merge(r[['record_id','published_partition']],on='record_id');q.published_partition=pd.to_numeric(q.published_partition).astype(int);x=q.merge(mem,on='sha256');p=[]
 for fold in range(5):
  t=x[x.published_partition!=fold].groupby(['pfam','label']).audit_cluster.nunique().unstack(fill_value=0).reindex(columns=[0,1],fill_value=0);v=x[x.published_partition==fold].groupby(['pfam','label']).audit_cluster.nunique().unstack(fill_value=0).reindex(index=t.index,columns=[0,1],fill_value=0).fillna(0)
  for k in [3,5]:p.append(dict(fold=fold,threshold=k,train_qualifying=int(((t[0]>=k)&(t[1]>=k)).sum()),evaluation_both=int(((t[0]>=k)&(t[1]>=k)&(v[0]>0)&(v[1]>0)).sum())))
 save(pd.DataFrame(p),'published_PSI_family_feasibility.tsv')
 # Named distributions enriched with full biological clan descriptions.
 names=read(CUR/'presentation_assets/pfam37_clan_names.tsv').set_index('clan').description
 dist=read(OUT/'family_clan_prevalence.tsv.gz');mask=dist.level=='clan';dist.loc[mask,'description']=dist.loc[mask,'accession'].map(names);dist['class_rank']=dist.groupby(['policy','level','label']).N.rank(method='min',ascending=False).astype(int);save(dist,'family_clan_prevalence.tsv.gz');save(dist[(dist.class_rank<=20)|(dist.combined_rank<=20)],'top20_named_families_clans.tsv')
 print(pd.DataFrame(summary).to_string(index=False))
if __name__=='__main__':main()
