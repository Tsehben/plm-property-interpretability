"""Descriptive diversity, overlap and partition feasibility; no prediction models."""
from audit import *
from collections import defaultdict,Counter

def main():
 r=read(OUT/'records.tsv.gz');r['label']=pd.to_numeric(r.label,errors='coerce');r['exp']=pd.to_numeric(r.exp,errors='coerce');r['sol']=pd.to_numeric(r.sol,errors='coerce')
 u=read(OUT/'sequences.tsv.gz');m=read(OUT/'policy_membership.tsv.gz');m['label']=pd.to_numeric(m.label,errors='coerce')
 edges=pd.read_csv(OUT/'homology_edges.tsv.gz',sep='\t',names=['query','target','identity','qcov','tcov','evalue','bits']);assert (edges.identity>=.30).all() and (edges.qcov>=.8).all() and (edges.tcov>=.8).all();assert edges.groupby('query').size().max()<10000
 parent={x:x for x in u.sha256}
 def root(x):
  while parent[x]!=x:parent[x]=parent[parent[x]];x=parent[x]
  return x
 for a,b in edges[['query','target']].itertuples(index=False,name=None):
  x,y=root(a),root(b)
  if x!=y:parent[max(x,y)]=min(x,y)
 clusters={s:root(s) for s in parent};u['audit_cluster']=u.sha256.map(clusters)
 # Label-independent prospective split for feasibility, NOT an untouched historical holdout.
 u['audit_fold']=u.audit_cluster.map(lambda s:int(sha('dataset-audit-v1:'+s)[:16],16)%5)
 save(u.drop(columns='sequence'),'sequence_groups.tsv.gz');m=m.merge(u[['sha256','audit_cluster','audit_fold']],on='sha256',validate='many_to_one');save(m,'policy_groups.tsv.gz')
 d=read(OUT/'domains.tsv.gz');clan=pd.read_csv(CUR/'pfam_annotation/Pfam-A.clans.tsv',sep='\t',header=None,names=['pfam','clan','clan_name','family_name','description'],keep_default_na=False)
 d=d.merge(clan,on='pfam',validate='many_to_one');assert ((d.ali_from>=1)&(d.ali_to<=d.sequence_length)).all()
 # Conservative uncertainty analysis: discard both families when different same-clan hits overlap.
 d['competing_overlap']=False
 for _,g in d[d.clan!=''].groupby(['sha256','clan']):
  ids=g.index.tolist()
  for i,j in itertools.combinations(ids,2):
   a,b=d.loc[i],d.loc[j]
   if a.pfam!=b.pfam and max(a.ali_from,b.ali_from)<=min(a.ali_to,b.ali_to):d.loc[[i,j],'competing_overlap']=True
 save(d.drop(columns=['raw_uid']),'domains_with_clans.tsv.gz')
 member=d[['sha256','pfam','clan','family_name','clan_name','description']].drop_duplicates(['sha256','pfam']);save(member,'family_memberships.tsv.gz')
 hit=set(member.sha256);u['pfam_annotated']=u.sha256.isin(hit);save(u.drop(columns='sequence'),'annotation_status.tsv.gz')
 sensitivity=[('all_GA',member),('no_competing_overlap',d[~d.competing_overlap][member.columns].drop_duplicates(['sha256','pfam'])),('HMM_coverage80',d[d.hmm_coverage>=.8][member.columns].drop_duplicates(['sha256','pfam']))]
 allsummary=[];allsupport=[];allclans=[];distribution=[];biochem=[];quality=[]
 kd=dict(zip(AA,[1.8,2.5,-3.5,-3.5,2.8,-.4,-3.2,4.5,-3.9,3.8,1.9,-3.5,-1.6,-3.5,-4.5,-.8,-.7,4.2,-.9,-1.3]))
 for a in AA:u['frac_'+a]=u.sequence.map(lambda s:s.count(a)/len(s) if set(s)<=set(AA) else np.nan)
 u['charge_proxy']=u.sequence.map(lambda s:(s.count('K')+s.count('R')-s.count('D')-s.count('E'))/len(s) if set(s)<=set(AA) else np.nan)
 u['hydropathy']=u.sequence.map(lambda s:np.mean([kd[a] for a in s]) if set(s)<=set(AA) else np.nan)
 u['terminal_his']=u.sequence.map(lambda s:'HHHHH' in s[:35] or 'HHHHH' in s[-35:])
 save(u.drop(columns='sequence'),'sequence_features.tsv.gz')
 for policy,q in m.groupby('policy',sort=False):
  print('Analyzing',policy,flush=True)
  isbinary=q.label.notna().all() and policy!='Price_usability';n=len(q);q=q.copy();q['label']=q.label.fillna(-1).astype(int)
  for filtername,mm in sensitivity:
   x=q.merge(mm,on='sha256',how='inner')
   f=x.groupby('pfam').agg(name=('family_name','first'),clan=('clan','first'),description=('description','first'),proteins=('record_id','nunique'),unique_sequences=('sha256','nunique'))
   def addcounts(df,value,prefix):
    tab=df.groupby(['pfam','label'])[value].nunique().unstack(fill_value=0).reindex(index=f.index,columns=[0,1],fill_value=0).fillna(0)
    for lab in [0,1]:f[prefix+str(lab)]=tab[lab].astype(int)
   binary=x[x.label.isin([0,1])]
   addcounts(binary,'record_id','n');f=f.rename(columns={'n0':'negative','n1':'positive'})
   addcounts(binary,'audit_cluster','groups_')
   labelcount=binary.groupby(['pfam','audit_cluster']).label.nunique().rename('cluster_labels').reset_index()
   ex=binary.merge(labelcount,on=['pfam','audit_cluster']);addcounts(ex[ex.cluster_labels==1],'audit_cluster','exclusive_groups_')
   f['shared_groups']=f.groups_0-f.exclusive_groups_0
   for split,df in [('train',binary[binary.audit_fold!=0]),('validation',binary[binary.audit_fold==0])]:
    addcounts(df,'record_id',split+'_n');addcounts(df,'audit_cluster',split+'_g')
   f=f.reset_index();f['policy']=policy;f['filter']=filtername;allsupport.extend(f.to_dict('records'))
   z=dict(policy=policy,filter=filtername,total=n,annotated=x.record_id.nunique(),annotation_percent=100*x.record_id.nunique()/n,families=x.pfam.nunique(),clans=x.loc[x.clan!='','clan'].nunique(),groups=q.audit_cluster.nunique(),unique_sequences=q.sha256.nunique(),multifamily_proteins=int((x.groupby('record_id').pfam.nunique()>1).sum()),mixed_families=int(((f.negative>0)&(f.positive>0)).sum()) if isbinary else np.nan)
   for k in [3,5]:
    adequate=(f.negative>=k)&(f.positive>=k);independent=(f.groups_0>=k)&(f.groups_1>=k);exclusive=(f.exclusive_groups_0>=k)&(f.exclusive_groups_1>=k)
    for label,mask in [('proteins',adequate),('groups',independent),('exclusive_groups',exclusive)]:
     z[f'families_{label}_ge{k}']=int(mask.sum()) if isbinary else np.nan
     z[f'proteins_in_{label}_ge{k}']=x.loc[x.pfam.isin(f.loc[mask,'pfam']),'record_id'].nunique() if isbinary else np.nan
    for criterion,trainmask in [('proteins',(f.train_n0>=k)&(f.train_n1>=k)),('groups',(f.train_g0>=k)&(f.train_g1>=k))]:
     z[f'train_{criterion}_ge{k}_eval_both']=int((trainmask&(f.validation_g0>=1)&(f.validation_g1>=1)).sum()) if isbinary else np.nan
     z[f'train_{criterion}_ge{k}_eval_ge{k}']=int((trainmask&(f.validation_g0>=k)&(f.validation_g1>=k)).sum()) if isbinary else np.nan
   allsummary.append(z)
   if filtername!='all_GA':continue
   # Each protein counted once per family/clan. Columns are prevalence, not a 100% partition.
   for kind,col,namecol in [('family','pfam','family_name'),('clan','clan','clan_name')]:
    y=x[x[col]!=''].drop_duplicates(['record_id',col]);totals=y.groupby(col).record_id.nunique().sort_values(ascending=False);ranks={v:i+1 for i,v in enumerate(totals.index)}
    for (cat,lab),g in y.groupby([col,'label']):
     den=int((q.label==lab).sum());distribution.append(dict(policy=policy,level=kind,accession=cat,name=g[namecol].iloc[0],description=g.description.iloc[0] if kind=='family' else '',label=lab,N=g.record_id.nunique(),denominator=den,percent=100*g.record_id.nunique()/den,combined_rank=ranks[cat],combined_N=totals[cat]))
  z=q.merge(u.drop(columns=['sequence','audit_cluster','audit_fold']),on='sha256');
  for lab,g in z.groupby('label'):
   quality.append(dict(policy=policy,label=lab,total=len(g),annotated=int(g.pfam_annotated.sum()),no_hit=int((~g.pfam_annotated).sum()),unknown_residues=int((g.nonstandard!='').sum()),beyond1022=int((g.length>1022).sum()),terminal_his=int(g.terminal_his.sum())))
   for feat in ['length','charge_proxy','hydropathy']+['frac_'+a for a in AA]:
    v=g[feat];biochem.append(dict(policy=policy,label=lab,feature=feat,N=int(v.notna().sum()),mean=v.mean(),q10=v.quantile(.1),median=v.median(),q90=v.quantile(.9)))
 save(pd.DataFrame(allsummary),'family_diversity.tsv');save(pd.DataFrame(allsupport),'family_support.tsv.gz');dist=pd.DataFrame(distribution);save(dist,'family_clan_prevalence.tsv.gz');save(dist[dist.combined_rank<=20],'top20_named_families_clans.tsv');save(pd.DataFrame(biochem),'biochemical_summary.tsv');save(pd.DataFrame(quality),'annotation_coverage_by_label.tsv')
 # Label-blind largest aligned-domain representative, for mutually exclusive composition tables.
 rep=d.assign(span=d.ali_to-d.ali_from+1).sort_values(['sha256','span','domain_score','pfam','ali_from'],ascending=[True,False,False,True,True]).drop_duplicates('sha256')
 reprcounts=[]
 for policy,q in m.groupby('policy',sort=False):
  q=q.copy();q['label']=q.label.fillna(-1).astype(int);x=q.merge(rep[['sha256','pfam','family_name','clan','clan_name']],on='sha256',how='left')
  for kind,col,namecol in [('family','pfam','family_name'),('clan','clan','clan_name')]:
   x[col]=x[col].fillna('NO_PFAM');x[namecol]=x[namecol].fillna('No Pfam assignment');x.loc[x[col]=='',col]='NO_CLAN';x.loc[x[col]=='NO_CLAN',namecol]='Annotated family without clan'
   for (lab,acc),g in x.groupby(['label',col]):reprcounts.append(dict(policy=policy,level=kind,label=lab,accession=acc,name=g[namecol].iloc[0],N=len(g),denominator=int((q.label==lab).sum()),percent=100*len(g)/(q.label==lab).sum()))
 save(pd.DataFrame(reprcounts),'representative_distributions.tsv')
 # Overlap: record counts in each side; sequences are not repeated experimental replicates.
 policies=['PSI_original','PSI_NetSolP','NESG_raw_scores','NESG_unique_E45_S45_vs0','NESG_NetSolP','UESolDS_development','UESolDS_DNASU_proxy'];overlap=[];conflicts=[];labelsets={}
 for pol in policies:
  q=m[m.policy==pol];labelsets[pol]=q.groupby('sha256').label.agg(lambda z:set(z.dropna().astype(int)))
 for a,b in itertools.combinations(policies,2):
  qa=m[m.policy==a];qb=m[m.policy==b];sa=set(qa.sha256);sb=set(qb.sha256);common=sa&sb
  z=dict(A=a,B=b,A_records=len(qa),B_records=len(qb),exact_sequences=len(common),A_exact_records=int(qa.sha256.isin(common).sum()),B_exact_records=int(qb.sha256.isin(common).sum()))
  conflictset={s for s in common if labelsets[a][s] and labelsets[b][s] and len(labelsets[a][s]|labelsets[b][s])>1};z['exact_conflicting_sequences']=len(conflictset)
  for s in sorted(conflictset):conflicts.append(dict(A=a,B=b,sha256=s,A_labels=','.join(map(str,sorted(labelsets[a][s]))),B_labels=','.join(map(str,sorted(labelsets[b][s])))))
  sharedgroups=set(qa.audit_cluster)&set(qb.audit_cluster);z['shared_components']=len(sharedgroups)
  cross=edges[edges['query'].isin(sa)&edges.target.isin(sb)]
  for cutoff in [.3,.9,.95]:
   hits=cross[cross.identity>=cutoff];tag=str(int(100*cutoff));z['A_direct_'+tag]=int(qa.sha256.isin(hits['query']).sum());z['B_direct_'+tag]=int(qb.sha256.isin(hits.target).sum())
   if cutoff>=.9:
    bad=0
    for x,y in hits[['query','target']].drop_duplicates().itertuples(index=False,name=None):
     la,lb=labelsets[a][x],labelsets[b][y]
     if la and lb and len(la|lb)>1:bad+=1
    z['discordant_pairs_'+tag]=bad
  overlap.append(z)
 save(pd.DataFrame(overlap),'dataset_overlap.tsv');save(pd.DataFrame(conflicts),'exact_conflicts.tsv')
 # Published PSI folds audit against direct local-MMseqs edges (not a rerun of ggsearch global identity).
 pub=r[r.dataset=='PSI_NetSolP'][['sha256','published_partition']].copy();pub['published_partition']=pd.to_numeric(pub.published_partition).astype(int);pub=pub.drop_duplicates()
 pe=edges.merge(pub.rename(columns={'sha256':'query','published_partition':'query_fold'}),on='query').merge(pub.rename(columns={'sha256':'target','published_partition':'target_fold'}),on='target');pe=pe[(pe.query_fold!=pe.target_fold)&(pe['query']<pe.target)];save(pe,'PSI_published_crossfold_edges.tsv')
 # Existing UESolDS complete train/validation manifest hashes; no other sequence files and never test.
 full=read(CUR/'results/cohort_manifest.tsv');assert set(full.split)<= {'train','validation'}
 fullsummary=[]
 for pol in policies[:5]:
  q=m[m.policy==pol];shared=set(q.sha256)&set(full.sha256);fullsummary.append(dict(policy=pol,candidate_records_with_exact_match=int(q.sha256.isin(shared).sum()),UESolDS_train_validation_records=int(full.sha256.isin(shared).sum()),unique_shared_sequences=len(shared)))
 save(pd.DataFrame(fullsummary),'full_training_validation_hash_overlap.tsv')
 # Direct validation of curated score mapping.
 raw=r[r.dataset=='NESG_raw'].set_index('sid');n=r[r.dataset=='NESG_NetSolP'];checks=[]
 for row in n.itertuples():
  v=raw.loc[row.sid];checks.append(dict(sid=row.sid,exact_sequence=row.sha256==v.sha256,E_agrees=row.exp==v.exp,S_agrees=row.sol==v.sol,label_rule_agrees=row.exp in [4,5] and row.sol in [0,4,5] and row.label==int(row.sol>=4)))
 save(pd.DataFrame(checks),'NetSolP_NESG_score_checks.tsv');assert pd.DataFrame(checks)[['exact_sequence','E_agrees','S_agrees']].all().all()
 # DNASU sequence-vs-identifier relation, without stripping sequences for equality tests.
 proxy=r[(r.dataset=='UESolDS_development')&(r.center=='DNASU_clone')];orig=r[r.dataset=='PSI_original'].set_index('sid');mapping=[]
 for row in proxy.itertuples():
  sid=row.sid.split('_')[0];z=dict(record_id=row.record_id,original_id=row.sid,matched_sid=sid,present_in_original=sid in orig.index)
  if sid in orig.index:
   o=orig.loc[sid];z.update(label_agrees=row.label==o.label,exact_sequence=row.sha256==o.sha256,cohort_sequence_in_original=row.sequence in o.sequence,original_length=len(o.sequence),cohort_length=len(row.sequence),vector=o.vector)
  mapping.append(z)
 save(pd.DataFrame(mapping),'DNASU_proxy_mapping.tsv')
 # Verify pET15 class meanings against raw outcome fields.
 status=pd.read_csv(RAW/'pET15_NESG.txt',sep='\t',header=None).rename(columns={0:'sid',5:'expression_status',6:'solubility_status',7:'purification_status'})
 selected=r[(r.dataset=='PSI_original')&(r.vector=='pET15')]
 status=status[status.sid.isin(selected.sid)]
 pp=selected.merge(status[['sid','expression_status','solubility_status','purification_status']].drop_duplicates(),on='sid',validate='one_to_one');save(pp.groupby(['label','expression_status','solubility_status','purification_status']).size().rename('N').reset_index(),'pET15_status_audit.tsv')
 run=dict(records=len(r),unique_sequences=len(u),groups=len(set(clusters.values())),edge_rows=len(edges),max_targets_per_query=int(edges.groupby('query').size().max()),pfam_domains=len(d),overlapping_domains=int(d.competing_overlap.sum()),published_crossfold_pairs=len(pe),fasta_sha256=digest(WORK/'union.fasta'),scripts={p.name:digest(p) for p in [ROOT/'audit.py',Path(__file__)]})
 (OUT/'analysis_run.json').write_text(json.dumps(run,indent=2)+'\n');print(pd.DataFrame(allsummary).query("filter=='all_GA'").to_string(index=False))
if __name__=='__main__':main()
