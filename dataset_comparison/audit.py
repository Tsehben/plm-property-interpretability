"""Experimental dataset audit. No model training or UESolDS test access."""
from pathlib import Path
import argparse,csv,gzip,hashlib,io,itertools,json,platform,subprocess,time,zipfile
import numpy as np
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
ROOT=Path(__file__).resolve().parent
CUR=ROOT.parent/'Violet969-PLM_Sol-3589a95/curation'
RAW=ROOT/'raw'; OUT=ROOT/'results'; WORK=Path('/tmp/uesolds-dataset-comparison')
AA='ACDEFGHIKLMNPQRSTVWY'
def sha(s):return hashlib.sha256(s.encode()).hexdigest()
def digest(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def read(p):return pd.read_csv(p,sep='\t',keep_default_na=False)
def save(d,n):d.to_csv(OUT/n,sep='\t',index=False)
def fasta(p):return pd.DataFrame([(r.id,str(r.seq)) for r in SeqIO.parse(p,'fasta')],columns=['sid','sequence'])
def prepare():
 OUT.mkdir(exist_ok=True);WORK.mkdir(exist_ok=True)
 rows=[];exclusions=[]
 def add(d,name):
  d=d.copy();d['dataset']=name;d['record_id']=[f'{name}:{i}' for i in range(len(d))];rows.append(d)
 c=pd.read_csv(RAW/'class.txt',sep='\t').rename(columns={'Accession':'sid'})
 a=pd.read_csv(RAW/'pET15_NESG.txt',sep='\t',header=None)[[0,4]];a.columns=['sid','dna'];a.dna='ATGGGCCATCACCATCACCATCACAGCCAT'+a.dna
 b=fasta(RAW/'pET21_NESG.fas').rename(columns={'sequence':'dna'})
 for v,d in [('pET15',a),('pET21',b)]:
  d=d.merge(c,on='sid').dropna();dups=d.duplicated(keep=False);exclusions.append(dict(stage=v+'_duplicate_join_records',N=int(dups.sum())));d=d[~dups]
  good=~d.dna.str.contains('[XZN*]',regex=True)&d.Class.isin([1,2]);exclusions.append(dict(stage=v+'_unmeasured_or_nonstandard_DNA',N=int((~good).sum())));d=d[good].copy()
  d['sequence']=d.dna.map(lambda s:str(Seq(s).translate()).removesuffix('*'));good=~d.sequence.str.contains(r'\*');exclusions.append(dict(stage=v+'_internal_stop',N=int((~good).sum())));d=d[good].copy();d['label']=d.Class-1;d['vector']=v
  add(d[['sid','sequence','label','vector','Class']],'PSI_original_'+v)
 psi=pd.concat(rows,ignore_index=True);psi['dataset']='PSI_original';psi.record_id=[f'PSI_original:{i}' for i in range(len(psi))];rows=[psi]
 assert len(psi)==12216 and psi.label.value_counts().to_dict()=={1:8238,0:3978}
 for filename,name,lab in [('PSI_Biology_solubility_trainset.csv','PSI_NetSolP','solublility|0=Insoluble|1=Soluble'),('NESG_testset.csv','NESG_NetSolP','solubility'),('Price_usability_trainset.csv','Price_usability','Usability|0=NotUsable|1=Usable')]:
  d=pd.read_csv(RAW/filename).rename(columns={'fasta':'sequence',lab:'label','Partition':'published_partition'});add(d,name)
 z=zipfile.ZipFile(RAW/'nesg.zip');scores=pd.read_csv(z.open('nesg.csv')).rename(columns={'id':'sid'});seq=fasta(io.StringIO(z.read('nesg.fasta').decode()));d=scores.merge(seq,on='sid',validate='one_to_one');assert len(d)==len(scores)==9703
 d['label']=np.nan;add(d,'NESG_raw')
 pd.crosstab(d.exp,d.sol,dropna=False).rename_axis('E').reset_index().to_csv(OUT/'NESG_E_by_S.tsv',sep='\t',index=False)
 u=read(CUR/'results/cohort_manifest.tsv');u=u[u.clean&u.cached].copy();assert len(u)==13880 and u.label.value_counts().to_dict()=={1:7996,0:5884}
 seqmap=read(CUR.parent/'esm2_probe/uesolds_pairs.tsv').set_index('key').sequence;u['sequence']=u.legacy_uid.map(seqmap);assert (u.sequence.map(sha)==u.sha256).all();u['sid']=u.original_id;add(u,'UESolDS_development')
 allr=pd.concat(rows,ignore_index=True);allr['sha256']=allr.sequence.map(sha);allr['length']=allr.sequence.str.len();allr['nonstandard']=allr.sequence.map(lambda s:''.join(sorted(set(s)-set(AA))));save(allr,'records.tsv.gz')
 unique=allr[['sha256','sequence','length','nonstandard']].drop_duplicates('sha256').sort_values('sha256');save(unique,'sequences.tsv.gz')
 with (WORK/'union.fasta').open('w') as f:
  for r in unique.itertuples():f.write(f'>{r.sha256}\n{r.sequence}\n')
 # Policy membership retains raw records; consistent duplicate policy picks deterministic ID only for counting.
 memberships=[];counts=[]
 def policy(name,q,parent,meaning='solubility'):
  q=q.copy();q['policy']=name;memberships.append(q[['policy','record_id','sha256','label']]);pos=int((q.label==1).sum());neg=int((q.label==0).sum());counts.append(dict(policy=name,total=len(q),soluble=pos if meaning=='solubility' else None,insoluble=neg if meaning=='solubility' else None,unclassified=len(q)-pos-neg,percent_soluble=100*pos/len(q) if len(q) and meaning=='solubility' else None,parent_N=parent,excluded=parent-len(q),outcome=meaning))
 for n in ['PSI_original','PSI_NetSolP','NESG_NetSolP','UESolDS_development']:
  q=allr[allr.dataset==n];policy(n,q,len(q))
 q=allr[allr.dataset=='Price_usability'];policy('Price_usability',q,len(q),'usability_NOT_solubility')
 q=allr[(allr.dataset=='UESolDS_development')&(allr.center=='DNASU_clone')];assert len(q)==460;policy('UESolDS_DNASU_proxy',q,len(q))
 d=allr[allr.dataset=='NESG_raw'].copy();d['sequence_repeated']=d.sha256.duplicated(False)
 g=d.groupby('sha256')[['exp','sol']].nunique();bad=set(g.index[g.max(axis=1)>1]);d['score_conflict']=d.sha256.isin(bad)
 save(d[d.sequence_repeated].drop(columns='sequence'),'NESG_repeated_records.tsv')
 policy('NESG_raw_scores',d,len(d),'E_and_S_ordinal_scores_NO_binary_label')
 uniq=d[~d.sequence_repeated];assert len(uniq)==9272;policy('NESG_unique_scores',uniq,len(d),'E_and_S_ordinal_scores_NO_binary_label')
 consistent=d[~d.score_conflict].sort_values('sid').drop_duplicates('sha256');policy('NESG_consistent_scores',consistent,len(d),'E_and_S_ordinal_scores_NO_binary_label')
 for prefix,pool in [('NESG_raw',d),('NESG_unique',uniq),('NESG_consistent',consistent)]:
  for tag,emin,smin,smaxneg in [('E45_S45_vs0',4,4,0),('E345_S45_vs0',3,4,0),('E45_S345_vs0',4,3,0),('E345_S345_vs0',3,3,0),('E45_S45_vs01',4,4,1)]:
   integer=pool.exp.isin(range(6))&pool.sol.isin(range(6));eligible=integer&(pool.exp>=emin)&((pool.sol>=smin)|(pool.sol<=smaxneg));q=pool[eligible].copy();q['label']=(q.sol>=smin).astype(int);policy(prefix+'_'+tag,q,len(pool))
  q=pool[pool.exp.isin([1,2,3,4,5])&pool.sol.isin(range(6))];policy(prefix+'_expressed_ordinal',q,len(pool),'S_ordinal_conditioned_on_E_positive')
 save(pd.concat(memberships,ignore_index=True),'policy_membership.tsv.gz');save(pd.DataFrame(counts),'dataset_counts.tsv');save(pd.DataFrame(exclusions),'PSI_reconstruction_exclusions.tsv')
 # Explicit construct mapping; sequence containment is descriptive, never edits source sequences.
 orig=allr[allr.dataset=='PSI_original'].set_index('sid');cur=allr[allr.dataset=='PSI_NetSolP'];links=[]
 for r in cur.itertuples():
  o=orig.loc[r.sid];assert isinstance(o,pd.Series)
  offset=o.sequence.find(r.sequence);links.append(dict(sid=r.sid,original_record=o.record_id,curated_record=r.record_id,original_sha256=o.sha256,curated_sha256=r.sha256,label_agrees=o.label==r.label,contained=offset>=0,N_removed=offset if offset>=0 else -1,C_removed=len(o.sequence)-offset-len(r.sequence) if offset>=0 else -1,vector=o.vector))
 save(pd.DataFrame(links),'PSI_construct_mapping.tsv')
 hashes={str(p.relative_to(ROOT.parent)):digest(p) for p in list(RAW.iterdir())+[CUR/'results/cohort_manifest.tsv',CUR.parent/'esm2_probe/uesolds_pairs.tsv',CUR/'pfam_annotation/domains.tsv',CUR/'pfam_annotation/searched_manifest.tsv'] if p.is_file()};(OUT/'input_checksums.json').write_text(json.dumps(hashes,indent=2)+'\n')
 print(pd.DataFrame(counts).to_string(index=False));print('union sequences',len(unique),flush=True)
def annotate():
 import pyhmmer
 assert pyhmmer.__version__=='0.11.2';ref=Path('/tmp/uesolds-pfam37');assert digest(ref/'Pfam-A.hmm.gz')=='66f57012fe5825c82f2e31447d7a8026367e1ca9e5e91b4b716c594c268db69b'
 import sys;sys.path.insert(0,str(CUR));from annotate_pfam import FIELDS
 old=read(CUR/'pfam_annotation/searched_manifest.tsv');oldhits=read(CUR/'pfam_annotation/domains.tsv');oldhits['sha256']=oldhits.raw_uid.map(old.set_index('raw_uid').sha256);oldhits=oldhits.drop_duplicates(['sha256','pfam','ali_from','ali_to']);oldhits['annotation_origin']='UESolDS_cached_Pfam37'
 u=read(OUT/'sequences.tsv.gz');new=u[~u.sha256.isin(old.sha256)];assert new.nonstandard.isin(['','X']).all(),new[new.nonstandard!=''][['sha256','nonstandard']]
 # X is passed unchanged to HMMER's documented unknown-residue alphabet; no amino-acid imputation.
 alphabet=pyhmmer.easel.Alphabet.amino();seqs=pyhmmer.easel.DigitalSequenceBlock(alphabet)
 for r in new.itertuples():seqs.append(pyhmmer.easel.TextSequence(name=r.sha256.encode(),sequence=r.sequence).digitize(alphabet))
 settings=dict(cpus=6,parallel='queries',bit_cutoffs='gathering',seed=42,Z=len(new),F1=.02,F2=.001,F3=.00001,bias_filter=True,null2=True);start=time.monotonic();n=0
 # Same implementation and domain inclusion criteria as completed UESolDS run; labels never passed.
 with gzip.open(OUT/'new_domains.tsv.gz','wt') as f,pyhmmer.plan7.HMMFile(ref/'Pfam-A.hmm.gz') as hmms:
  writer=csv.DictWriter(f,fieldnames=FIELDS,delimiter='\t');writer.writeheader()
  for hits in pyhmmer.hmmer.hmmsearch(hmms,seqs,**settings):
   m=hits.query;acc=m.accession.decode() if isinstance(m.accession,bytes) else m.accession;name=m.name.decode() if isinstance(m.name,bytes) else m.name;ga=m.cutoffs.gathering
   for hit in hits:
    if not hit.included:continue
    for dom in hit.domains:
     if not dom.included:continue
     a=dom.alignment;uid=hit.name.decode() if isinstance(hit.name,bytes) else hit.name
     assert hit.score>=ga[0] and dom.score>=ga[1]
     writer.writerow(dict(raw_uid=uid,pfam_versioned=acc,pfam=acc.split('.')[0],name=name,model_length=m.M,sequence_length=hit.length,sequence_score=hit.score,sequence_evalue=hit.evalue,domain_score=dom.score,conditional_evalue=dom.c_evalue,independent_evalue=dom.i_evalue,hmm_from=a.hmm_from,hmm_to=a.hmm_to,ali_from=a.target_from,ali_to=a.target_to,env_from=dom.env_from,env_to=dom.env_to,hmm_coverage=(a.hmm_to-a.hmm_from+1)/m.M,sequence_ga=ga[0],domain_ga=ga[1]))
   n+=1
   if n%1000==0:print('Pfam profiles',n,'seconds',round(time.monotonic()-start),flush=True)
 assert n==21979
 d=read(OUT/'new_domains.tsv.gz');d['sha256']=d.raw_uid;d['annotation_origin']='new_Pfam37';d=pd.concat([oldhits,d],ignore_index=True);d=d[d.sha256.isin(u.sha256)];save(d,'domains.tsv.gz')
 (OUT/'pfam_run.json').write_text(json.dumps(dict(settings=settings,release='37.0',pyhmmer=pyhmmer.__version__,new_sequences=len(new),reused_sequences=len(u)-len(new),profiles=n,elapsed_seconds=time.monotonic()-start,reference_sha256=digest(ref/'Pfam-A.hmm.gz'),domains_sha256=digest(OUT/'domains.tsv.gz')),indent=2)+'\n')
def homology():
 binary='/tmp/mmseqs/bin/mmseqs';base=['--min-seq-id','0.30','-c','0.80','--cov-mode','0','-s','7.5','-e','0.001','--alignment-mode','3','--threads','4']
 # Direct all-vs-all edges close leaks missed by representative clustering. Connected components thereafter.
 cmd=[binary,'easy-search',str(WORK/'union.fasta'),str(WORK/'union.fasta'),str(OUT/'homology_edges.tsv'),str(WORK/'search')]+base+['--max-seqs','10000','--format-output','query,target,fident,qcov,tcov,evalue,bits']
 with (OUT/'mmseqs.log').open('w') as f:subprocess.run(cmd,stdout=f,stderr=subprocess.STDOUT,check=True)
 with (OUT/'homology_edges.tsv').open('rb') as inp,gzip.GzipFile(filename=str(OUT/'homology_edges.tsv.gz'),mode='wb',mtime=0) as out:
  import shutil;shutil.copyfileobj(inp,out)
 (OUT/'homology_edges.tsv').unlink()
 (OUT/'homology_run.json').write_text(json.dumps(dict(command=cmd,version=subprocess.check_output([binary,'version'],text=True).strip(),binary_sha256=digest(binary),fasta_sha256=digest(WORK/'union.fasta')),indent=2)+'\n')
if __name__=='__main__':
 a=argparse.ArgumentParser();a.add_argument('stage',choices=['prepare','annotate','homology']);args=a.parse_args();globals()[args.stage]()
