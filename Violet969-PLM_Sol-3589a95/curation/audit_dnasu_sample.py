#!/usr/bin/env python3
"""Parse the fixed 12-record public DNASU page sample; no label rewriting.
Run with --pages /tmp/uesolds-source-audit after downloading clone URLs in output table.
"""
from pathlib import Path
import argparse,re,html,hashlib,json,zipfile
import pandas as pd
ROOT=Path(__file__).resolve().parent;OUT=ROOT/'source_audit_results'
def plain(s):return re.sub(r'\s+',' ',html.unescape(re.sub('<[^>]+>',' ',s))).strip()
def main(pages):
 d=pd.read_csv(ROOT/'pfam_annotation/annotated_cohort.tsv',sep='\t',keep_default_na=False);sample=[]
 for _,g in d[d.center=='DNASU_clone'].groupby(['split','label']):
  g=g.sort_values(['length','raw_uid']);sample.append(g.iloc[[0,len(g)//2,len(g)-1]])
 s=pd.concat(sample);s['clone_id']=s.original_id.str.extract(r'CD0*(\d+)',expand=False)
 s.to_csv(OUT/'dnasu_recovery_sample.tsv',sep='\t',index=False)
 seq=dict(pd.read_csv(ROOT.parent/'esm2_probe/uesolds_pairs.tsv',sep='\t').itertuples(index=False,name=None));rows=[]
 for r in s.itertuples():
  p=pages/f'dnasu_{r.clone_id}.html';base=dict(raw_uid=r.raw_uid,original_id=r.original_id,split=r.split,label=r.label,
    clone_id=r.clone_id,url=f'https://dnasu.org/DNASU/GetCloneDetail.do?cloneid={r.clone_id}',cohort_sha256=r.sha256)
  if pages.is_file():
   with zipfile.ZipFile(pages) as z:h=z.read(p.name).decode()
  elif p.exists():h=p.read_text()
  else:rows.append(dict(**base,status='page_unavailable'));continue
  prot=re.search(r'id="protseq".*?<pre>(.*?)</pre>',h,re.S);txt=plain(h)
  protein=''.join(plain(prot[1]).split()) if prot else ''
  outcome=re.search(r'SolubleProtein:\s*(\S+)',txt);expression=re.search(r'ProteinExpressed:\s*(\S+)',txt)
  source=re.search(r'Source:\s*(.*?)\s*Description:',txt);vector=re.search(r'Vector:\s*(.*?)\s*Additional Annotations',txt)
  original=seq[r.legacy_uid];assert hashlib.sha256(original.encode()).hexdigest()==r.sha256
  agree=bool(outcome and outcome[1]=={0:'Tested_Not_Soluble',1:'Protein_Soluble'}[r.label])
  ix=original.find(protein) if protein else -1
  rows.append(dict(**base,status='retrieved',page_sha256=hashlib.sha256(h.encode()).hexdigest(),
    clone_name_present=r.original_id.split('_')[0] in txt,source=source[1] if source else '',vector=vector[1] if vector else '',
    solubility_status=outcome[1] if outcome else '',expression_status=expression[1] if expression else '',published_label_agrees=agree,
    listed_protein_length=len(protein),cohort_length=len(original),exact_listed_protein_match=protein==original,
    listed_protein_sha256=hashlib.sha256(protein.encode()).hexdigest() if protein else '',
    listed_protein_contained=ix>=0,cohort_extra_prefix=original[:ix] if ix>=0 else '',
    cohort_extra_suffix=original[ix+len(protein):] if ix>=0 else '',
    limitation='Listed protein is translated insert; tag-bearing assayed construct and original assay conditions still require reconciliation'))
 pd.DataFrame(rows).to_csv(OUT/'dnasu_sample_matches.tsv',sep='\t',index=False)
 print(pd.DataFrame(rows)[['raw_uid','label','solubility_status','published_label_agrees','exact_listed_protein_match','listed_protein_contained','cohort_extra_prefix','cohort_extra_suffix']].to_string(index=False))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--pages',type=Path,required=True);main(p.parse_args().pages)
