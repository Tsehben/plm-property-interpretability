#!/usr/bin/env python3
"""Recover fixed development sample from the public 2017 TargetTrack archive.
No UESolDS test reads. Does not infer outcomes or modify eligibility.
First run source_validity_audit.py --stage diversity, then:
python recover_targettrack_sample.py --xml /tmp/uesolds-source-audit/tt.xml.gz
"""
from pathlib import Path
import argparse,collections,gzip,hashlib,json,xml.etree.ElementTree as ET
import pandas as pd
ROOT=Path(__file__).resolve().parent;OUT=ROOT/'source_audit_results'
def text(e,path):return (e.findtext(path) or '').strip()
def sequence(e):return ''.join(text(e,'oneLetterCode').split()).upper()
def digest(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(1048576),b''):h.update(b)
 return h.hexdigest()
def main(xml):
 sample=pd.read_csv(OUT/'recovery_sample.tsv',sep='\t',keep_default_na=False)
 pairs=pd.read_csv(ROOT.parent/'esm2_probe/uesolds_pairs.tsv',sep='\t',keep_default_na=False)
 seq=dict(pairs.itertuples(index=False,name=None));byid={k:g for k,g in sample.groupby('original_id')}
 evidence=OUT/'targettrack_sample';evidence.mkdir(exist_ok=True)
 protocols={};counts=collections.Counter();found={};n=0
 with gzip.open(xml,'rb') as f:
  for event,e in ET.iterparse(f,events=['end']):
   if e.tag=='protocol':
    pid=e.attrib.get('id','');center=pid.split('-')[0];kind=text(e,'protocolType');counts[(center,kind)]+=1
    if center in ['JCSG','MCSG']:
     protocols[pid]=dict(protocol_id=pid,name=text(e,'protocolName'),type=kind,description=text(e,'protocolDescription'),text=text(e,'protocolText'))
    e.clear()
   elif e.tag=='target':
    tid=e.attrib.get('id','');n+=1
    if tid in byid:
     found[tid]=ET.tostring(e,encoding='unicode')
     (evidence/(tid+'.xml')).write_text(found[tid])
    e.clear()
 rows=[];trialrows=[];needed=set()
 for r in sample.itertuples():
  base=dict(raw_uid=r.raw_uid,original_id=r.original_id,center=r.center,label=r.label,split=r.split,cohort_sha256=r.sha256)
  s=seq[r.legacy_uid];assert hashlib.sha256(s.encode()).hexdigest()==r.sha256
  if r.original_id not in found:
   rows.append(dict(**base,recovered=False,interpretation='Target ID not found; outcome unknown'));continue
  e=ET.fromstring(found[r.original_id]);targets=e.findall('./targetSequenceList/targetSequence');trials=e.findall('./trialList/trial')
  exact_targets=[x.attrib.get('id','') for x in targets if sequence(x)==s]
  exact_trials=[];status=[]
  for trial in trials:
   sequences=trial.findall('./trialSequenceList/trialSequence');match=any(sequence(x)==s for x in sequences)
   refs=[]
   for ref in trial.findall('./trialProtocolList/protocolRef'):
    pid=text(ref,'protocolId');full=pid if pid in protocols else r.center+'-'+pid
    refs.append(full);needed.add(full)
   histories=[text(h,'status') for h in trial.findall('./statusHistoryList/statusHistory')]
   # Store status text literally; no negative class is inferred from pipeline cessation.
   trialrows.append(dict(**base,trial_id=trial.attrib.get('id',''),exact_sequence_match=match,
                         trial_sequence_hashes=';'.join(hashlib.sha256(sequence(x).encode()).hexdigest() for x in sequences),
                         current_status=text(trial,'status'),history=';'.join(histories),protocol_ids=';'.join(refs),
                         remarks=' | '.join((x.text or '').strip() for x in trial.iter() if x.tag in ['remark','stopDetails','failureReason','experimentDescription','result','experimentalResult'])))
   if match:exact_trials.append(trial.attrib.get('id',''));status.extend(histories+[text(trial,'status')])
  rows.append(dict(**base,recovered=True,target_status=text(e,'status'),target_sequence_count=len(targets),
                   exact_target_sequence_ids=';'.join(exact_targets),trial_count=len(trials),exact_trial_ids=';'.join(exact_trials),
                   exact_trial_statuses=';'.join(sorted(set(status))),target_remark=text(e,'remark'),
                   interpretation='Exact trial sequence linkage; outcome still requires experiment/protocol adjudication' if exact_trials else 'Target identity/sequence alone does not establish outcome for released construct'))
 pd.DataFrame(rows).to_csv(OUT/'targettrack_sample_matches.tsv',sep='\t',index=False)
 pd.DataFrame(trialrows).to_csv(OUT/'targettrack_sample_trials.tsv',sep='\t',index=False)
 pd.DataFrame([dict(center=c,protocol_type=k,n=n) for (c,k),n in sorted(counts.items())]).to_csv(OUT/'targettrack_protocol_counts.tsv',sep='\t',index=False)
 (evidence/'referenced_protocols.json').write_text(json.dumps({k:protocols.get(k,{'unresolved_reference':True}) for k in sorted(needed)},indent=2)+'\n')
 run=dict(archive_url='https://zenodo.org/records/821654',xml_sha256=digest(xml),sample_sha256=digest(OUT/'recovery_sample.tsv'),
          script_sha256=digest(Path(__file__)),targets_scanned=n,sample_records=len(sample),target_ids_found=len(found),referenced_protocols=len(needed))
 (OUT/'targettrack_recovery_run.json').write_text(json.dumps(run,indent=2)+'\n')
 print(json.dumps(run,indent=2))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--xml',type=Path,required=True);main(p.parse_args().xml)
