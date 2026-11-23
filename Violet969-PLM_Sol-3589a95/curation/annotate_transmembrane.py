#!/usr/bin/env python3
"""Resumable uniform DeepTMHMM2 batch runner, after the successful pilot.

CPU throughput was ~10–11 hours/cohort on this host. No full run is implicit in
the main annotation script. Checkpoints and temporary FASTAs stay outside git.
"""
import argparse,json,os,shutil,subprocess,time
from pathlib import Path
from annotate_biology import OUT,IDENTITY,load,save,sha,write_json,check_inputs

def run(out,work,executable,source,device,threads,batch_records,max_shards):
    d,seqs=load();check_inputs(out)
    pilot=json.loads((out/'transmembrane_pilot.json').read_text())
    assert pilot['exit_code']==0
    for name,h in pilot['checkpoint_hashes'].items():assert sha(source/'checkpoints'/name)==h
    work.mkdir(parents=True,exist_ok=True);destination=out/'transmembrane_shards';destination.mkdir(exist_ok=True)
    settings=dict(device=device,threads=threads,batch_size=1,source_commit=pilot['source_commit'],checkpoint_hashes=pilot['checkpoint_hashes'])
    completed=0
    for first in range(0,len(d),batch_records):
        group=d.iloc[first:first+batch_records];dest=destination/f'{first:05d}'
        fingerprint=dict(settings=settings,sequences=group[['raw_uid','sha256']].to_dict('records'))
        if (dest/'completed.json').exists():
            old=json.loads((dest/'completed.json').read_text());assert old['fingerprint']==fingerprint
            for name,h in old['outputs'].items():assert sha(dest/name)==h
            continue
        if max_shards is not None and completed>=max_shards:break
        tmp=work/f'{first:05d}';tmp.mkdir(exist_ok=True)
        fasta=tmp/'input.fasta'
        with fasta.open('w') as f:
            for u in group.raw_uid:f.write('>'+u+'\n'+seqs[u]+'\n')
        cmd=[str(executable),str(fasta),str(tmp/'results'),'--model-dir',str(source/'checkpoints'),
             '--device',device,'--threads',str(threads),'--batch-size','1','--simplify-io']
        start=time.monotonic()
        with (tmp/'run.log').open('w') as f:subprocess.run(cmd,stdout=f,stderr=subprocess.STDOUT,check=True)
        lines=(tmp/'results/predicted_topologies.3line').read_text().splitlines();assert len(lines)==3*len(group)
        observed=set()
        for h,s,top in zip(lines[::3],lines[1::3],lines[2::3]):
            uid=h[1:].split(' | ')[0];assert uid not in observed and s==seqs[uid] and len(s)==len(top);observed.add(uid)
        assert observed==set(group.raw_uid)
        dest.mkdir(exist_ok=True)
        for p in (tmp/'results').iterdir():
            if p.is_file():shutil.copyfile(p,dest/p.name)
        shutil.copyfile(tmp/'run.log',dest/'run.log')
        write_json(dest/'completed.json',dict(fingerprint=fingerprint,elapsed_seconds=time.monotonic()-start,
            outputs={p.name:sha(p) for p in dest.iterdir() if p.name!='completed.json'}))
        completed+=1;print(f'Completed {first+len(group)}/{len(d)}',flush=True)
    check_inputs(out)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--out',type=Path,default=OUT);p.add_argument('--work',type=Path,default=Path('/tmp/uesolds-tm-batch'))
    p.add_argument('--executable',type=Path,default=Path('/tmp/uesolds-tm-venv/bin/dtm2'));p.add_argument('--source',type=Path,default=Path('/tmp/uesolds-deeptmhmm2'))
    p.add_argument('--device',default='cpu',choices=['cpu','cuda']);p.add_argument('--threads',type=int,default=4)
    p.add_argument('--batch-records',type=int,default=128);p.add_argument('--max-shards',type=int)
    a=p.parse_args();run(a.out,a.work,a.executable,a.source,a.device,a.threads,a.batch_records,a.max_shards)
