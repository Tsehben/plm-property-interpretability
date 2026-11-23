#!/usr/bin/env python3
"""Export only the frozen pilot to temporary storage for DeepTMHMM2 feasibility.

The predictor uses its own ESM final-layer inference; no existing ESM/SAE probe
is trained or altered. Model/runtime provenance is captured after the run.
"""
import argparse,json,subprocess,time
from pathlib import Path
from annotate_biology import OUT,load,read,save,sha,write_json

def run(out,work,python,source):
    d,s=load();pilot=read(out/'pilot_manifest.tsv')
    work.mkdir(exist_ok=True,parents=True)
    fasta=work/'pilot.fasta'
    with fasta.open('w') as f:
        for u in pilot.raw_uid:f.write('>'+u+'\n'+s[u]+'\n')
    cmd=[str(python.parent/'dtm2'),str(fasta),str(work/'results'),'--model-dir',str(source/'checkpoints'),
         '--device','cpu','--threads','4','--batch-size','1','--simplify-io']
    start=time.monotonic()
    with (work/'run.log').open('w') as f:r=subprocess.run(cmd,stdout=f,stderr=subprocess.STDOUT)
    model_files=list((source/'checkpoints').glob('*.ckpt'))
    backbone=work.parent/'uesolds-tm-models/hub/checkpoints/esm2_t33_650M_UR50D.pt'
    result=dict(command=cmd,exit_code=r.returncode,elapsed_seconds=time.monotonic()-start,
         source_commit='526f6c60f93c0868c7315413da8e669d12373cc3',pilot_records=len(pilot),
         checkpoint_hashes={p.name:sha(p) for p in model_files},backbone_sha256=sha(backbone) if backbone.exists() else '',
         log=(work/'run.log').read_text(),script_sha256=sha(Path(__file__)))
    write_json(out/'transmembrane_pilot.json',result)
    print(json.dumps({k:v for k,v in result.items() if k not in ['log','checkpoint_hashes']},indent=2))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--out',type=Path,default=OUT);p.add_argument('--work',type=Path,default=Path('/tmp/uesolds-tm-pilot'))
    p.add_argument('--python',type=Path,default=Path('/tmp/uesolds-tm-venv/bin/python'));p.add_argument('--source',type=Path,default=Path('/tmp/uesolds-deeptmhmm2'))
    a=p.parse_args();run(a.out,a.work,a.python,a.source)
