#!/usr/bin/env python3
"""Track the authorized background run and refresh merged outputs on completion."""
import argparse,json,os,subprocess,sys,time
from pathlib import Path
from annotate_biology import ROOT,OUT,write_json

def main(out):
    job=json.loads((out/'transmembrane_job.json').read_text());pid=job['pid']
    while True:
        records=set()
        for p in (out/'transmembrane_shards').glob('*/completed.json'):
            try:records.update(x['raw_uid'] for x in json.loads(p.read_text())['fingerprint']['sequences'])
            except (OSError,json.JSONDecodeError):pass  # atomic completion publication may still be in flight
        try:os.kill(pid,0);alive=True
        except ProcessLookupError:alive=False
        state='completed' if len(records)==13880 else 'running' if alive else 'stopped_incomplete'
        status=dict(state=state,completed_records=len(records),total_records=13880,pid=pid,
                    updated_utc=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()),log=job['log'])
        write_json(out/'transmembrane_progress.json',status)
        if state=='completed':
            with (out/'transmembrane_refresh.log').open('w') as f:
                for command in [
                    [sys.executable,str(ROOT/'annotate_biology.py'),'build','--out',str(out)],
                    [sys.executable,str(ROOT/'summarize_biology.py'),'--out',str(out)],
                    [sys.executable,str(ROOT/'render_biological_report.py'),'--out',str(out)]]:
                    subprocess.run(command,stdout=f,stderr=subprocess.STDOUT,check=True)
            status['merged_outputs_refreshed']=True;write_json(out/'transmembrane_progress.json',status)
            return
        if not alive:return
        time.sleep(30)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--out',type=Path,default=OUT);a=p.parse_args();main(a.out)
