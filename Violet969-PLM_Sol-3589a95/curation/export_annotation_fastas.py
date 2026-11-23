#!/usr/bin/env python3
"""Server upload packets outside the checkout, with sequence/hash/ID mappings.

All development records are exported for NetSurfP. SignalP settings use only
exact-match species consensus lineage; unresolved records get a separate file.
"""
import argparse,json,zipfile
from pathlib import Path
import pandas as pd
from annotate_biology import OUT,load,read,save,sha,write_json

def export(out,destination):
    d,seqs=load();t=read(out/'taxonomy.tsv')
    d=d.merge(t[['raw_uid','taxonomy_status','taxonomy_lineage']],on='raw_uid',validate='one_to_one')
    d['upload_id']=d.raw_uid.map(lambda x:'uesolds_'+x.replace(':','_'))
    def kingdom(r):
        if r.taxonomy_status!='assigned_exact_species_consensus':return 'unresolved'
        if ':Eukaryota' in r.taxonomy_lineage:return 'eukarya'
        if ':Bacteria' in r.taxonomy_lineage or ':Archaea' in r.taxonomy_lineage:return 'other'
        return 'unresolved'
    d['signalp_setting']=[kingdom(r) for r in d.itertuples()]
    destination.mkdir(parents=True,exist_ok=True);manifest=[];files=[]
    def fasta(g,path,tool,setting):
        path.parent.mkdir(exist_ok=True,parents=True)
        with path.open('w') as f:
            for r in g.itertuples():
                f.write('>'+r.upload_id+'\n'+seqs[r.raw_uid]+'\n')
                manifest.append(dict(tool=tool,file=str(path.relative_to(destination)),upload_id=r.upload_id,raw_uid=r.raw_uid,
                    sha256=r.sha256,length=r.length,split=r.split,setting=setting))
        files.append(dict(file=str(path.relative_to(destination)),sha256=sha(path),proteins=len(g),residues=int(g.length.sum()),tool=tool,setting=setting))
    fasta(d,destination/'all_development.fasta','archive','all_13880_do_not_submit_as_one_server_job')
    for split,g in d.groupby('split',sort=True):
        for length_class,h,size in [('le1022',g[(g.length>=10)&(g.length<=1022)],1000),('gt1022',g[(g.length>1022)&(g.length<=5000)],40)]:
            for first in range(0,len(h),size):
                fasta(h.iloc[first:first+size],destination/'netsurfp'/f'{split}_{length_class}_{first//size+1:02d}.fasta','NetSurfP3.0','full_sequence')
        unsupported=g[(g.length<10)|(g.length>5000)]
        if len(unsupported):fasta(unsupported,destination/'netsurfp_unsupported'/f'{split}.fasta','NetSurfP3.0','outside_documented_server_length_range')
        for setting,h in g.groupby('signalp_setting',sort=True):
            if setting=='unresolved':
                fasta(h,destination/'signalp_unresolved'/f'{split}.fasta','SignalP6.0','unresolved_do_not_guess')
            else:
                for first in range(0,len(h),100):
                    fasta(h.iloc[first:first+100],destination/'signalp'/f'{split}_{setting}_{first//100+1:02d}.fasta','SignalP6.0',setting)
    save(pd.DataFrame(manifest),destination/'upload_manifest.tsv')
    save(pd.DataFrame(manifest),out/'server_upload_manifest.tsv')
    write_json(out/'server_upload_files.json',files)
    counts=d.signalp_setting.value_counts().to_dict()
    text=f'''# UESolDS development annotation server uploads

Exactly {len(d):,} original proteins; 9,880 train + 4,000 validation. No independent test.
FASTA headers are neutral stable upload IDs. Sequences are unchanged and full length.
`upload_manifest.tsv` maps each upload ID to original record, split and SHA-256.
Repeated appearances across tools are intentional; do not deduplicate primary records.

## NetSurfP-3.0
Upload each file in `netsurfp/` separately at https://services.healthtech.dtu.dk/services/NetSurfP-3.0/ .
Use full sequences and default prediction settings. Download the complete per-residue CSV
output, preserve the original filename/job ID and submission date. Keep all records,
including any reported failures. Do not remove proteins or trim them to 1,022 residues.
The server documents a long-sequence memory issue: our >1,022-residue batches have
at most 40 proteins; the other batches have at most 1,000. These are operational
upload batches, not new biological cohorts or train/validation splits.

## SignalP-6.0
Upload each file in `signalp/` separately at https://services.healthtech.dtu.dk/services/SignalP-6.0/ .
Select **Eukarya** for filenames containing `eukarya`; select **Other** for `other`
(only reference-confirmed Bacteria/Archaea). Use **Slow** mode for region boundaries
and **Long output**. Batches have at most 100 proteins to reduce long-output timeouts.
Download full JSON, GFF3, prediction summary and input/processed-entry FASTA when
available. Retain original submitted sequences: a processed FASTA may have signal
peptides removed and must never replace the cohort sequences.

Settings supported by exact archive species consensus: {counts}.
Files in `signalp_unresolved/` are retained for accounting, **not ready for a
kingdom-specific submission**. Do not select a kingdom from centre, host or name.
Verified archived taxon does not establish a native/full-length construct.

## Audit of returned results
Save each job's raw outputs in a separate folder together with the input batch name,
server name/version displayed, date, chosen settings and job identifier. Remote model
checkpoint hashes are not exposed by these servers, so returned predictions must be
described as archived server outputs, not a byte-reproducible pinned local model.
Missing results are failed/missing annotations, never negatives. We will validate
IDs, returned sequence, one-based positions and full-length coverage before merging.

`all_development.fasta` is a convenient complete archive, **not one valid upload batch**.
This packet lives outside the git checkout so temporary FASTAs do not follow branches.
'''
    (destination/'README.md').write_text(text)
    archive=destination.with_suffix('.zip')
    with zipfile.ZipFile(archive,'w',compression=zipfile.ZIP_DEFLATED) as z:
        for p in sorted(destination.rglob('*')):
            if p.is_file():z.write(p,p.relative_to(destination))
    write_json(out/'server_upload_export.json',dict(destination=str(destination),archive=str(archive),archive_sha256=sha(archive),
        proteins=len(d),signalp_settings=counts,files=len(files),script_sha256=sha(Path(__file__))))
    print(json.dumps(dict(archive=str(archive),proteins=len(d),signalp_settings=counts,files=len(files)),indent=2))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--out',type=Path,default=OUT);p.add_argument('--destination',type=Path,default=Path('/tmp/uesolds-server-upload'));a=p.parse_args();export(a.out,a.destination)
