#!/usr/bin/env python3
"""Uniform full-sequence Pfam37.0 search of the frozen 13,880 cached records.
No test path, model scores, annotation preselection, tag stripping or truncation.
"""
import argparse, csv, gzip, hashlib, json, platform, shutil, time
from pathlib import Path
import pandas as pd
import pyhmmer
ROOT=Path(__file__).resolve().parent
BASE=ROOT.parent
URL='https://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam37.0/'
FILES=['Pfam-A.hmm.gz','Pfam.version.gz','Pfam-A.clans.tsv.gz']
FIELDS=['raw_uid','pfam_versioned','pfam','name','model_length','sequence_length','sequence_score','sequence_evalue','domain_score','conditional_evalue','independent_evalue','hmm_from','hmm_to','ali_from','ali_to','env_from','env_to','hmm_coverage','sequence_ga','domain_ga']
def digest(p,kind='sha256'):
    h=hashlib.new(kind)
    with p.open('rb') as f:
        for chunk in iter(lambda:f.read(8*1024*1024),b''):h.update(chunk)
    return h.hexdigest()
def text(x):return x.decode() if isinstance(x,bytes) else x

def run(ref,out,cpus):
    assert pyhmmer.__version__=='0.11.2'
    out.mkdir(parents=True,exist_ok=True)
    expected={line.split()[-1].lstrip('*'):line.split()[0] for line in (ref/'md5_checksums').read_text().splitlines() if line.strip()}
    refs={}
    for name in FILES:
        assert digest(ref/name,'md5')==expected[name],f'Checksum failure: {name}'
        refs[name]={'url':URL+name,'md5':expected[name],'sha256':digest(ref/name)}
    version=gzip.open(ref/'Pfam.version.gz','rt').read();assert '37.0' in version
    (out/'Pfam.version.txt').write_text(version)
    shutil.copyfile(ref/'md5_checksums',out/'release_md5_checksums.txt')
    with gzip.open(ref/'Pfam-A.clans.tsv.gz','rb') as f:(out/'Pfam-A.clans.tsv').write_bytes(f.read())
    manifest_path=ROOT/'results/cohort_manifest.tsv';pairs_path=BASE/'esm2_probe/uesolds_pairs.tsv'
    d=pd.read_csv(manifest_path,sep='\t',keep_default_na=False)
    d=d[d.clean & d.cached].copy().sort_values('raw_uid')
    assert len(d)==13880 and d.raw_uid.is_unique
    assert d.split.value_counts().to_dict()=={'train':9880,'validation':4000}
    pairs=pd.read_csv(pairs_path,sep='\t',keep_default_na=False)
    seqmap=dict(pairs.iloc[:,:2].itertuples(index=False,name=None))
    alphabet=pyhmmer.easel.Alphabet.amino();seqs=pyhmmer.easel.DigitalSequenceBlock(alphabet)
    for r in d.itertuples():
        s=seqmap[r.legacy_uid];assert hashlib.sha256(s.encode()).hexdigest()==r.sha256
        assert len(s)==r.length
        seqs.append(pyhmmer.easel.TextSequence(name=r.raw_uid.encode(),sequence=s).digitize(alphabet))
    d.to_csv(out/'searched_manifest.tsv',sep='\t',index=False)
    settings=dict(Pfam_release='37.0',PyHMMER=pyhmmer.__version__,cpus=cpus,parallel='queries',bit_cutoffs='gathering',seed=42,Z=len(d),F1=0.02,F2=0.001,F3=0.00001,bias_filter=True,null2=True)
    inputs={str(p.relative_to(BASE)):digest(p) for p in [manifest_path,pairs_path]}
    fingerprint={'inputs':inputs,'references':refs,'settings':settings}
    (out/'search_settings.json').write_text(json.dumps(fingerprint,indent=2)+'\n')
    start=time.monotonic();n=0;domains=0
    tmp=out/'domains.partial.tsv';models=[]
    with tmp.open('w') as f, pyhmmer.plan7.HMMFile(ref/'Pfam-A.hmm.gz') as hmms:
        writer=csv.DictWriter(f,fieldnames=FIELDS,delimiter='\t');writer.writeheader()
        for hits in pyhmmer.hmmer.hmmsearch(hmms,seqs,cpus=cpus,parallel='queries',bit_cutoffs='gathering',seed=42,Z=len(d),F1=0.02,F2=0.001,F3=0.00001,bias_filter=True,null2=True):
            model=hits.query;acc=text(model.accession);name=text(model.name);ga=model.cutoffs.gathering
            assert ga is not None
            models.append(dict(pfam_versioned=acc,pfam=acc.split('.')[0],name=name,model_length=model.M,sequence_ga=ga[0],domain_ga=ga[1]))
            for hit in hits:
                if not hit.included:continue
                for dom in hit.domains:
                    if not dom.included:continue
                    a=dom.alignment
                    row=dict(raw_uid=text(hit.name),pfam_versioned=acc,pfam=acc.split('.')[0],name=name,model_length=model.M,sequence_length=hit.length,sequence_score=hit.score,sequence_evalue=hit.evalue,domain_score=dom.score,conditional_evalue=dom.c_evalue,independent_evalue=dom.i_evalue,hmm_from=a.hmm_from,hmm_to=a.hmm_to,ali_from=a.target_from,ali_to=a.target_to,env_from=dom.env_from,env_to=dom.env_to,hmm_coverage=(a.hmm_to-a.hmm_from+1)/model.M,sequence_ga=ga[0],domain_ga=ga[1])
                    assert hit.score>=ga[0] and dom.score>=ga[1]
                    writer.writerow(row);domains+=1
            n+=1
            if n%100==0:
                f.flush(); print(f'Profiles {n}; included domains {domains}; elapsed {time.monotonic()-start:.1f}s',flush=True)
    assert n==21979,(n,'unexpected Pfam37 profile count')
    pd.DataFrame(models).to_csv(out/'profile_thresholds.tsv',sep='\t',index=False)
    tmp.replace(out/'domains.tsv')
    assert all(digest(BASE/p)==h for p,h in inputs.items())
    result=dict(**fingerprint,n_sequences=len(d),n_profiles=n,n_domains=domains,elapsed_seconds=time.monotonic()-start,python=platform.python_version(),platform=platform.platform(),script_sha256=digest(Path(__file__)),domains_sha256=digest(out/'domains.tsv'),completed=True)
    (out/'search_run.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:result[k] for k in ['n_sequences','n_profiles','n_domains','elapsed_seconds','completed']},indent=2),flush=True)
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--reference',type=Path,required=True);ap.add_argument('--out',type=Path,default=ROOT/'pfam_annotation');ap.add_argument('--cpus',type=int,default=6);a=ap.parse_args();run(a.reference,a.out,a.cpus)
