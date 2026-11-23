#!/usr/bin/env python3
"""Uniform exact-sequence UniParc lookup; immutable response snapshots, no guessing.

Does not use centre, expression host, labels, names or legacy taxonomic hints.
"""
import argparse, concurrent.futures, gzip, hashlib, json, time
from pathlib import Path
from urllib.parse import urlencode
import requests
import pandas as pd
from annotate_biology import OUT,load,read,save,sha,write_json,check_inputs
PINNED_RELEASE='2026_04'

def snapshot(url,folder):
    key=hashlib.sha256(url.encode()).hexdigest();path=folder/(key+'.json.gz')
    if path.exists():
        with gzip.open(path,'rt') as f:return json.load(f)
    response=None
    for attempt in range(3):
        try:
            r=requests.get(url,headers={'Accept':'application/json','User-Agent':'UESolDS-academic-annotation/1.0'},timeout=60)
            r.raise_for_status();data=r.json()
            release=r.headers.get('X-UniProt-Release')
            if release and release!=PINNED_RELEASE:raise RuntimeError('Reference release changed: '+release)
            response=dict(url=url,retrieved_utc=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()),status='ok',
                headers={k:v for k,v in r.headers.items() if k.lower() in ['x-uniprot-release','x-uniprot-release-date','x-total-results','link']},data=data)
            if 'next' in r.headers.get('Link',''):raise RuntimeError('Unexpected pagination: refusing incomplete no-match calls')
            break
        except Exception as e:
            response=dict(url=url,status='failed',error=type(e).__name__+':'+str(e))
            time.sleep(attempt+1)
    # Failed calls are evidence, but not reusable successful no-match responses.
    if response['status']=='ok':
        with path.open('wb') as f:
            with gzip.GzipFile(filename='',mode='wb',fileobj=f,mtime=0) as z:z.write(json.dumps(response,sort_keys=True).encode())
    return response

def run(out,pilot):
    d,seqs=load();check_inputs(out)
    if pilot:d=d[d.raw_uid.isin(read(out/'pilot_manifest.tsv').raw_uid)]
    elif not (out/'taxonomy_pilot.json').exists():raise RuntimeError('Run --pilot first')
    folder=out/'taxonomy_responses';folder.mkdir(exist_ok=True)
    ids=d.raw_uid.tolist();md5={u:hashlib.md5(seqs[u].encode()).hexdigest().upper() for u in ids}
    batches=[ids[i:i+40] for i in range(0,len(ids),40)]
    def request_batch(group):
        query=' OR '.join('checksum:'+md5[u] for u in group)
        url='https://rest.uniprot.org/uniparc/search?'+urlencode(dict(query='('+query+')',format='json',size=500))
        return group,snapshot(url,folder)
    matches={};errors={};responses=[]
    with concurrent.futures.ThreadPoolExecutor(max_workers=3) as pool:
        for i,(group,res) in enumerate(pool.map(request_batch,batches)):
            responses.append({k:v for k,v in res.items() if k!='data'})
            if res['status']!='ok':
                for u in group:errors[u]=res.get('error','request_failed')
                continue
            lookup={e['sequence']['value']:e for e in res['data']['results']}
            for u in group:
                e=lookup.get(seqs[u]);matches[u]=e
                if e is not None:assert e['sequence']['md5'].upper()==md5[u]
            if (i+1)%20==0 or pilot:print(f'taxonomy sequence batches {i+1}/{len(batches)}',flush=True)
    taxids=sorted({int(x['commonTaxonId']) for e in matches.values() if e for x in e.get('commonTaxons',[]) if x.get('commonTaxonId') is not None})
    taxa={};tax_responses=[]
    # Lineage and rank queries are also batched; no lineage inferred from strings.
    urls=['https://rest.uniprot.org/taxonomy/search?'+urlencode(dict(query='('+' OR '.join('id:'+str(x) for x in taxids[i:i+40])+')',format='json',size=500)) for i in range(0,len(taxids),40)]
    with concurrent.futures.ThreadPoolExecutor(max_workers=3) as pool:
        for res in pool.map(lambda url:snapshot(url,folder),urls):
            tax_responses.append({k:v for k,v in res.items() if k!='data'})
            if res['status']=='ok':
                for t in res['data']['results']:taxa[int(t['taxonId'])]=t
    rows=[]
    for r in d.itertuples():
        u=r.raw_uid;e=matches.get(u)
        row=dict(raw_uid=u,sha256=r.sha256,taxonomy_status='failed' if u in errors else 'no_exact_match',
                 reference_accession='',match_identity='',match_coverage='',organism='',taxonomy_id='',taxonomy_lineage='',
                 taxonomy_confidence='',taxonomy_candidates='',taxonomy_error=errors.get(u,''),reference_sequence_sha256='',
                 taxonomy_mapping_method='exact_full_sequence_only')
        if e:
            candidates=e.get('commonTaxons',[])
            row.update(reference_accession=e['uniParcId'],match_identity=1.,match_coverage=1.,
                reference_sequence_sha256=hashlib.sha256(e['sequence']['value'].encode()).hexdigest(),
                taxonomy_candidates=json.dumps(candidates,sort_keys=True),taxonomy_status='exact_match_ambiguous_taxonomy')
            if len(candidates)==1 and candidates[0].get('commonTaxonId') is not None:
                tid=int(candidates[0]['commonTaxonId']);t=taxa.get(tid)
                if t is None:row['taxonomy_status']='exact_match_lineage_unavailable'
                else:
                    lineage=t.get('lineage',[])
                    # A species-level common ancestor supports species consensus,
                    # but a genus-level common ancestor does not identify an organism.
                    species_rank=t.get('rank')=='species' or any(x.get('rank')=='species' for x in lineage)
                    artificial=any(x.get('scientificName','').lower() in ['synthetic construct','artificial sequences'] for x in [t,*lineage])
                    row['taxonomy_lineage']=';'.join(f"{x.get('taxonId')}:{x.get('scientificName')}" for x in [*lineage,t])
                    if artificial:row['taxonomy_status']='exact_match_engineered_taxon'
                    elif species_rank:
                        row.update(taxonomy_status='assigned_exact_species_consensus',organism=t['scientificName'],taxonomy_id=tid,
                                   taxonomy_confidence='exact_archived_sequence_species_consensus_not_native_construct_verification')
        rows.append(row)
    result=pd.DataFrame(rows);save(result,out/('taxonomy_pilot.tsv' if pilot else 'taxonomy.tsv'))
    write_json(out/('taxonomy_pilot.json' if pilot else 'taxonomy_run.json'),dict(records=len(d),statuses=result.taxonomy_status.value_counts().to_dict(),
        passed=not errors,reference_release=PINNED_RELEASE,sequence_requests=responses,lineage_requests=tax_responses,
        script_sha256=sha(Path(__file__)),response_files={p.name:sha(p) for p in sorted(folder.glob('*.json.gz'))}))
    check_inputs(out);print(result.taxonomy_status.value_counts().to_string(),flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--out',type=Path,default=OUT);p.add_argument('--pilot',action='store_true');a=p.parse_args();run(a.out,a.pilot)
