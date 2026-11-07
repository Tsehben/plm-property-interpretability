import json,os,sys,time,urllib.request
OUT=sys.argv[1]
ids=json.load(open(os.path.join(OUT,'taxids.json')))
fp=os.path.join(OUT,'taxnames.json')
db=json.load(open(fp)) if os.path.exists(fp) else {}
todo=[i for i in ids if i not in db]
for k in range(0,len(todo),40):
    chunk=todo[k:k+40]
    q='+OR+'.join(f'tax_id:{c}' for c in chunk)
    url=f"https://rest.uniprot.org/taxonomy/search?query=({q})&size=100&fields=id,scientific_name,rank,lineage"
    try:
        with urllib.request.urlopen(url,timeout=90) as r: d=json.loads(r.read().decode())
    except Exception as e:
        print("ERR",k,e); time.sleep(3); continue
    for t in d.get('results',[]):
        sp=None
        for anc in (t.get('lineage') or []):
            if anc.get('rank')=='species': sp=(str(anc['taxonId']),anc['scientificName'])
        if t.get('rank')=='species': sp=(str(t['taxonId']),t['scientificName'])
        db[str(t['taxonId'])]={'name':t.get('scientificName'),'rank':t.get('rank'),
                               'species_id':sp[0] if sp else str(t['taxonId']),
                               'species_name':sp[1] if sp else t.get('scientificName')}
    print(f"resolved {min(k+40,len(todo))}/{len(todo)}",flush=True); time.sleep(0.3)
json.dump(db,open(fp,'w'),indent=0)
print("missing:",[i for i in ids if i not in db][:20],"total missing",len([i for i in ids if i not in db]))
