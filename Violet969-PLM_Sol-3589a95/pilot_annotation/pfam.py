import json,os,sys,time,urllib.request
OUT=sys.argv[1]
P=json.load(open(os.path.join(OUT,'uniparc_parsed.json')))
accs=sorted({a for e in P.values() for a in e['up_acc']})
print("accessions to fetch:",len(accs))
fp=os.path.join(OUT,'uniprot_entries.json')
db=json.load(open(fp)) if os.path.exists(fp) else {}
todo=[a for a in accs if a not in db]
for k in range(0,len(todo),25):
    chunk=todo[k:k+25]
    q='+OR+'.join(f'accession:{c}' for c in chunk)
    url=f"https://rest.uniprot.org/uniprotkb/search?query=({q})&size=100&fields=accession,organism_name,organism_id,xref_pfam,protein_name,reviewed,length"
    try:
        with urllib.request.urlopen(url,timeout=90) as r: d=json.loads(r.read().decode())
    except Exception as e:
        print("ERR",k,e); time.sleep(4); continue
    for t in d.get('results',[]):
        pf=[(x['id'],{p['key']:p['value'] for p in x.get('properties',[])}.get('EntryName'))
            for x in t.get('uniProtKBCrossReferences',[]) if x['database']=='Pfam']
        db[t['primaryAccession']]={'organism':t['organism']['scientificName'],
            'taxid':str(t['organism']['taxonId']),'pfam':pf,
            'name':(t.get('proteinDescription',{}).get('recommendedName',{}) or {}).get('fullName',{}).get('value'),
            'reviewed':t.get('entryType','').startswith('UniProtKB reviewed'),
            'length':t.get('sequence',{}).get('length')}
    print(f"  {min(k+25,len(todo))}/{len(todo)}",flush=True); time.sleep(0.3)
json.dump(db,open(fp,'w'),indent=0)
miss=[a for a in accs if a not in db]
print("fetched:",len(db),"| unresolved accessions (obsolete/deleted):",len(miss),miss[:10])
n_pf=sum(1 for a in db if db[a]['pfam'])
print("accessions with >=1 Pfam:",n_pf)
