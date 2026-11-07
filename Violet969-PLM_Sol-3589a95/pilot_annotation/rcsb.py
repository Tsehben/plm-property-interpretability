import json,os,sys,time,urllib.request
OUT=sys.argv[1]
P=json.load(open(os.path.join(OUT,'uniparc_parsed.json')))
pdbs=sorted({p.lower() for e in P.values() for p in e['pdb']})
print("PDB entries to query:",len(pdbs))
fp=os.path.join(OUT,'rcsb.json'); db=json.load(open(fp)) if os.path.exists(fp) else {}
Q="""query($ids:[String!]!){entries(entry_ids:$ids){rcsb_id
 polymer_entities{rcsb_id
  rcsb_polymer_entity_annotation{type annotation_id name}
  rcsb_entity_source_organism{ncbi_scientific_name ncbi_taxonomy_id}
  rcsb_entity_host_organism{ncbi_scientific_name ncbi_taxonomy_id}
  rcsb_polymer_entity_container_identifiers{uniprot_ids}}}}"""
todo=[p for p in pdbs if p not in db]
for k in range(0,len(todo),20):
    chunk=todo[k:k+20]
    body=json.dumps({'query':Q,'variables':{'ids':[c.upper() for c in chunk]}}).encode()
    req=urllib.request.Request("https://data.rcsb.org/graphql",data=body,
        headers={'Content-Type':'application/json'},method='POST')
    try:
        with urllib.request.urlopen(req,timeout=90) as r: d=json.loads(r.read().decode())
    except Exception as e: print("ERR",e); time.sleep(3); continue
    for en in (d.get('data',{}).get('entries') or []):
        db[en['rcsb_id'].lower()]=en
    print(f"  {min(k+20,len(todo))}/{len(todo)}",flush=True); time.sleep(0.3)
json.dump(db,open(fp,'w'))
print("entries fetched:",len(db))
npf=sum(1 for e in db.values() for pe in e['polymer_entities']
        for a in (pe.get('rcsb_polymer_entity_annotation') or []) if a['type']=='Pfam')
print("Pfam annotations found across entries:",npf)
