import json,os,sys,collections
OUT=sys.argv[1]; CACHE=os.path.join(OUT,'uniparc_cache')
sample=json.load(open(os.path.join(OUT,'sample_100.json')))
UPDB={'UniProtKB/Swiss-Prot','UniProtKB/TrEMBL','UniProtKB/Swiss-Prot protein isoforms'}
out={}
st=collections.Counter()
for r in sample:
    fp=os.path.join(CACHE,r['id'].replace('/','_').replace(':','_')+'.json')
    d=json.load(open(fp)); st[d['status']]+=1
    e={'upi':None,'up_acc':[],'up_tax':set(),'all_tax':set(),'pdb':[],'status':d['status'],
       'up_active':[], 'n_xref':0}
    res=d.get('result')
    if res:
        e['upi']=res.get('accession')
        for x in res.get('dbReference',[]):
            e['n_xref']+=1
            props={p['type']:p['value'] for p in x.get('property',[])}
            tax=props.get('NCBI_taxonomy_id')
            if tax: e['all_tax'].add(tax)
            if x['type'] in UPDB:
                e['up_acc'].append(x['id'])
                if x.get('active')=='Y':
                    e['up_active'].append(x['id'])
                    if tax: e['up_tax'].add(tax)
            if x['type']=='PDB': e['pdb'].append(x['id'])
    e['up_tax']=sorted(e['up_tax']); e['all_tax']=sorted(e['all_tax'])
    e['up_acc']=sorted(set(e['up_acc'])); e['up_active']=sorted(set(e['up_active']))
    e['pdb']=sorted(set(e['pdb']))
    out[r['id']]=e
json.dump(out,open(os.path.join(OUT,'uniparc_parsed.json'),'w'),indent=1)
print("statuses:",dict(st))
m=[e for e in out.values() if e['upi']]
print(f"\nwith UniParc hit: {len(m)}/100")
print(f"  with >=1 active UniProtKB acc: {sum(1 for e in m if e['up_active'])}")
print(f"  with >=1 any UniProtKB acc:    {sum(1 for e in m if e['up_acc'])}")
print(f"  with >=1 PDB xref:             {sum(1 for e in m if e['pdb'])}")
tc=collections.Counter(len(e['up_tax']) for e in m)
print("  #distinct taxa among ACTIVE UniProtKB xrefs:",dict(sorted(tc.items())))
tc2=collections.Counter(len(e['all_tax']) for e in m)
print("  #distinct taxa among ALL xrefs:",dict(sorted(tc2.items())))
allt=set()
for e in m: allt|=set(e['all_tax'])
print("\ntotal distinct taxids to resolve:",len(allt))
json.dump(sorted(allt),open(os.path.join(OUT,'taxids.json'),'w'))
