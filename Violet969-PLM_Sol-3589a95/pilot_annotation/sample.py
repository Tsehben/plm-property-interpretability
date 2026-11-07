import re, random, collections, json, os, sys
BASE=sys.argv[1]; OUT=sys.argv[2]
SEED=42; N=100; EXCLUDE={'PDB_230301_28958'}
def parse(p):
    recs=[];h=None;s=[]
    for l in open(p):
        l=l.strip()
        if l.startswith('>'):
            if h is not None: recs.append((h,''.join(s)))
            h=l[1:];s=[]
        else: s.append(l)
    if h is not None: recs.append((h,''.join(s)))
    return recs
def idclass(i):
    if re.match(r'^PDB_230301_\d+$',i): return 'PDB_serial'
    if re.match(r'^(NP|XP|YP|WP)_\d+',i): return 'RefSeq'
    if re.match(r'^[A-Za-z]{2}CD\d+_pET',i): return 'DNASU_clone'
    if i.startswith('TB-Rv'): return 'TB_Rv'
    if re.match(r'^(JCSG|MCSG|NESG|NYSGRC|SSGCID|SGPP|SECSG|EFI|BSGI|BSGC|MPP|CESG|SGC|APC)',i): return 'SG_target'
    return 'other'
recs=parse(os.path.join(BASE,'embedding_dataset','Train_dataset.fasta'))
rows=[]
for h,s in recs:
    p=h.split(); pid=p[1] if len(p)==3 else p[0]
    if pid in EXCLUDE: continue
    lab=p[-1].split('-')[-1]
    rows.append({'id':pid,'label':lab,'seq':s,'len':len(s),'cls':idclass(pid),
                 'center':pid.split('-')[0] if '-' in pid else ''})
print("pool after exclusion:",len(rows),"(dropped",len(recs)-len(rows),")")
strata=collections.defaultdict(list)
for r in rows: strata[(r['cls'],r['label'])].append(r)
total=len(rows)
# largest-remainder proportional allocation
raw={k:len(v)/total*N for k,v in strata.items()}
alloc={k:int(v) for k,v in raw.items()}
rem=N-sum(alloc.values())
for k,_ in sorted(raw.items(), key=lambda kv:(-(kv[1]-int(kv[1])), kv[0]))[:rem]: alloc[k]+=1
rng=random.Random(SEED)
sample=[]
for k in sorted(strata):
    pool=sorted(strata[k], key=lambda r:r['id'])
    n=min(alloc[k],len(pool))
    if n: sample+=rng.sample(pool,n)
sample.sort(key=lambda r:(r['cls'],r['label'],r['id']))
print(f"\n{'id_class':<14}{'label':<7}{'pool':>8}{'prop%':>8}{'sampled':>9}")
for k in sorted(strata):
    print(f"{k[0]:<14}{k[1]:<7}{len(strata[k]):>8}{len(strata[k])/total*100:>7.2f}%{alloc[k]:>9}")
print(f"{'TOTAL':<21}{total:>8}{100.0:>7.2f}%{sum(alloc.values()):>9}  | sampled={len(sample)}")
with open(os.path.join(OUT,'sample_100.json'),'w') as f: json.dump(sample,f,indent=1)
with open(os.path.join(OUT,'sample_100.fasta'),'w') as f:
    for r in sample: f.write(f">{r['id']} A-{r['label']}\n{r['seq']}\n")
print("\nwrote sample_100.json / sample_100.fasta")
