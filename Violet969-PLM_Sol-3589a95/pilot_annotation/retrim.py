import json,os,sys,re,time,urllib.request,urllib.error
OUT=sys.argv[1]; CACHE=os.path.join(OUT,'uniparc_trim_cache'); os.makedirs(CACHE,exist_ok=True)
S={r['id']:r for r in json.load(open(os.path.join(OUT,'sample_100.json')))}
P=json.load(open(os.path.join(OUT,'uniparc_parsed.json')))
URL="https://www.ebi.ac.uk/proteins/api/uniparc/sequence?rfActive=true"
def post(seq):
    req=urllib.request.Request(URL,data=seq.encode(),
        headers={'Accept':'application/json','Content-Type':'text/plain'},method='POST')
    with urllib.request.urlopen(req,timeout=90) as r: return json.loads(r.read().decode())
TEV=re.compile(r'ENLYFQ[GS]?')
def cands(s):
    out=[]
    def add(x,why):
        if x and len(x)>=30 and (x,why) not in out: out.append((x,why))
    # trailing vector/His
    t=s
    t=re.sub(r'(LE)?H{5,}(\*)?$','',t); add(t,'trim_C_His')
    t2=re.sub(r'LE$','',t); add(t2,'trim_C_His+LE')
    # leading His / vector
    for base,lbl in [(t,'trim_C_His'),(t2,'trim_C_His+LE'),(s,'orig')]:
        h=re.sub(r'^M?[A-Z]{0,12}?H{5,}[A-Z]{0,25}?(?=[A-Z]{30})','',base,count=1)
        add(h,lbl+'+trim_N_His')
        m=TEV.search(base[:60])
        if m: add(base[m.end():],lbl+'+trim_N_TEV')
        if base.startswith('M'): add(base[1:],lbl+'+trim_Met')
    return out[:6]
res={}
nm=[i for i,e in P.items() if not e['upi']]
for n,i in enumerate(nm):
    fp=os.path.join(CACHE,i.replace('/','_').replace(':','_')+'.json')
    if os.path.exists(fp): res[i]=json.load(open(fp)); continue
    rec={'id':i,'hit':None}
    for c,why in cands(S[i]['seq']):
        try:
            r=post(c); rec['hit']={'why':why,'trimmed_len':len(c),'result':r}; break
        except urllib.error.HTTPError as e:
            if e.code!=404: time.sleep(2)
        except Exception: time.sleep(2)
        time.sleep(0.3)
    json.dump(rec,open(fp,'w')); res[i]=rec
    print(f"{n+1:3d}/{len(nm)} {i:<26} {'HIT '+rec['hit']['why'] if rec['hit'] else 'still_no_match'}",flush=True)
print("\nrecovered:",sum(1 for r in res.values() if r['hit']),"/",len(nm))
