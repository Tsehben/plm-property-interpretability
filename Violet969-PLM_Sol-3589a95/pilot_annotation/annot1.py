import json,os,sys,time,urllib.request,urllib.error
OUT=sys.argv[1]; CACHE=os.path.join(OUT,'uniparc_cache')
os.makedirs(CACHE,exist_ok=True)
sample=json.load(open(os.path.join(OUT,'sample_100.json')))
URL="https://www.ebi.ac.uk/proteins/api/uniparc/sequence?rfActive=true"
def post(seq):
    req=urllib.request.Request(URL,data=seq.encode(),
        headers={'Accept':'application/json','Content-Type':'text/plain'},method='POST')
    with urllib.request.urlopen(req,timeout=90) as r: return json.loads(r.read().decode())
for i,r in enumerate(sample):
    fp=os.path.join(CACHE,r['id'].replace('/','_').replace(':','_')+'.json')
    if os.path.exists(fp): continue
    rec={'id':r['id']}
    for attempt in range(3):
        try:
            rec['result']=post(r['seq']); rec['status']='ok'; break
        except urllib.error.HTTPError as e:
            if e.code==404: rec['status']='no_match'; rec['result']=None; break
            rec['status']=f'http_{e.code}'; time.sleep(3)
        except Exception as e:
            rec['status']='error:'+type(e).__name__; time.sleep(3)
    json.dump(rec,open(fp,'w'))
    print(f"{i+1:3d}/100 {r['id']:<28} {rec['status']}", flush=True)
    time.sleep(0.4)
print("done")
