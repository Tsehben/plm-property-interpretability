#!/usr/bin/env python3
"""Fetch Pfam entry-type metadata from InterPro. ACCESSIONS AND METADATA ONLY.

No protein sequence is transmitted. The request is a paginated listing of the Pfam
member-database entry catalogue; the query string contains no sequence, no protein
accession, and no identifying information. Results are cached so the fetch runs once.

Usage: fetch_pfam_entry_types.py <out_dir>
"""
import json, os, sys, time, urllib.request, urllib.error

OUT = sys.argv[1]
CACHE = os.path.join(OUT, 'cache', 'pfam_entry_types.json')
os.makedirs(os.path.dirname(CACHE), exist_ok=True)
UA = 'plm-property-interpretability/1.0 (academic research; entry metadata only)'

db = json.load(open(CACHE)) if os.path.exists(CACHE) else {}
url = 'https://www.ebi.ac.uk/interpro/api/entry/pfam/?page_size=200'
page, t0 = 0, time.time()
while url:
    for attempt in range(5):
        try:
            req = urllib.request.Request(url, headers={'User-Agent': UA, 'Accept': 'application/json'})
            with urllib.request.urlopen(req, timeout=60) as r:
                d = json.loads(r.read().decode())
            break
        except Exception as e:
            if attempt == 4: raise
            time.sleep(2 ** attempt)
    for e in d['results']:
        m = e['metadata']
        db[m['accession']] = {'name': m['name'], 'type': m['type'],
                              'integrated_interpro': m.get('integrated')}
    page += 1
    url = d.get('next')
    if page % 25 == 0 or not url:
        json.dump(db, open(CACHE, 'w'))
        print(f'  page {page}: {len(db):,} entries cached  [{time.time()-t0:,.0f}s]', flush=True)
    time.sleep(0.15)
json.dump(db, open(CACHE, 'w'))
import collections
print(f'done: {len(db):,} Pfam entries')
print('type distribution:', dict(collections.Counter(v['type'] for v in db.values())))
