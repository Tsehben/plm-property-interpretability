import json,os,sys,re,collections
OUT=sys.argv[1]
S=json.load(open(os.path.join(OUT,'sample_100.json')))
P=json.load(open(os.path.join(OUT,'uniparc_parsed.json')))
U=json.load(open(os.path.join(OUT,'uniprot_entries.json')))
T=json.load(open(os.path.join(OUT,'taxnames.json')))
R=json.load(open(os.path.join(OUT,'rcsb.json')))
TRIM={}
d=os.path.join(OUT,'uniparc_trim_cache')
for f in os.listdir(d) if os.path.isdir(d) else []:
    r=json.load(open(os.path.join(d,f)))
    if r.get('hit'): TRIM[r['id']]=r
JUNK={'32630','32644','6666666','2','81077'}  # synthetic construct, unidentified, invalid, generic
def sp(t):
    i=T.get(t); return (i['species_id'],i['species_name']) if i else (t,'taxid:'+t)
# prefix organism hints (naming convention, UNVERIFIED)
HINT={'Lmaj':'Leishmania major','Tcru':'Trypanosoma cruzi','Pfal':'Plasmodium falciparum',
 'Tbru':'Trypanosoma brucei','Pviv':'Plasmodium vivax','Pfu':'Pyrococcus furiosus',
 'BobuA':'Borrelia/Borreliella sp.','BaheA':'Bartonella henselae','BuxeA':'Burkholderia xenovorans',
 'CrneC':'Cryptococcus neoformans','LebrA':'Leishmania braziliensis','LyraC':'Leishmania?'}
def hint(pid):
    m=re.match(r'^SGPP-([A-Za-z]{4})',pid)
    if m and m.group(1) in HINT: return HINT[m.group(1)]
    m=re.match(r'^SSGCID-([A-Za-z]{5})\.',pid)
    if m and m.group(1) in HINT: return HINT[m.group(1)]
    m=re.match(r'^SECSG-([A-Za-z]{3})-',pid)
    if m and m.group(1) in HINT: return HINT[m.group(1)]
    return ''
rows=[];stats=collections.Counter()
for r in S:
    pid=r['id']; e=P[pid]
    used_trim=pid in TRIM
    if used_trim:
        tr=TRIM[pid]['hit']; res=tr['result']
        e=dict(e); e['upi']=res.get('accession'); e['up_acc']=[];e['up_active']=[];e['up_tax']=set();e['all_tax']=set();e['pdb']=[]
        for x in res.get('dbReference',[]):
            pr={p['type']:p['value'] for p in x.get('property',[])}
            if pr.get('NCBI_taxonomy_id'): e['all_tax'].add(pr['NCBI_taxonomy_id'])
            if x['type'].startswith('UniProtKB'):
                e['up_acc'].append(x['id'])
                if x.get('active')=='Y':
                    e['up_active'].append(x['id'])
                    if pr.get('NCBI_taxonomy_id'): e['up_tax'].add(pr['NCBI_taxonomy_id'])
            if x['type']=='PDB': e['pdb'].append(x['id'])
        e['up_acc']=sorted(set(e['up_acc']));e['up_active']=sorted(set(e['up_active']))
        e['up_tax']=sorted(e['up_tax']);e['all_tax']=sorted(e['all_tax']);e['pdb']=sorted(set(e['pdb']))
    # --- ORGANISM ---
    ev=[];org='';status='';accs=e['up_active'] or e['up_acc']
    if not e['upi']:
        status='no_exact_match'; org=''
        ev.append('none')
    else:
        upsp={sp(t) for t in e['up_tax'] if t not in JUNK}
        if len(upsp)==1:
            org=list(upsp)[0][1]; status='confident'; ev.append('UniParc_exact->UniProtKB_active')
        elif len(upsp)>1:
            org='; '.join(sorted(n for _,n in upsp)); status='ambiguous_multi_organism'
            ev.append('UniParc_exact->UniProtKB_active(multi)')
        else:
            # no UniProtKB: fall back to RCSB source organism of matching entries
            rs=set()
            for p in e['pdb']:
                en=R.get(p.lower())
                if not en: continue
                pes=en['polymer_entities'] or []
                for pe in pes:
                    uids=(pe.get('rcsb_polymer_entity_container_identifiers') or {}).get('uniprot_ids') or []
                    if len(pes)==1 or (set(uids)&set(accs)):
                        for so in (pe.get('rcsb_entity_source_organism') or []):
                            tid=str(so.get('ncbi_taxonomy_id') or '')
                            if not tid or tid in JUNK: continue
                            if tid in T: rs.add(sp(tid))
                            else: rs.add((tid, so.get('ncbi_scientific_name') or ('taxid:'+tid)))
            allsp={sp(t) for t in e['all_tax'] if t not in JUNK}
            if len(rs)==1: org=list(rs)[0][1]; status='confident'; ev.append('UniParc_exact->PDB_source_organism')
            elif len(rs)>1: org='; '.join(sorted(n for _,n in rs)); status='ambiguous_multi_organism'; ev.append('PDB_source_organism(multi)')
            elif len(allsp)==1: org=list(allsp)[0][1]; status='confident'; ev.append('UniParc_exact->EMBL/other_xref')
            elif len(allsp)>1:
                org='; '.join(sorted(n for _,n in allsp)[:6])+(f' (+{len(allsp)-6} more)' if len(allsp)>6 else '')
                status='ambiguous_multi_organism'; ev.append(f'UniParc_exact->xrefs({len(allsp)}_species)')
            else: status='match_no_taxonomy'; ev.append('UniParc_exact_no_tax')
    # --- EXPRESSION HOST (kept separate) ---
    hosts=set()
    for p in e['pdb']:
        en=R.get(p.lower())
        if not en: continue
        pes=en['polymer_entities'] or []
        for pe in pes:
            uids=(pe.get('rcsb_polymer_entity_container_identifiers') or {}).get('uniprot_ids') or []
            if len(pes)==1 or (set(uids)&set(accs)):
                for ho in (pe.get('rcsb_entity_host_organism') or []):
                    if ho.get('ncbi_scientific_name'): hosts.add(ho['ncbi_scientific_name'])
    vec=''
    m=re.search(r'_(pET\d+)$',pid)
    if m: vec=m.group(1)
    # --- PFAM ---
    pf={}
    for a in accs:
        for pid2,nm in (U.get(a,{}).get('pfam') or []): pf[pid2]=nm
    pf_src='UniProt_xref' if pf else ''
    if not pf:
        for p in e['pdb']:
            en=R.get(p.lower())
            if not en: continue
            pes=en['polymer_entities'] or []
            for pe in pes:
                uids=(pe.get('rcsb_polymer_entity_container_identifiers') or {}).get('uniprot_ids') or []
                if len(pes)==1 or (set(uids)&set(accs)):
                    for a2 in (pe.get('rcsb_polymer_entity_annotation') or []):
                        if a2['type']=='Pfam': pf[a2['annotation_id']]=a2.get('name')
        if pf: pf_src='RCSB_PDB_entity_annotation'
    if pf: ev.append('Pfam:'+pf_src)
    stats[status]+=1
    rows.append({'original_id':pid,'id_class':r['cls'],'raw_label':r['label'],'seq_len':r['len'],
      'uniparc_upi':e['upi'] or '','matched_accessions':';'.join(accs) or '',
      'n_uniprot_xrefs':len(e['up_acc']),'pdb_xrefs':';'.join(e['pdb'][:6])+('...' if len(e['pdb'])>6 else ''),
      'source_organism':org,'organism_status':status,
      'pfam_accessions':';'.join(sorted(pf)) or '','pfam_names':';'.join(pf[k] or '' for k in sorted(pf)) or '',
      'expression_host':';'.join(sorted(hosts)) or '','expression_vector_from_id':vec,
      'organism_hint_from_id_unverified':hint(pid),
      'match_method':'tag_trimmed_exact:'+TRIM[pid]['hit']['why'] if used_trim else ('exact_full_sequence' if e['upi'] else 'none'),
      'evidence_source':'|'.join(ev)})
cols=list(rows[0].keys())
tsv=os.path.join(OUT,'uesolds_pilot_annotations.tsv')
with open(tsv,'w') as f:
    f.write('\t'.join(cols)+'\n')
    for r in rows: f.write('\t'.join(str(r[c]).replace('\t',' ') for c in cols)+'\n')
print("wrote",tsv)
print("\norganism_status:",dict(stats))
print("with >=1 Pfam:",sum(1 for r in rows if r['pfam_accessions']))
json.dump(rows,open(os.path.join(OUT,'rows.json'),'w'),indent=1)
