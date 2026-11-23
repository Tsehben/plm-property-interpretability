#!/usr/bin/env python3
"""Additive biological annotation of the frozen development cohort; no test IO.

Commands: prepare, disorder-pilot, disorder, build. Remote taxonomy is separate.
"""
import argparse, hashlib, importlib.metadata, json, os, platform, time
from pathlib import Path
os.environ.setdefault('MPLCONFIGDIR','/tmp/uesolds-matplotlib')
import numpy as np
import pandas as pd
from biological_features import (AA,KD,A3V,GROUPS,HST,WINDOWS,canonical,biochemistry,
                                aggregation,charge_profile,regions,clip_regions,regional_summary,membrane_summary)
ROOT=Path(__file__).resolve().parent
BASE=ROOT.parent
OUT=ROOT/'biological_annotation'
IDENTITY=['raw_uid','legacy_uid','original_id','split','label','sha256','center','id_class','length','cluster']

def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda:f.read(8*1024*1024),b''): h.update(block)
    return h.hexdigest()

def read(path): return pd.read_csv(path,sep='\t',keep_default_na=False,low_memory=False)
def save(d,path): d.to_csv(path,sep='\t',index=False,float_format='%.12g',na_rep='NA')
def write_json(path,obj): path.write_text(json.dumps(obj,indent=2,sort_keys=True,allow_nan=False)+'\n')

def load():
    d=read(ROOT/'pfam_annotation/annotated_cohort.tsv').sort_values('raw_uid').reset_index(drop=True)
    m=read(ROOT/'results/cohort_manifest.tsv');m=m[m.clean&m.cached].sort_values('raw_uid').reset_index(drop=True)
    pd.testing.assert_frame_equal(d[IDENTITY],m[IDENTITY])
    assert len(d)==13880 and d.raw_uid.is_unique and d.split.value_counts().to_dict()=={'train':9880,'validation':4000}
    p=read(BASE/'esm2_probe/uesolds_pairs.tsv'); pairs=dict(p.iloc[:,:2].itertuples(index=False,name=None))
    sequences={r.raw_uid:pairs[r.legacy_uid] for r in d.itertuples()}
    for r in d.itertuples():
        s=sequences[r.raw_uid]
        assert len(s)==r.length and hashlib.sha256(s.encode()).hexdigest()==r.sha256
    return d,sequences

def prepare(out):
    out.mkdir(parents=True,exist_ok=True)
    d,s=load()
    paths=[ROOT/'results/cohort_manifest.tsv',ROOT/'results/records.tsv',
           ROOT/'results/diversity_cluster.tsv',BASE/'esm2_probe/uesolds_pairs.tsv',
           *[ROOT/'pfam_annotation'/f for f in ['annotated_cohort.tsv','domains.tsv','domains_flagged.tsv','memberships.tsv','search_run.json','Pfam-A.clans.tsv']]]
    inputs={str(p.relative_to(BASE)):sha(p) for p in paths}
    if (out/'input_hashes.json').exists(): assert json.loads((out/'input_hashes.json').read_text())==inputs
    write_json(out/'input_hashes.json',inputs)
    save(d[IDENTITY],out/'identity_manifest.tsv')
    # Representative operational pilot only. Never changes primary membership.
    chosen=set(d.sort_values('sha256').groupby(['split','center','label']).head(1).raw_uid)
    for q in np.linspace(0,1,11):
        chosen.add(d.iloc[(d.length-d.length.quantile(q)).abs().argmin()].raw_uid)
    for flag in ['terminal_his','terminal_cleavage','terminal_linker','esm_truncated']:
        chosen.update(d[d[flag]].sort_values('sha256').head(2).raw_uid)
    pilot=d[d.raw_uid.isin(chosen)].copy()
    save(pilot[IDENTITY+['terminal_his','terminal_cleavage','terminal_linker']],out/'pilot_manifest.tsv')
    write_json(out/'protocol.json',dict(version='1.0',cohort_size=13880,esm_residue_limit=1022,
        aa_order=AA,kd_scale=KD,aggregation_scale=A3V,aggregation_threshold=HST,aggregation_min_region=5,
        aggregation_windows={'le75':5,'le175':7,'le300':9,'gt300':11},aggregation_terminal_policy='paper; reversed-order example as sensitivity',
        groups=GROUPS,local_windows=WINDOWS,ph=[7.,7.4],pi='Biopython 1.86 Bjellqvist, bisection [0,14]',
        disorder_package='metapredict==3.1.1',disorder_network='V3',disorder_threshold=.5,
        disorder_segments='consecutive score >= 0.5, no gap closure, minimum 1',
        visible_regional='clip predictions made in full-sequence context; no predictor rerun on prefixes',
        positions='1-based inclusive; numpy offsets 0-based half-open',pilot_size=len(pilot)))
    print(f'Frozen {len(d)} proteins; pilot {len(pilot)}; {sum(map(len,s.values()))} residues',flush=True)

def check_inputs(out):
    for p,h in json.loads((out/'input_hashes.json').read_text()).items(): assert sha(BASE/p)==h,p

def disorder(out,pilot=False):
    import torch,metapredict
    assert importlib.metadata.version('metapredict')=='3.1.1'
    torch.set_num_threads(1);torch.manual_seed(42);torch.use_deterministic_algorithms(True)
    d,seqs=load();check_inputs(out)
    chosen=read(out/'pilot_manifest.tsv').raw_uid.tolist() if pilot else d.raw_uid.tolist()
    if not pilot: assert json.loads((out/'disorder_pilot.json').read_text())['passed']
    folder=out/'disorder_shards';folder.mkdir(exist_ok=True)
    start=time.monotonic();logs=[];maxdiff=0.
    for first in range(0,len(chosen),128):
        ids=chosen[first:first+128]; dest=folder/f'{"pilot" if pilot else "cohort"}_{first:05d}.npz'
        if dest.exists() and not pilot:
            z=np.load(dest,allow_pickle=False);assert z['raw_uid'].tolist()==ids
            assert z['sequence_sha256'].tolist()==[hashlib.sha256(seqs[u].encode()).hexdigest() for u in ids]
            continue
        arrays=[];offsets=[0];statuses=[]
        for uid in ids:
            s=seqs[uid];t=time.monotonic()
            if not canonical(s): v=np.full(len(s),np.nan);status='unsupported_input'
            else:
                # Explicit CPU, single-sequence inference eliminates batch-padding differences.
                try:
                    v=np.asarray(metapredict.predict(s,version='V3',use_device='cpu',normalized=True,round_values=False),dtype=np.float32)
                    assert v.shape==(len(s),) and np.isfinite(v).all() and ((v>=0)&(v<=1)).all()
                    if pilot:
                        repeat=np.asarray(metapredict.predict(s,version='V3',use_device='cpu',normalized=True,round_values=False))
                        maxdiff=max(maxdiff,float(np.max(np.abs(v-repeat))))
                    status='ok' if (v>=.5).any() else 'zero_features'
                except Exception as e:
                    v=np.full(len(s),np.nan,dtype=np.float32);status='failed:'+type(e).__name__
            arrays.append(v);offsets.append(offsets[-1]+len(s));statuses.append(status)
            logs.append(dict(raw_uid=uid,length=len(s),status=status,seconds=time.monotonic()-t))
        np.savez_compressed(dest,raw_uid=np.array(ids),sequence_sha256=np.array([hashlib.sha256(seqs[u].encode()).hexdigest() for u in ids]),
                            offsets=np.array(offsets),scores=np.concatenate(arrays),status=np.array(statuses))
        print(f'disorder {"pilot" if pilot else "cohort"}: {min(first+128,len(chosen))}/{len(chosen)}; {time.monotonic()-start:.1f}s',flush=True)
    package=Path(metapredict.__file__).parent
    weights={str(p.relative_to(package)):sha(p) for p in package.rglob('*') if p.suffix in ('.pt','.pth','.ckpt','.pickle','.pkl')}
    provenance=dict(package_version=importlib.metadata.version('metapredict'),reported_internal_version=metapredict.__version__,network='V3',
                    torch=torch.__version__,device='cpu',threads=1,weights=weights,code_sha256=sha(Path(__file__)),
                    elapsed_seconds=time.monotonic()-start,records=len(chosen),max_repeat_difference=maxdiff,
                    passed=all(not x['status'].startswith(('failed','unsupported')) for x in logs) and maxdiff<1e-6)
    write_json(out/('disorder_pilot.json' if pilot else 'disorder_run.json'),provenance)
    save(pd.DataFrame(logs),out/('disorder_pilot_timings.tsv' if pilot else 'disorder_run_timings.tsv'))
    check_inputs(out)

def load_disorder(out):
    result={}
    for p in sorted((out/'disorder_shards').glob('cohort_*.npz')):
        with np.load(p,allow_pickle=False) as archive:z={k:archive[k] for k in archive.files}
        for i,u in enumerate(z['raw_uid']):
            assert u not in result
            result[str(u)]=(z['scores'][z['offsets'][i]:z['offsets'][i+1]],str(z['status'][i]),str(z['sequence_sha256'][i]))
    return result

def build(out):
    d,seqs=load();check_inputs(out);disorder_data=load_disorder(out)
    domains=read(ROOT/'pfam_annotation/domains_flagged.tsv')
    domain_map={u:g.sort_values(['ali_from','ali_to','pfam_versioned']) for u,g in domains.groupby('raw_uid')}
    memberships=read(ROOT/'pfam_annotation/memberships.tsv')
    cluster_labels=d.groupby('cluster').label.nunique().to_dict()
    cluster_families=memberships.groupby('cluster').pfam.nunique().to_dict()
    tax_path=out/'taxonomy.tsv';tax=read(tax_path).set_index('raw_uid').to_dict('index') if tax_path.exists() else {}
    tm={}
    for completed in sorted((out/'transmembrane_shards').glob('*/completed.json')):
        meta=json.loads(completed.read_text())
        for name,h in meta['outputs'].items():assert sha(completed.parent/name)==h
        for x in meta['fingerprint']['sequences']:
            assert hashlib.sha256(seqs[x['raw_uid']].encode()).hexdigest()==x['sha256']
        for x in json.loads((completed.parent/'predictions.json').read_text()):
            if 'id' in x:
                assert x['id'] not in tm and len(x['topology_string'])==len(seqs[x['id']]);tm[x['id']]=x
    rows=[];statuses=[];intervals=[];agg_scores=[];offsets=[0]
    for idx,r in enumerate(d.itertuples()):
        uid=r.raw_uid;s=seqs[uid];n=len(s);visible=min(n,1022)
        row=r._asdict();row.pop('Index',None)
        row['pfam_release']=str(r.pfam_release)  # version string, never %.12g -> "37"
        # Old pilot organism flag is metadata, not uniform taxonomy success.
        row['legacy_pilot_organism_assigned']=row.pop('organism_assigned')
        for scope,seq in [('full',s),('esm',s[:1022])]:
            for key,value in biochemistry(seq).items():
                if scope=='esm' and key in ['pi_estimate','net_charge_ph7p0','net_charge_ph7p4','charge_density_ph7p0','charge_density_ph7p4']:
                    key='hypothetical_fragment_'+key
                row[scope+'_'+key]=value
        for ph in (7.,7.4):
            row['esm_original_context_charge_ph'+str(ph).replace('.','p')]=float(charge_profile(s,ph)[:visible].sum())
        def status(concept,value,tool,evidence,reason=''):
            row[concept+('_annotation_status' if concept=='pfam' else '_status')]=value
            statuses.append(dict(raw_uid=uid,concept=concept,status=value,tool=tool,evidence_type=evidence,reason=reason))
        status('biochemistry','ok' if canonical(s) else 'unsupported_input','biological_features-v1/Biopython1.86','direct_sequence_derived')
        av,ar,ast=aggregation(s);agg_scores.append(av);offsets.append(offsets[-1]+n)
        alternate,alternate_regions,_=aggregation(s,terminal_policy='published_example')
        row['aggregation_terminal_sensitivity_regions_changed']=ar!=alternate_regions
        row['aggregation_alternate_n_regions']=len(alternate_regions) if np.isfinite(alternate).all() else np.nan
        row['aggregation_alternate_fraction']=regional_summary(alternate,alternate_regions)['fraction']
        row['aggregation_length_extrapolation']=n>2000
        status('aggregation',ast,'AGGRESCAN-2007-independent-reimplementation-v1','biophysical_sequence_propensity')
        dv,dst,dsha=disorder_data.get(uid,(np.full(n,np.nan),'not_run',r.sha256))
        assert len(dv)==n and dsha==r.sha256
        dr=regions(dv>=.5) if np.isfinite(dv).all() else []
        status('disorder',dst,'metapredict3.1.1/networkV3','pretrained_model_prediction')
        for concept,v,rr in [('aggregation',av,ar),('disorder',dv,dr)]:
            for scope,limit in [('full',n),('esm',visible)]:
                for key,val in regional_summary(v[:limit],clip_regions(rr,limit)).items():row[f'{scope}_{concept}_{key}']=val
            for a,b in rr:
                intervals.append(dict(raw_uid=uid,concept=concept,start=a,end=b,mean_score=float(v[a-1:b].mean()),
                                      esm_start=a if a<=visible else np.nan,esm_end=min(b,visible) if a<=visible else np.nan,
                                      crosses_esm_boundary=a<=visible<b))
        dg=domain_map.get(uid)
        row['pfam_architecture_status']='no_hit' if dg is None else 'ambiguous_overlap' if dg.clan_competing_overlap.any() else 'assigned'
        row['full_domain_architecture']='' if dg is None else ';'.join(dg.pfam_versioned)
        row['full_domain_architecture_coordinates']='' if dg is None else ';'.join(f'{h.pfam_versioned}@{h.ali_from}-{h.ali_to}' for h in dg.itertuples())
        vg=None if dg is None else dg[dg.ali_from<=visible]
        row['esm_domain_architecture']='' if vg is None else ';'.join(vg.pfam_versioned)
        row['full_domain_count']=0 if dg is None else len(dg)
        row['esm_domain_count']=0 if vg is None else len(vg)
        row['esm_partial_domains']=0 if vg is None else int((vg.ali_to>visible).sum())
        row['domain_competing_overlap']=False if dg is None else bool(dg.clan_competing_overlap.any())
        domain_mask=np.zeros(n,dtype=bool)
        if dg is not None:
            for h in dg.itertuples():domain_mask[h.ali_from-1:h.ali_to]=True
        row['full_domain_covered_fraction']=domain_mask.mean()
        row['esm_domain_covered_fraction']=domain_mask[:visible].mean()
        status('pfam','no_hit' if dg is None else 'ok','Pfam37.0/PyHMMER0.11.2','external_biological_annotation')
        row['cluster_label_count']=cluster_labels[r.cluster];row['cluster_family_count']=cluster_families.get(r.cluster,0)
        row['cluster_spans_labels']=row['cluster_label_count']>1
        row['cluster_spans_families']=row['cluster_family_count']>1
        status('homology','ok','existing_MMseqs2_30pct_80pct','external_biological_annotation')
        tr=tax.get(uid,dict(taxonomy_status='not_run'))
        row.update({k:v for k,v in tr.items() if k not in IDENTITY})
        lineage=tr.get('taxonomy_lineage','')
        row['signalp_organism_setting']='unknown'
        if tr['taxonomy_status']=='assigned_exact_species_consensus':
            if ':Eukaryota' in lineage:row['signalp_organism_setting']='eukarya'
            elif ':Bacteria' in lineage or ':Archaea' in lineage:row['signalp_organism_setting']='other'
        status('taxonomy',tr['taxonomy_status'],'UniParc_exact_full_sequence','external_biological_annotation')
        # Explicit unavailability is not a negative biological prediction.
        for concept,tool,reason in [
            ('secondary_structure','NetSurfP3.0','registered_package_and_checkpoint_unavailable'),
            ('accessibility','NetSurfP3.0','registered_package_and_checkpoint_unavailable'),
            ('signal_peptide','SignalP6.0','registered_package_unavailable;kingdom_must_be_verified')]:
            status(concept,'unavailable',tool,'pretrained_model_prediction',reason)
        for key in ['helix_fraction','sheet_fraction','coil_fraction','rsa_mean','exposed_fraction','tm_segment_count','tm_fraction','membrane_associated_fraction','signal_peptide_presence','signal_peptide_confidence','signal_peptide_cleavage_site']:
            row['full_'+key]=row['esm_'+key]=np.nan
        row['signal_peptide_type']='';row['tm_classification']=''
        prediction=tm.get(uid)
        if prediction is None:
            status('transmembrane','not_run_compute','DeepTMHMM2/526f6c60','pretrained_model_prediction','Uniform resumable CPU job; record not yet in a completed shard. See transmembrane_progress.json.')
        else:
            segments=prediction['segments'];trr=[(a,b) for name,a,b in segments if name in ('TMhelix','Beta sheet')]
            row['tm_classification']=prediction['type']
            for scope,limit in [('full',n),('esm',visible)]:
                for k,v in membrane_summary(segments,limit).items():row[scope+'_'+k]=v
            for name,a,b in segments:
                if name not in ('TMhelix','Beta sheet','reentrant','interfacial','signal','transit_peptide'):continue
                intervals.append(dict(raw_uid=uid,concept='DeepTMHMM2_'+name,start=a,end=b,mean_score=np.nan,
                    esm_start=a if a<=visible else np.nan,esm_end=min(b,visible) if a<=visible else np.nan,crosses_esm_boundary=a<=visible<b))
            status('transmembrane','ok' if trr else 'zero_features','DeepTMHMM2/526f6c60','pretrained_model_prediction')
        rows.append(row)
        if (idx+1)%2000==0:print(f'biochemistry/merge {idx+1}/{len(d)}',flush=True)
    result=pd.DataFrame(rows);pd.testing.assert_frame_equal(result[IDENTITY],d[IDENTITY])
    save(result,out/'protein_annotations.tsv')
    save(pd.DataFrame(statuses),out/'annotation_status.tsv')
    reg=pd.DataFrame(intervals);assert (reg.start>=1).all()
    assert (reg.end<=reg.raw_uid.map(d.set_index('raw_uid').length)).all()
    save(reg,out/'regions.tsv')
    # Reuse domain table without precision loss in original scores.
    raw=read(ROOT/'pfam_annotation/domains.tsv')
    assert raw.raw_uid.tolist()==domains.raw_uid.tolist()
    raw['clan']=domains.clan;raw['clan_competing_overlap']=domains.clan_competing_overlap
    raw['esm_visible']=raw.ali_from<=1022;raw['esm_partially_visible']=(raw.ali_from<=1022)&(raw.ali_to>1022)
    raw.to_csv(out/'pfam_domains.tsv',sep='\t',index=False,float_format='%.17g')
    np.savez_compressed(out/'aggregation_profiles.npz',raw_uid=d.raw_uid.to_numpy(dtype=str),sequence_sha256=d.sha256.to_numpy(dtype=str),offsets=np.array(offsets),scores=np.concatenate(agg_scores))
    # Sequence code is a compact, lossless positional record for KD/charge profiles.
    np.savez_compressed(out/'sequence_profiles.npz',raw_uid=d.raw_uid.to_numpy(dtype=str),sequence_sha256=d.sha256.to_numpy(dtype=str),offsets=np.array(offsets),aa_ascii=np.frombuffer(''.join(seqs[u] for u in d.raw_uid).encode(),dtype=np.uint8))
    check_inputs(out)
    write_json(out/'build_run.json',dict(records=len(result),splits=result.split.value_counts().to_dict(),residues=int(d.length.sum()),
        nonstandard_records=int((result.full_nonstandard_count>0).sum()),esm_truncated_records=int((d.length>1022).sum()),
        input_identity_checks_passed=True,coordinates_checked=True,script_sha256=sha(Path(__file__)),features_sha256=sha(ROOT/'biological_features.py'),
        python=platform.python_version(),packages={k:importlib.metadata.version(k) for k in ['numpy','pandas','scipy','biopython','metapredict','torch']},
        transmembrane_cohort_predictions=len(tm),missing_tools=['NetSurfP3.0','SignalP6.0']+(['DeepTMHMM2_uniform'] if len(tm)<len(d) else []),all_requested_annotations_complete=False))
    print('Annotation table written; unavailable predictors remain explicitly missing.',flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('command',choices=['prepare','disorder-pilot','disorder','build']);p.add_argument('--out',type=Path,default=OUT);a=p.parse_args()
    if a.command=='prepare':prepare(a.out)
    elif a.command.startswith('disorder'):disorder(a.out,pilot=a.command=='disorder-pilot')
    else:build(a.out)
