#!/usr/bin/env python3
"""Check full-run accounting and one-thread reproducibility on first 50 profiles."""
import argparse,itertools,json
from pathlib import Path
import numpy as np
import pandas as pd
import pyhmmer
from annotate_pfam import ROOT,BASE,digest,text

def run(ref,out):
    a=pd.read_csv(out/'searched_manifest.tsv',sep='\t',keep_default_na=False)
    pairs=pd.read_csv(BASE/'esm2_probe/uesolds_pairs.tsv',sep='\t',keep_default_na=False)
    seqmap=dict(pairs.iloc[:,:2].itertuples(index=False,name=None))
    alphabet=pyhmmer.easel.Alphabet.amino()
    block=pyhmmer.easel.DigitalSequenceBlock(alphabet,[pyhmmer.easel.TextSequence(name=r.raw_uid.encode(),sequence=seqmap[r.legacy_uid]).digitize(alphabet) for r in a.itertuples()])
    # Read the full HMM iterator to avoid closing gzip midstream; retain first 50 only.
    profiles=[]
    with pyhmmer.plan7.HMMFile(ref/'Pfam-A.hmm.gz') as f:
        for i,h in enumerate(f):
            if i<50:profiles.append(h)
    rows=[]
    for hits in pyhmmer.hmmer.hmmsearch(profiles,block,cpus=1,parallel='queries',bit_cutoffs='gathering',seed=42,Z=13880,F1=.02,F2=.001,F3=.00001,bias_filter=True,null2=True):
        for hit in hits:
            if not hit.included:continue
            for dom in hit.domains:
                if dom.included:rows.append(dict(raw_uid=text(hit.name),pfam_versioned=text(hits.query.accession),ali_from=dom.alignment.target_from,ali_to=dom.alignment.target_to,domain_score=dom.score))
    fields=['raw_uid','pfam_versioned','ali_from','ali_to'];cols=fields+['domain_score']
    observed=pd.DataFrame(rows,columns=cols).sort_values(fields).reset_index(drop=True)
    expected=pd.read_csv(out/'domains.tsv',sep='\t')
    expected=expected[expected.pfam_versioned.isin([text(h.accession) for h in profiles])][cols].sort_values(fields).reset_index(drop=True)
    assert observed[fields].equals(expected[fields])
    np.testing.assert_allclose(observed.domain_score,expected.domain_score,rtol=0,atol=1e-10)
    d=pd.read_csv(out/'annotated_cohort.tsv',sep='\t',keep_default_na=False)
    m=pd.read_csv(out/'memberships.tsv',sep='\t',keep_default_na=False)
    assert len(d)==13880 and d.raw_uid.is_unique and d.annotation_attempted.all()
    assert set(d[d.pfam_status=='assigned'].raw_uid)==set(m.raw_uid)
    assert not m.duplicated(['raw_uid','pfam']).any()
    s=pd.read_csv(out/'feasibility_summary.tsv',sep='\t');s=s[s.policy=='GA_all']
    def counts(mm):return pd.crosstab(mm.pfam,mm.label).reindex(columns=[0,1],fill_value=0)
    pooled=counts(m);tr=counts(m[m.split=='train']);va=counts(m[m.split=='validation'])
    for r in s.itertuples():
        train=set(tr.index[(tr[0]>=r.k)&(tr[1]>=r.k)]);valid=set(va.index[(va[0]>=r.k)&(va[1]>=r.k)])
        both=set(va.index[(va[0]>0)&(va[1]>0)])
        assert len(train)==r.training_eligible
        assert len(train&valid)==r.both_splits_threshold
        assert len(train&both)==r.training_eligible_val_both
        assert int(((pooled[0]>=r.k)&(pooled[1]>=r.k)).sum())==r.pooled_qualifying
    meta=json.loads((out/'search_run.json').read_text())
    assert all(digest(BASE/p)==h for p,h in meta['inputs'].items())
    result=dict(first_50_profiles_single_thread_domains=len(observed),single_thread_scores_match=True,all_13880_accounted=True,membership_uniqueness=True,independent_crosstab_counts_match=True,input_hashes_unchanged=True,validation_script_sha256=digest(Path(__file__)))
    (out/'validation_checks.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--reference',type=Path,required=True);ap.add_argument('--out',type=Path,default=ROOT/'pfam_annotation');a=ap.parse_args();run(a.reference,a.out)
