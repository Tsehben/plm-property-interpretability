#!/usr/bin/env python3
"""Additional descriptive coverage and evidence flags; no model fits or cohort changes."""
from pathlib import Path
import hashlib,json
import numpy as np
import pandas as pd
from sklearn.preprocessing import StandardScaler
from sklearn.neighbors import NearestNeighbors
from threadpoolctl import threadpool_limits
from source_validity_audit import ROOT,OUT,read,save,AA

def main():
 d=read(ROOT/'pfam_annotation/annotated_cohort.tsv').sort_values('raw_uid').reset_index(drop=True)
 rec=read(ROOT/'results/records.tsv').set_index('raw_uid').loc[d.raw_uid].reset_index()
 assert np.array_equal(rec.sha256,d.sha256)
 x=np.c_[np.log(d.length),rec[['frac_'+a for a in AA]].to_numpy()]
 scale=StandardScaler().fit(x[d.split=='train']);x=scale.transform(x)
 rows=[]
 with threadpool_limits(limits=2):
  # Max cluster size + 1 guarantees at least one different detected group in k neighbours.
  k=int(d.cluster.value_counts().max())+1
  dist,ix=NearestNeighbors(n_neighbors=k,algorithm='brute').fit(x).kneighbors(x)
  same=d.cluster.to_numpy()[ix]==d.cluster.to_numpy()[:,None]
  base=np.where(same,np.inf,dist).min(1);assert np.isfinite(base).all()
  for source in sorted(d.center.unique()):
   selected=np.flatnonzero(d.center.eq(source));retained=np.flatnonzero(d.center.ne(source))
   ds,js=NearestNeighbors(n_neighbors=k,algorithm='brute').fit(x[retained]).kneighbors(x[selected])
   same=d.cluster.to_numpy()[retained[js]]==d.cluster.to_numpy()[selected,None]
   after=np.where(same,np.inf,ds).min(1);assert np.isfinite(after).all()
   increase=after-base[selected];assert increase.min()>-1e-6
   rows.append(dict(center=source,n=len(selected),baseline_nearest_other_group_median=np.median(base[selected]),
                    retained_nearest_other_group_median=np.median(after),increase_median=np.median(increase),increase_p95=np.quantile(increase,.95),
                    fraction_distance_increased=float((increase>1e-6).mean())))
 save(pd.DataFrame(rows),'source_biochemical_neighbour_loss.tsv')
 # Original source labels remain frozen; expose subprefixes hidden inside 'other'.
 other=d[d.center=='other'].copy();other['header_subprefix']=other.original_id.str.split('-').str[0]
 save(other.groupby(['header_subprefix','split','label']).size().reset_index(name='n'),'other_header_subprefixes.tsv')
 # Separate recovered records/statuses from eligibility for a reliable binary study.
 manifest=read(OUT/'policy_manifest.tsv')
 tt=read(OUT/'targettrack_sample_matches.tsv').set_index('raw_uid')
 trials=read(OUT/'targettrack_sample_trials.tsv');dn=read(OUT/'dnasu_sample_matches.tsv').set_index('raw_uid')
 manifest['record_evidence']='not individually audited'
 manifest['exact_archive_trial_sequence']=False;manifest['archive_explicit_soluble_status']=False
 manifest['archive_explicit_negative_measurement']=False;manifest['dnasu_explicit_status_matches_label']=False
 for i,r in manifest.iterrows():
  if r.raw_uid in tt.index:
   q=tt.loc[r.raw_uid];exact=bool(str(q.exact_trial_ids));hist=str(q.exact_trial_statuses).split(';')
   manifest.loc[i,'record_evidence']='TargetTrack ID and target sequence recovered; '+('exact trial sequence' if exact else 'trial sequence differs')
   manifest.loc[i,'exact_archive_trial_sequence']=exact
   manifest.loc[i,'archive_explicit_soluble_status']='soluble' in hist
  if r.raw_uid in dn.index:
   q=dn.loc[r.raw_uid];manifest.loc[i,'record_evidence']='DNASU explicit status recovered; assayed fusion versus released construct unresolved'
   manifest.loc[i,'dnasu_explicit_status_matches_label']=bool(q.published_label_agrees)
 manifest['property_eligibility_reason']='No validated comparable binary assay cohort established; recovered statuses are evidence tiers, not eligibility certification'
 save(manifest,'policy_manifest.tsv')
 (OUT/'supplement_run.json').write_text(json.dumps(dict(script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    features='full sequence log length and 20 AA fractions; train-only StandardScaler',distance='Euclidean; nearest protein outside query MMseqs group',
    note='Descriptive coverage sensitivity, not a population support bound; no labels used for geometry'),indent=2)+'\n')
 print('Feature-space loss and record-evidence flags complete')
if __name__=='__main__':main()
