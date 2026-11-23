#!/usr/bin/env python3
"""Additional uncertainty for the fixed cross-centre transfer matrix; no new fits."""
from source_validity_audit import OUT,read,save,binary_summary,bootstrap_indices,multi
import numpy as np
import pandas as pd

def main():
 p=read(OUT/'solubility_predictions.tsv');v=read(OUT/'evaluation_manifest.tsv');rows=[]
 for (protocol,model),q in p[p.protocol.str.startswith('only_')].groupby(['protocol','model']):
  pp=q.set_index('raw_uid').loc[v.raw_uid].probability.to_numpy()
  for center,g in v.groupby('center'):
   ix=g.index.to_numpy();rows.append(dict(protocol=protocol,model=model,evaluation=center,**binary_summary(g,pp[ix])))
 save(pd.DataFrame(rows),'cross_centre_uncertainty.tsv')
 # Use one shared validation-group bootstrap for each candidate model's accuracy.
 # Source metrics already include all-class uncertainty; conditional family point
 # summaries are not independent tests and are not used to select cohorts.
 print('Cross-centre uncertainty complete')
if __name__=='__main__':main()


def conditional_source_uncertainty():
 from source_validity_audit import ROOT
 support=read(OUT/'source_within_family_support.tsv');v=read(OUT/'evaluation_manifest.tsv')
 members=read(ROOT/'pfam_annotation/memberships.tsv')
 probs={name:read(OUT/f'source_predictions_{name}.tsv').set_index('raw_uid').loc[v.raw_uid] for name in ['training_prior','composition','Pfam','ESM_layer8','composition_plus_ESM']}
 rows=[]
 for fam,g in support.groupby('pfam'):
  centres=sorted(g.center);uids=set(members[members.pfam==fam].raw_uid)
  ix=np.flatnonzero(v.raw_uid.isin(uids)&v.center.isin(centres));vv=v.iloc[ix];y=np.array([centres.index(c) for c in vv.center]);draws=bootstrap_indices(vv.cluster)
  prior=g.set_index('center').loc[centres].train_n.to_numpy(float);prior/=prior.sum()
  for name in list(probs)+['family_training_prior']:
   pp=np.tile(prior,(len(ix),1)) if name=='family_training_prior' else probs[name].iloc[ix][centres].to_numpy(float)
   pp=pp/pp.sum(1,keepdims=True)
   point=multi(y,pp,len(centres));boot=[multi(y[j],pp[j],len(centres)) for j in draws if len(set(y[j]))==len(centres)]
   for metric in ['accuracy','balanced_accuracy','macro_ovr_auroc','logloss']:
    lo,hi=np.quantile([b[metric] for b in boot],[.025,.975]) if len(boot)>=20 else [np.nan,np.nan]
    rows.append(dict(pfam=fam,model=name,metric=metric,value=point[metric],low=lo,high=hi,valid_bootstraps=len(boot),n=len(ix)))
 save(pd.DataFrame(rows),'source_within_family_uncertainty.tsv')
 print('Conditional source uncertainty complete')

if __name__=='__main__':conditional_source_uncertainty()
