#!/usr/bin/env python3
"""Render measured tables and a per-centre plot after pipeline.py finish."""
from pathlib import Path
import os
os.environ.setdefault('MPLCONFIGDIR','/tmp/uesolds-matplotlib')
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
P=Path(__file__).resolve().parent
D=P/'results'
def table(d):
    def val(v):
        if pd.isna(v): return 'undefined'
        if isinstance(v,float): return f'{v:.4f}'
        return str(v)
    return '| '+' | '.join(d.columns)+' |\n| '+' | '.join(['---']*len(d.columns))+' |\n'+'\n'.join('| '+' | '.join(map(val,r))+' |' for r in d.itertuples(index=False,name=None))+'\n'
b=pd.read_csv(D/'baseline_results.tsv',sep='\t')
c=pd.read_csv(D/'cohort_statistics.tsv',sep='\t')
import json
s=json.loads((D/'summary.json').read_text())
text='\n## Measured before/after comparison\n\n'
text+=f"The audit found {s['exact_duplicate_records']} exact-duplicate records and {s['conflicting_records']} records with identical-sequence label conflicts across train/development. There are {s['repeated_identifier_records']} records with reused identifiers. MMseqs2 found {s['cross_split_pairs']} cross-split pairs involving {s['train_with_hits']} training and {s['validation_with_hits']} validation records. All development records occupy {s['clusters']} detected clusters. The cache mapping detects {s['cache_index_displacements']} displaced legacy indices.\n\n"
text+='Full-pool cohort composition (cluster counts refer to the same all-development clustering):\n\n'+table(c[c.scope=='full'][['cohort','split','n','insoluble','soluble','unique_sequences','clusters','centers']])
text+='\nMatched cached-embedding comparison sizes:\n\n'+table(c[c.scope=='cached'][['cohort','split','n','insoluble','soluble']])
text+='\nAll three models are refit per cohort on exactly those cached training rows, with training-only preprocessing. No ESM weights are loaded or updated. The original row below is the previously frozen 10,000-training comparison, not a new full-70,030-training ESM experiment.\n\n'
v=b[b.group.isin(['POOLED','MACRO_training_selected_centers'])].pivot(index=['cohort','group'],columns='model',values='auroc').reset_index();v.columns.name=None
text+=table(v)
text+='\nPer-centre source-controlled validation results (single-class AUROC is undefined):\n\n'
v=b[(b.cohort=='source_controlled')&~b.group.isin(['POOLED','MACRO_training_selected_centers'])]
a=v.pivot(index=['group','n0','n1'],columns='model',values='auroc').reset_index();a.columns.name=None
a[['n0','n1']]=a[['n0','n1']].astype(int)
text+=table(a)
text+='\nThe reported training-selected macro averages **six defined centres out of eight**. NESG and NYSGRC have no positive validation records; SSGCID has only one negative. Its large score contribution makes the six-centre macro unstable; use the historical five-centre sensitivity and per-centre counts to interpret it. The TSV reports every source for every cohort, along with AP and sample sizes. The macro averages only centres with both validation labels. Small strata have high uncertainty; no significance claim or cohort selection is based on these scores. Different pooled populations are not paired comparisons.\n\n'
# Historical five centres on each cohort's same evaluation rows.
five=['EFI','JCSG','MCSG','SECSG','SGPP']
f=b[b.group.isin(five)].groupby(['cohort','model']).auroc.mean().unstack().reset_index();f.columns.name=None
text+='Historical five-centre macro for continuity (same named centres, refit training per cohort):\n\n'+table(f)
text+='\nThe original ESM five-centre mean reproduces to about 0.0001 across the old and new software environments (0.5174 historically versus 0.5173 here); environment versions are recorded. This numerical difference is not a biological change.\n'
text+='\n![Per-centre source-controlled probes](results/within_center_probes.png)\n'
fig,ax=plt.subplots(figsize=(10,4)); plot=a.set_index('group')[['source','composition','ESM_layer8']];plot.plot.bar(ax=ax,color=['#777777','#d98268','#4e8ca8']);ax.axhline(.5,color='black',ls='--',lw=1);ax.set_ylim(0,1);ax.set_ylabel('Development AUROC');ax.set_xlabel('Centre (undefined single-class scores omitted)');fig.tight_layout();fig.savefig(D/'within_center_probes.png',dpi=160);plt.close(fig)
text+='''
## Recommendation and remaining work

Use the **cached source-controlled cohort after homology filtering** for a bounded exploratory SAE comparison, retaining the original and cleaned-all-source results as diagnostics. This choice follows provenance support and leakage controls, not whichever AUROC is largest. Use the exact `raw_uid`/`legacy_uid` mapping and unchanged sequences to reuse ESM representations. Report each centre, its class counts, and the macro denominator; preserve the historical five-centre score as a secondary continuity measure.

Do **not** label this cohort experimentally verified or use it alone to claim direct solubility representations. Source control reduces a measured shortcut and train/development homology filtering addresses detectable leakage, but neither repairs workflow-derived labels. Strict-label comparison is explicitly unavailable (zero verified rows, no model fit). Recover construct-matched DNASU statuses or original quantitative eSOL measurements before making stronger solubility claims. Lack of matches or scores must remain unknown. This is the main scientific blocker, not insufficient probe AUROC.

The six-week scope supports one frozen cohort, one layer and paired concept/property analyses with these limitations. Avoid a large annotation/relabeling project or a model/threshold sweep. Family and organism controls remain limited by pilot coverage, and repeated constructs and assay heterogeneity remain possible biases. Future uncertainty estimates should resample homology clusters within centres; the reported point estimates alone do not establish an ESM advantage. Earlier statements that homology would inflate every model equally are not justified: different feature sets may exploit relatives differently.

Validation checks passed for original input preservation, sequence/hash cache alignment, finite model inputs, optimizer convergence and cleaned split separation by exact sequence, shared identifier, detected cross-split hit and cluster. No independent test file was opened. Generated inputs are reproducible and ignored by Git; manifests, compact evidence tables, plots and run provenance are retained.
'''
report=P/'EDA_REPORT.md'
original=report.read_text().split('\n## Measured before/after comparison')[0]
report.write_text(original+text)
print(text[:1500])

# Freeze generated evidence and implementation provenance after rendering.
import hashlib
(D/'output_hashes.json').write_text(json.dumps({f.name:hashlib.sha256(f.read_bytes()).hexdigest() for f in sorted(D.iterdir()) if f.is_file() and f.name not in ['output_hashes.json','implementation_hashes.json']},indent=2)+'\n')
(D/'implementation_hashes.json').write_text(json.dumps({f.name:hashlib.sha256(f.read_bytes()).hexdigest() for f in [P/'pipeline.py',P/'summarize.py',P/'README.md',P/'EDA_REPORT.md']},indent=2)+'\n')
