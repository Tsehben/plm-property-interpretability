#!/usr/bin/env python3
"""Print figures for LOCALISATION_FINAL_CHECKIN.pdf. Saved results only."""
import csv, json, collections, os
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

OUT='report_build/figs'; os.makedirs(OUT, exist_ok=True)
S  = dict(surface='#ffffff', ink='#0b0b0b', ink2='#52514e', grid='#e3e2de',
          s1='#2a78d6', s2='#eb6834', s3='#1baf7a')
def style(ax):
    ax.set_facecolor(S['surface']); ax.grid(True, color=S['grid'], lw=.7, zorder=0)
    ax.set_axisbelow(True)
    for sp in ('top','right'): ax.spines[sp].set_visible(False)
    for sp in ('left','bottom'): ax.spines[sp].set_color(S['grid']); ax.spines[sp].set_linewidth(.9)
    ax.tick_params(colors=S['ink2'], labelsize=8)
plt.rcParams.update({'font.size':9,'axes.titlesize':10,'axes.labelsize':9})

# ---------------------------------------------------------------- F1: layer curve
r=list(csv.DictReader(open('layer_sweep/v1/p3_layer_curves.tsv'),delimiter='\t'))
L=sorted({int(x['layer']) for x in r})
g=lambda fs:[float([z for z in r if int(z['layer'])==l and z['feature_set']==fs][0]['auroc_P3']) for l in L]
comp,esm,both=g('comp'),g('esm'),g('comp_plus_esm')
fig,ax=plt.subplots(figsize=(7.0,3.5),dpi=300); fig.patch.set_facecolor(S['surface']); style(ax)
ax.plot(L,comp,color=S['s3'],lw=2,ls=(0,(5,3)),zorder=3,label='composition only (41 features)')
ax.plot(L,esm,color=S['s1'],lw=1.6,marker='o',ms=3.2,zorder=4,
        markeredgecolor=S['surface'],markeredgewidth=.8,label='ESM only (2,560)')
ax.plot(L,both,color=S['s2'],lw=2.2,marker='s',ms=3.4,zorder=5,
        markeredgecolor=S['surface'],markeredgewidth=.8,
        label='composition + ESM (2,601) — selection curve')
for x in (8,26): ax.axvline(x,color=S['ink2'],lw=.9,ls=':',zorder=1,alpha=.7)
i8,i26=L.index(8),L.index(26)
ax.annotate(f'layer 8\n{both[i8]:.4f}',(8,both[i8]),textcoords='offset points',xytext=(-2,-30),
            ha='center',fontsize=7.5,color=S['ink2'])
ax.annotate(f'layer 26 selected\n{both[i26]:.4f}',(26,both[i26]),textcoords='offset points',
            xytext=(0,12),ha='center',fontsize=7.5,color=S['ink'],fontweight='bold')
ax.set_xlabel('ESM-2 encoder layer   (hidden_states[i] = output of block i)',color=S['ink2'])
ax.set_ylabel('AUROC on P3  (n = 2,548)',color=S['ink2'])
ax.set_xticks([1,4,8,12,16,20,24,26,29,33]); ax.set_xlim(0,34); ax.set_ylim(.940,.9925)
lg=ax.legend(loc='lower right',frameon=False,fontsize=7.5)
for t in lg.get_texts(): t.set_color(S['ink2'])
fig.tight_layout(); fig.savefig(f'{OUT}/f1_layer_curve.png',facecolor=S['surface']); plt.close(fig)
print('F1 ok | comp flat at',round(comp[0],6),'| L26',round(both[i26],6),'| L8',round(both[i8],6))

# ---------------------------------------------------------------- F2: P4 comparisons
p4=json.load(open('layer_sweep/v1/p4_results.json'))
c=[x for x in p4['comparisons'] if x['subset']=='P4_all' and x['comparison']=='layer26_esm_vs_comp'][0]
fam=json.load(open('family_novelty_followup/subgroup_results.json'))['difference']
chg=json.load(open('charge_band_followup/charge_band_results.json'))['primary']
hyd=json.load(open('hydrophobicity_band_followup/hydro_band_results.json'))['primary']
rows=[('All P4 proteins\n2,366 (571 / 1,795) · 1,547 clusters',c['delta_auroc'],c['auroc_ci95_lo'],c['auroc_ci95_hi']),
      ('Unseen recorded Pfam family\n443 (183 / 260) · 234 clusters',fam['delta_auroc'],fam['auroc_ci95_lo'],fam['auroc_ci95_hi']),
      ('Pairs in the same charge band\n218,976 pairs · 1,547 clusters',chg['delta_auroc_esm_minus_composition'],chg['ci95_lo'],chg['ci95_hi']),
      ('Pairs in the same hydrophobicity band\n90,815 pairs · 1,546 clusters',hyd['delta_auroc_esm_minus_composition'],hyd['ci95_lo'],hyd['ci95_hi'])]
fig,ax=plt.subplots(figsize=(7.0,2.9),dpi=300); fig.patch.set_facecolor(S['surface']); style(ax)
y=np.arange(len(rows))[::-1]
ax.axvline(0,color=S['ink2'],lw=1.1,zorder=2)
for i,(lbl,d,lo,hi) in zip(y,rows):
    col=S['s1'] if lo>0 else S['s2']
    ax.plot([lo,hi],[i,i],color=col,lw=2.4,solid_capstyle='round',zorder=3)
    ax.plot([d],[i],'o',ms=7,color=col,markeredgecolor=S['surface'],markeredgewidth=1.2,zorder=4)
    ax.annotate(f'{d:+.4f}',(d,i),textcoords='offset points',xytext=(0,9),ha='center',
                fontsize=7.5,color=S['ink'],fontweight='bold')
ax.set_yticks(y); ax.set_yticklabels([r[0] for r in rows],fontsize=7.4,color=S['ink2'])
ax.set_xlabel('Δ AUROC  (layer-26 ESM-only  −  composition baseline)',color=S['ink2'])
ax.set_xlim(-0.02,0.11); ax.set_ylim(-0.6,len(rows)-0.3)
ax.annotate('interval crosses zero',(fam['auroc_ci95_lo'],y[1]),textcoords='offset points',
            xytext=(-4,-13),ha='left',fontsize=7,color=S['s2'],style='italic')
fig.tight_layout(); fig.savefig(f'{OUT}/f2_p4_comparisons.png',facecolor=S['surface']); plt.close(fig)
print('F2 ok |',[(r[0].split(chr(10))[0],round(r[1],6)) for r in rows])

# ---------------------------------------------------------------- F3: feature 1446, four groups
sig={(x['feature_id'],x['group']):x for x in csv.DictReader(
      open('candidate_discovery/v1/activation_by_signal_status.tsv'),delimiter='\t')}
cls={(x['feature_id'],x['group']):x for x in csv.DictReader(
      open('candidate_discovery/v1/activation_by_class.tsv'),delimiter='\t')}
F='1446'
ann_all=sig[(F,'annotated_present')]; nota=sig[(F,'not_annotated')]
cy_ann=sig[(F,'Cytoplasm_and_annotated_present')]; ex_not=sig[(F,'Extracellular_and_not_annotated')]
n_ann,n_not=int(ann_all['n']),int(nota['n'])
n_cy_ann,n_ex_not=int(cy_ann['n']),int(ex_not['n'])
n_ex_ann=n_ann-n_cy_ann; n_cy_not=n_not-n_ex_not
m_ex_ann=(n_ann*float(ann_all['mean'])-n_cy_ann*float(cy_ann['mean']))/n_ex_ann
m_cy_not=(n_not*float(nota['mean'])-n_ex_not*float(ex_not['mean']))/n_cy_not
p_ex_ann=(n_ann*float(ann_all['prevalence_nterm'])-n_cy_ann*float(cy_ann['prevalence_nterm']))/n_ex_ann
p_cy_not=(n_not*float(nota['prevalence_nterm'])-n_ex_not*float(ex_not['prevalence_nterm']))/n_cy_not
G=[('Extracellular\nannotated',n_ex_ann,m_ex_ann,p_ex_ann,S['s1']),
   ('Cytoplasmic\nannotated',n_cy_ann,float(cy_ann['mean']),float(cy_ann['prevalence_nterm']),S['s2']),
   ('Extracellular\nnot annot.',n_ex_not,float(ex_not['mean']),float(ex_not['prevalence_nterm']),S['s1']),
   ('Cytoplasmic\nnot annot.',n_cy_not,m_cy_not,p_cy_not,S['s2'])]
tot=sum(g[1] for g in G)
fig,axs=plt.subplots(1,2,figsize=(7.0,3.1),dpi=300)
fig.patch.set_facecolor(S['surface'])
ax=axs[0]; style(ax)
x=np.arange(4)
bars=ax.bar(x,[g[2] for g in G],color=[g[4] for g in G],width=.62,zorder=3,
            edgecolor=S['surface'],linewidth=1.4)
for xi,g in zip(x,G):
    ax.annotate(f'{g[2]:.2f}',(xi,g[2]),textcoords='offset points',xytext=(0,4),ha='center',
                fontsize=7.5,color=S['ink'],fontweight='bold')
ax.set_xticks(x)
ax.set_xticklabels([f'{g[0]}\nn = {g[1]:,}' for g in G],fontsize=6.4,color=S['ink2'])
ax.set_ylabel('mean activation of feature 1446\n(first 30 residues)',color=S['ink2'],fontsize=8)
ax.set_ylim(0,5.4)
ax=axs[1]; style(ax)
ax.bar(x,[100*g[3] for g in G],color=[g[4] for g in G],width=.62,zorder=3,
       edgecolor=S['surface'],linewidth=1.4)
for xi,g in zip(x,G):
    ax.annotate(f'{100*g[3]:.1f}%',(xi,100*g[3]),textcoords='offset points',xytext=(0,4),
                ha='center',fontsize=7.5,color=S['ink'],fontweight='bold')
ax.set_xticks(x)
ax.set_xticklabels([f'{g[0]}\nn = {g[1]:,}' for g in G],fontsize=6.4,color=S['ink2'])
ax.set_ylabel('% of proteins with any activity',color=S['ink2'],fontsize=8); ax.set_ylim(0,112)
fig.tight_layout(); fig.savefig(f'{OUT}/f3_feature1446_groups.png',facecolor=S['surface'],bbox_inches='tight'); plt.close(fig)
print(f'F3 ok | groups {[(g[0].replace(chr(10)," "),g[1],round(g[2],3),round(g[3],4)) for g in G]} | total {tot}')
json.dump({'groups':[{'group':g[0].replace('\n',' '),'n':g[1],'mean':g[2],'prevalence':g[3]} for g in G],
           'total':tot}, open(f'{OUT}/f3_values.json','w'),indent=1)

# ---------------------------------------------------------------- F4: 24-protein ablation
pp=list(csv.DictReader(open('probe_ablation_pilot/pilot_per_protein.tsv'),delimiter='\t'))
pp.sort(key=lambda z:(int(z['n_positions_active_1446'])==0, -float(z['paired_candidate_minus_control'])))
fig,ax=plt.subplots(figsize=(7.0,3.9),dpi=300); fig.patch.set_facecolor(S['surface']); style(ax)
xs=np.arange(len(pp))
vals=[float(z['paired_candidate_minus_control']) for z in pp]
abl=[float(z['mean_d_ablate']) for z in pp]
act=[int(z['n_positions_active_1446'])>0 for z in pp]
cols=[(S['s1'] if z['localisation_class']=='Extracellular' else S['s2']) for z in pp]
ax.axhline(0,color=S['ink'],lw=1.2,zorder=2)
ax.bar(xs,vals,color=[c if a else '#c9c8c3' for c,a in zip(cols,act)],width=.68,zorder=3,
       edgecolor=S['surface'],linewidth=.8)
ax.plot(xs,abl,'D',ms=3.4,color=S['ink2'],zorder=5,markeredgecolor=S['surface'],
        markeredgewidth=.6,label='ablation − clean')
n_act=sum(act)
ax.axvline(n_act-0.5,color=S['ink2'],lw=.9,ls='--',alpha=.8,zorder=2)
ax.annotate(f'{n_act} proteins with 1446 active at a tested position',(n_act/2-0.5,0.0665),
            ha='center',fontsize=7,color=S['ink2'])
ax.annotate(f'{len(pp)-n_act} inactive at all three positions:\nexact zero edit, exact zero change',
            ((n_act+len(pp))/2-0.5,0.0575),ha='center',fontsize=7,color=S['ink2'])
ax.set_xticks(xs); ax.set_xticklabels([z['ACC'] for z in pp],rotation=90,fontsize=5.6,color=S['ink2'])
ax.set_ylabel('Δ extracellular log-odds\n(mean over positions 10, 20, 30)',color=S['ink2'],fontsize=8)
ax.set_ylim(-0.02,0.098); ax.set_xlim(-0.8,len(pp)-0.2)
from matplotlib.patches import Patch
h=[Patch(facecolor=S['s1'],label='Extracellular  (bar: ablation − random control)'),
   Patch(facecolor=S['s2'],label='Cytoplasmic  (bar: ablation − random control)'),
   Patch(facecolor='#c9c8c3',label='1446 inactive at all three positions'),
   plt.Line2D([],[],marker='D',ls='',ms=4,color=S['ink2'],label='ablation − clean')]
lg=ax.legend(handles=h,loc='upper center',bbox_to_anchor=(0.5,-0.33),ncol=2,
             frameon=False,fontsize=7)
for t in lg.get_texts(): t.set_color(S['ink2'])
fig.tight_layout(); fig.savefig(f'{OUT}/f4_ablation_24.png',bbox_inches='tight',facecolor=S['surface']); plt.close(fig)
print(f'F4 ok | n={len(pp)} active={n_act} median_paired={np.median(vals):.6f} mean_abl={np.mean(abl):.6f}')

# ---------------------------------------------------------------- F5: masked-token panel
mt=list(csv.DictReader(open('layer_sweep/v1/intervention_cluster_sensitivity.tsv'),delimiter='\t'))
fig,ax=plt.subplots(figsize=(7.0,2.2),dpi=300); fig.patch.set_facecolor(S['surface']); style(ax)
y=np.arange(len(mt))[::-1]
ax.axvline(0,color=S['ink'],lw=1.1,zorder=2)
for i,m in zip(y,mt):
    lo,hi=float(m['cluster_ci9875_lo']),float(m['cluster_ci9875_hi']); est=float(m['estimate'])
    ax.plot([lo*1e5,hi*1e5],[i,i],color=S['s2'],lw=2.4,solid_capstyle='round',zorder=3)
    ax.plot([est*1e5],[i],'o',ms=6.5,color=S['s2'],markeredgecolor=S['surface'],markeredgewidth=1.2,zorder=4)
ax.set_yticks(y)
ax.set_yticklabels([f"candidate {m['candidate']}, ×{float(m['multiplier']):g}\n"
                    f"n = {m['n_proteins']} · {m['n_clusters_id30']} clusters" for m in mt],
                   fontsize=7,color=S['ink2'])
ax.set_xlabel('paired KL difference  (×10⁻⁵ nats)   ·   98.75 % cluster-bootstrap interval',color=S['ink2'])
fig.tight_layout(); fig.savefig(f'{OUT}/f5_masked_token.png',facecolor=S['surface']); plt.close(fig)
print('F5 ok | adjusted intervals excluding zero:',sum(int(m['cluster_ci9875_excludes_zero']) for m in mt),'of',len(mt))
