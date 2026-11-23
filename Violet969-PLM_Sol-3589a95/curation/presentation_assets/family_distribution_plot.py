#!/usr/bin/env python3
"""Named Pfam family and clan panels for the frozen 13,880-protein operational cohort.

Run from the repository root:
  python Violet969-PLM_Sol-3589a95/curation/presentation_assets/family_distribution_plot.py
No annotation/model/sequence/test execution or dataset writes.
"""
import argparse
from fractions import Fraction
import hashlib
import json
import os
from pathlib import Path

os.environ.setdefault('MPLCONFIGDIR', '/tmp/uesolds-matplotlib')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = Path(__file__).resolve().parent
GROUPS = [(1, 'Soluble', 7996), (0, 'Insoluble', 5884)]
NAVY, MUTED, PURPLE, TRACK = '#18243B', '#566176', '#8060B5', '#ECE7F3'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_validated():
    names = ['results/cohort_manifest.tsv', 'pfam_annotation/searched_manifest.tsv',
             'pfam_annotation/annotated_cohort.tsv', 'pfam_annotation/domains.tsv',
             'pfam_annotation/memberships.tsv', 'pfam_annotation/search_run.json']
    paths = [ROOT/name for name in names]
    hashes = dict(zip(names, map(sha, paths)))
    read = lambda p: pd.read_csv(p, sep='\t', keep_default_na=False, low_memory=False)
    original = read(paths[0])
    d = original.loc[original.clean & original.cached].sort_values('raw_uid').reset_index(drop=True)
    searched, annotated, domains, saved_memberships = map(read, paths[1:5])
    run = json.loads(paths[5].read_text())
    assert run['completed'] and run['domains_sha256'] == hashes['pfam_annotation/domains.tsv']
    identity = ['raw_uid', 'legacy_uid', 'original_id', 'split', 'label', 'sha256', 'length', 'cluster']
    assert len(d) == 13880 and d.raw_uid.is_unique
    for other in [searched, annotated]:
        assert len(other) == 13880 and other.raw_uid.is_unique
        pd.testing.assert_frame_equal(d[identity], other.sort_values('raw_uid').reset_index(drop=True)[identity])
    assert d.split.value_counts().to_dict() == {'train': 9880, 'validation': 4000}
    assert d.groupby(['split', 'label']).size().to_dict() == {
        ('train', 0): 3884, ('train', 1): 5996, ('validation', 0): 2000, ('validation', 1): 2000}
    assert d.label.value_counts().to_dict() == {1: 7996, 0: 5884}
    assert set(annotated.pfam_release.astype(str)) == {'37.0'}
    assert set(domains.raw_uid) <= set(d.raw_uid)
    assert (domains.sequence_score >= domains.sequence_ga).all()
    assert (domains.domain_score >= domains.domain_ga).all()
    # Repeated domains count once per protein–Pfam pair. Different families remain
    # distinct; neither clans nor MMseqs clusters replace Pfam accessions.
    pairs = domains[['raw_uid', 'pfam']].drop_duplicates()
    assert set(map(tuple, pairs.values)) == set(map(tuple, saved_memberships[['raw_uid', 'pfam']].drop_duplicates().values))
    membership = pairs.merge(d[['raw_uid', 'label', 'split']], on='raw_uid', validate='many_to_one')
    families_by_protein = pairs.groupby('raw_uid').pfam.agg(set).to_dict()
    ann = annotated.set_index('raw_uid')
    for r in d.itertuples():
        families = families_by_protein.get(r.raw_uid, set())
        assert families == set(filter(None, ann.loc[r.raw_uid, 'pfam_accessions'].split(';')))
        assert bool(families) == (ann.loc[r.raw_uid, 'pfam_status'] == 'assigned')
    return d, membership, families_by_protein, hashes


# Short biological display names; exact release descriptions remain in statistics.
ALIASES = {
 'PF00069': 'Protein kinase', 'PF13489': 'Methyltransferase',
 'PF00076': 'RNA recognition motif', 'PF04055': 'Radical SAM',
 'PF13416': 'Extracellular solute-binding', 'PF00583': 'GNAT acetyltransferase',
 'PF03466': 'LysR substrate-binding',
 'CL0063': 'FAD/NAD(P)-binding Rossmann fold', 'CL0123': 'Helix-turn-helix',
 'CL0023': 'P-loop NTP hydrolase', 'CL0036': 'Phosphate-binding TIM barrel',
 'CL0016': 'Protein kinase', 'CL0020': 'Tetratricopeptide repeat',
 'CL0177': 'Periplasmic binding protein',
 'other_family': 'Other families', 'other_clan': 'Other clans',
 'no_clan': 'Representative family has no clan', 'no_hit': 'No Pfam assignment'}
PALETTE = ['#829FAF', '#B69FBC', '#8FB3A6', '#C8AE85', '#99A2C6', '#C28F89', '#ABC08B']


def union_length(intervals):
    """Union length of 1-based inclusive aligned intervals, without double counting."""
    total, end = 0, 0
    for a, b in sorted(intervals):
        assert 1 <= a <= b
        total += max(0, b-max(a, end+1)+1)
        end = max(end, b)
    return total


def assignments(d, domains, ref):
    clan_map = ref.set_index('pfam').clan.to_dict()
    rows = []
    for uid, hits in domains.groupby('raw_uid', sort=True):
        candidates = []
        for fam, h in hits.groupby('pfam', sort=True):
            candidates.append((fam, union_length(zip(h.ali_from, h.ali_to)),
                               int((h.ali_to-h.ali_from+1).max())))
        ranked = sorted(candidates, key=lambda v: (-v[1], -v[2], v[0]))
        primary = ranked[0][0]
        score = hits.sort_values(['domain_score','ali_from','pfam'], ascending=[False,True,True]).iloc[0].pfam
        first = hits.sort_values(['ali_from','ali_to','pfam']).iloc[0].pfam
        rows.append(dict(raw_uid=uid, representative_family=primary, representative_clan=clan_map[primary] or 'no_clan',
                         representative_coverage=ranked[0][1], family_count=len(ranked), domain_hit_count=len(hits),
                         coverage_tie=sum(v[1] == ranked[0][1] for v in ranked)>1,
                         score_family=score, score_clan=clan_map[score] or 'no_clan',
                         first_family=first, first_clan=clan_map[first] or 'no_clan',
                         clan_count=len({clan_map[v[0]] for v in ranked if clan_map[v[0]]}),
                         no_clan_despite_other_clan=not clan_map[primary] and any(clan_map[v[0]] for v in ranked)))
    a = d.merge(pd.DataFrame(rows), on='raw_uid', how='left', validate='one_to_one')
    for col in ['representative_family','representative_clan','score_family','score_clan','first_family','first_clan']:
        a[col] = a[col].fillna('no_hit')
    for col in ['family_count','domain_hit_count','clan_count','representative_coverage']:
        a[col] = a[col].fillna(0).astype(int)
    for col in ['coverage_tie','no_clan_despite_other_clan']:
        a[col] = a[col].fillna(False).astype(bool)
    return a


def main(out):
    out.mkdir(parents=True, exist_ok=True)
    d, membership, protein_families, hashes = load_validated()
    read = lambda p: pd.read_csv(p, sep='\t', keep_default_na=False)
    domains = read(ROOT/'pfam_annotation/domains.tsv')
    ref = pd.read_csv(ROOT/'pfam_annotation/Pfam-A.clans.tsv', sep='\t', header=None,
                      names=['pfam','clan','clan_name','family_name','description'], keep_default_na=False)
    clans = read(OUT/'pfam37_clan_names.tsv')
    flagged = read(ROOT/'pfam_annotation/domains_flagged.tsv')
    for p in [ROOT/'pfam_annotation/Pfam-A.clans.tsv', ROOT/'pfam_annotation/domains_flagged.tsv', OUT/'pfam37_clan_names.tsv']:
        hashes[str(p.relative_to(ROOT))] = sha(p)
    assert ref.pfam.is_unique and clans.clan.is_unique
    assert set(ref.clan)-{''} <= set(clans.clan)
    assert set(domains.pfam) <= set(ref.pfam)
    lengths = d.set_index('raw_uid').length
    assert (domains.ali_from >= 1).all() and (domains.ali_to >= domains.ali_from).all()
    assert (domains.ali_to <= domains.raw_uid.map(lengths)).all()
    assert union_length([(1,5),(3,7),(10,10)]) == 8
    a = assignments(d, domains, ref)
    # The representative rule is independent of input row order and all labels.
    reranked = assignments(d, domains.sample(frac=1, random_state=42), ref)
    pd.testing.assert_frame_equal(a, reranked)
    mixed = set(membership.groupby('pfam').label.nunique().loc[lambda x:x==2].index)
    assert len(mixed) == 1660 and membership.pfam.nunique() == 6043
    a['mixed_family'] = a.raw_uid.map(lambda u: bool(protein_families.get(u,set()) & mixed))
    a['beyond_1022'] = a.length > 1022
    a['competing_overlap'] = a.raw_uid.isin(flagged.loc[flagged.clan_competing_overlap, 'raw_uid'])
    family_desc = ref.set_index('pfam').description.to_dict()
    clan_desc = clans.set_index('clan').description.to_dict()
    definitions = {**family_desc, **clan_desc,
        'other_family': 'Representative Pfam is outside the seven displayed families.',
        'other_clan': 'Representative-family clan is outside the seven displayed clans.',
        'no_clan': 'Representative Pfam has no clan in Pfam 37.0; another hit may have a clan.',
        'no_hit': 'No GA-passing Pfam 37.0 assignment.',
        'mixed_family': 'At least one original Pfam family has both labels in the operational cohort.',
        'beyond_1022': 'Original full sequence length exceeds 1,022 residues.'}
    top, categories = {}, {}
    for level in ['family','clan']:
        counts = a[f'representative_{level}'].value_counts().to_dict()
        top[level] = sorted((k for k in counts if k not in ['no_hit','no_clan']), key=lambda k:(-counts[k],k))[:7]
        categories[level] = top[level]+[f'other_{level}']+(['no_clan'] if level=='clan' else [])+['no_hit']
        for method, prefix in [('primary','representative'),('highest_score','score'),('first_domain','first')]:
            a[f'{method}_{level}_category'] = a[f'{prefix}_{level}'].map(
                lambda k:k if k in top[level] or k in ['no_hit','no_clan'] else f'other_{level}')
    rows = []
    def add(metric, category, count, n, group, label, level='', method='primary', description=''):
        pct = 100*int(count)/int(n) if n else None
        rows.append(dict(metric=metric, granularity=level, assignment_method=method, category=category,
                         name=ALIASES.get(category, definitions.get(category,category)),
                         definition=description or definitions.get(category,''), group=group, label=label,
                         count=int(count), denominator=n or '', percent=pct,
                         exact_percent_fraction=str(Fraction(100*int(count),int(n))) if n else '',
                         displayed_percent=f'{pct:.1f}%' if n else ''))
    groups = GROUPS+[('all','Combined',13880)]
    for label, group, n in groups:
        g = a if label=='all' else a[a.label==label]
        assert len(g)==n
        for level in categories:
            for method in ['primary','highest_score','first_domain']:
                counts = g[f'{method}_{level}_category'].value_counts()
                assert sum(counts)==n and set(counts.index)<=set(categories[level])
                for cat in categories[level]:
                    add('stacked_distribution',cat,counts.get(cat,0),n,group,label,level,method)
            for method, prefix in [('highest_score','score'),('first_domain','first')]:
                add('sensitivity_changed_assignment',level,(g[f'{prefix}_{level}']!=g[f'representative_{level}']).sum(),n,group,label,level,method)
                add('sensitivity_changed_plotted_category',level,(g[f'{method}_{level}_category']!=g[f'primary_{level}_category']).sum(),n,group,label,level,method)
        for metric in ['mixed_family','beyond_1022']:
            add('proportion_bar',metric,g[metric].sum(),n,group,label)
        for metric, mask in [('multiple_distinct_families',g.family_count>1),('multiple_domain_hits',g.domain_hit_count>1),
                             ('multiple_clans',g.clan_count>1),('coverage_tie',g.coverage_tie),
                             ('competing_clan_overlap',g.competing_overlap),('no_clan_despite_other_clan',g.no_clan_despite_other_clan)]:
            add('ambiguity',metric,mask.sum(),n,group,label)
    # All-hit membership prevalence is an examination of families, not a partition.
    full = membership.merge(ref[['pfam','clan']], on='pfam', validate='many_to_one')
    for level,col in [('family','pfam'),('clan','clan')]:
        pairs = full[full[col]!=''][['raw_uid','label',col]].drop_duplicates()
        for label,group,n in groups:
            selected=pairs if label=='all' else pairs[pairs.label==label]
            for cat,count in selected[col].value_counts().items():
                add('all_hit_membership',cat,count,n,group,label,level,'all_hits')
    add('distinct_families','mixed_label_families',len(mixed),None,'Combined','all','family','all_hits')
    stats = pd.DataFrame(rows)
    stats.to_csv(out/'family_distribution_statistics.tsv',sep='\t',index=False,float_format='%.17g',na_rep='')
    cols=['raw_uid','split','label','sha256','length','cluster','representative_family','representative_clan',
          'representative_coverage','family_count','domain_hit_count','clan_count','coverage_tie','competing_overlap',
          'no_clan_despite_other_clan','score_family','score_clan','first_family','first_clan',
          'primary_family_category','primary_clan_category','mixed_family','beyond_1022']
    a[cols].to_csv(out/'family_distribution_assignments.tsv',sep='\t',index=False)
    for level in categories:
        plot_named(a,level,categories[level],out)
    write_named_summary(a,stats,top,categories,hashes,out)
    assert {name:sha(ROOT/name) for name in hashes} == hashes
    print(stats[(stats.assignment_method=='primary') & stats.metric.isin(['stacked_distribution','proportion_bar']) & (stats.label!='all')][
        ['granularity','group','category','count','denominator','displayed_percent']].to_string(index=False))
    print('All frozen cohort, GA, coordinate, partition, order-invariance and input-immutability checks passed.')


def plot_named(a,level,categories,out):
    plt.rcParams.update({'font.family':'sans-serif','font.sans-serif':['Arial','DejaVu Sans'],
                        'svg.fonttype':'none','svg.hashsalt':'uesolds-named-pfam-v1','text.color':NAVY})
    fig=plt.figure(figsize=(6.6,6.5),facecolor='white')
    ax=fig.add_axes([0,0,1,1]);ax.set_xlim(0,660);ax.set_ylim(650,0);ax.axis('off')
    colors=dict(zip(categories[:7],PALETTE))
    colors.update({f'other_{level}':'#C5CBD3','no_clan':'#D8CDBB','no_hit':'#EAECEF'})
    x,width,height=24,520,17
    for (label,group,n),y in zip(GROUPS,[22,244]):
        g=a[a.label==label]
        ax.text(x,y,f'{group} ({n:,})',fontsize=14,weight='bold',va='top')
        bar_y=y+53
        title='Protein family distribution' if level=='family' else 'Protein superfamily distribution (Pfam clans)'
        ax.text(x,bar_y-10,title,fontsize=10.5,va='bottom')
        counts=g[f'primary_{level}_category'].value_counts();offset=0
        for cat in categories:
            segment=width*counts.get(cat,0)/n
            ax.add_patch(Rectangle((x+offset,bar_y),segment,height,facecolor=colors[cat],linewidth=0))
            offset+=segment
        assert abs(offset-width)<1e-8
        for metric,title,yy in [('mixed_family','In a mixed-label family',y+108),('beyond_1022','Beyond 1,022 residues',y+163)]:
            pct=100*g[metric].sum()/n
            ax.text(x,yy-10,title,fontsize=10.5,va='bottom')
            ax.add_patch(Rectangle((x,yy),width,height,facecolor=TRACK,linewidth=0))
            ax.add_patch(Rectangle((x,yy),width*pct/100,height,facecolor=PURPLE,linewidth=0))
            ax.text(x+width+16,yy+height/2,f'{pct:.1f}%',fontsize=12,weight='bold',color=PURPLE,va='center')
    # Column-major legend: named categories precede the residual categories.
    import textwrap
    for i,cat in enumerate(categories):
        lx=24+322*(i//5);ly=467+34*(i%5)
        ax.add_patch(Rectangle((lx,ly+2),10,10,facecolor=colors[cat],linewidth=0))
        label=ALIASES.get(cat,cat)
        if cat.startswith(('PF','CL')): label+=f' ({cat})'
        ax.text(lx+17,ly,'\n'.join(textwrap.wrap(label,43)),fontsize=8.5,color=MUTED,va='top',linespacing=1.15)
    stem=f'pfam_{level}_distribution'
    fig.savefig(out/f'{stem}.svg',facecolor='white',metadata={'Date':None,'Description':f'Pfam 37.0 {level} composition of UESolDS development cohort; one coverage-based representative family per protein.'})
    fig.savefig(out/f'{stem}.png',dpi=400,facecolor='white')
    plt.close(fig)


def write_named_summary(a,stats,top,categories,hashes,out):
    lines=['# Named Pfam composition of the UESolDS operational cohort','',
      'Verified: **13,880 proteins = 9,880 training + 4,000 validation**; soluble (label 1) **7,996**, insoluble (label 0) **5,884**.',
      'Both candidates use exactly this cohort and completed uniform Pfam 37.0 gathering-threshold annotations. No sequences, labels, splits, clusters or membership were changed; no independent-test files were accessed.','',
      '[Family PNG](pfam_family_distribution.png) · [Family SVG](pfam_family_distribution.svg) · [Clan PNG](pfam_clan_distribution.png) · [Clan SVG](pfam_clan_distribution.svg)',
      '[Statistics](family_distribution_statistics.tsv) · [Protein assignments](family_distribution_assignments.tsv) · [Script](family_distribution_plot.py)','',
      '## Protein-level assignment and category selection','',
      'For each protein, combine all aligned residue intervals assigned to each Pfam accession (1-based inclusive); rank families by total union coverage, then greatest single aligned-domain span, then ascending accession. The top family is the representative. Repeated/overlapping hits to the same family do not double-count residues. Coverage favors large domains or repeated domains: this is a reproducible display convention, not a biological claim that one domain defines the entire protein.',
      'Map that same representative family to its Pfam 37.0 clan for the broader candidate. A family without a clan remains explicitly “Representative family has no clan,” even if another domain has a clan. Clans are evolutionary groupings of domain families, not MMseqs clusters or complete protein architectures.',
      'For each granularity, select seven named categories with the most representative assignments in the **combined** cohort, breaking count ties by accession. This selection and the representative rule use no solubility labels, sources or performance. Every protein contributes exactly once to each stacked bar; residual categories are exhaustive and disjoint. “No Pfam assignment” means no GA-passing hit, not absence of biological domains.',
      'Display names are shortened biological descriptions, with accessions retained to distinguish similarly named families. Exact release descriptions are in the TSV.','']
    for level in ['family','clan']:
        lines += [f'## {level.title()} candidate: exact plotted counts','',
                  '| Category | Soluble (7,996) | Insoluble (5,884) | Combined (13,880) |','|---|---:|---:|---:|']
        for cat in categories[level]:
            cells=[]
            for label,_,_ in GROUPS+[('all','Combined',13880)]:
                r=stats[(stats.metric=='stacked_distribution')&(stats.assignment_method=='primary')&(stats.granularity==level)&(stats.category==cat)&(stats.label==label)].iloc[0]
                cells.append(f"{r['count']:,} ({r.percent:.6f}%)")
            lines.append('| '+ALIASES.get(cat,cat)+(f' ({cat})' if cat.startswith(('PF','CL')) else '')+' | '+' | '.join(cells)+' |')
        named=a[f'representative_{level}'].isin(top[level]).sum()
        lines += ['',f'The seven named {level} categories cover **{named:,}/13,880 ({100*named/13880:.2f}%)** of the cohort.','']
    lines += ['## Additional purple rows','', '| Quantity | Soluble | Insoluble |','|---|---:|---:|']
    for metric,title in [('mixed_family','In a mixed-label family'),('beyond_1022','Beyond 1,022 residues')]:
        cells=[]
        for label,_,n in GROUPS:
            count=int(a.loc[a.label==label,metric].sum());cells.append(f'{count:,}/{n:,} ({100*count/n:.6f}%; displayed {100*count/n:.1f}%)')
        lines.append('| '+title+' | '+' | '.join(cells)+' |')
    lines += ['', 'There are **6,043 distinct Pfam families**, including **1,660 mixed-label families** (at least one protein of each label). The mixed-family row uses **all original family assignments**, not only representatives. It is descriptive and uses both training and validation labels; it is not a held-out predictor or a minimum-support family contrast. Length is the original full-sequence length, and the threshold is strictly >1,022.','',
              '## All-hit prevalence before selecting representatives','',
              'These counts deduplicate protein–family or protein–clan pairs but allow proteins in multiple rows; they are **not** parts of a 100% bar. Top-seven lists below rank by combined prevalence independently of labels. Representative ranking differs because proteins can have multiple, sometimes overlapping alternatives; choosing a top-listed family whenever present would artificially privilege the displayed categories.','']
    for level in ['family','clan']:
        ranked=stats[(stats.metric=='all_hit_membership')&(stats.granularity==level)&(stats.label=='all')].sort_values(['count','category'],ascending=[False,True]).head(7)
        lines += [f'### Most prevalent {level} memberships','', '| Accession | Exact biological description | Proteins | Cohort percent |','|---|---|---:|---:|']
        for r in ranked.itertuples(): lines.append(f'| {r.category} | {r.definition} | {r.count:,} | {r.percent:.6f}% |')
        lines.append('')
    lines += ['## Multidomain and overlapping-assignment ambiguity','',
              'Multiple distinct Pfam assignments may describe distinct domains, repeated regions, or overlapping alternative profiles. They must not all be called separate biological domains. The existing same-clan competing-overlap flag is retained without resolving or rerunning annotations.','',
              '| Diagnostic | Soluble | Insoluble | Combined |','|---|---:|---:|---:|']
    for metric in stats.loc[stats.metric=='ambiguity','category'].unique():
        sub=stats[(stats.metric=='ambiguity')&(stats.category==metric)]
        lines.append('| '+metric+' | '+' | '.join(f"{r['count']:,} ({r.percent:.2f}%)" for _,r in sub.iterrows())+' |')
    lines += ['', 'Alternative label-independent rules: highest domain bit score (ties: earliest start, accession), or earliest domain (ties: earliest end, accession). Bit scores are not uniformly calibrated across Pfam models, so that alternative is only a sensitivity diagnostic. Keep the primary top-seven categories fixed when comparing plotted distributions.','',
              '| Alternative | Level | Changed representative, combined | Changed plotted category, combined | Largest category shift, soluble / insoluble |','|---|---|---|---|']
    for method,prefix in [('highest_score','score'),('first_domain','first')]:
        for level in ['family','clan']:
            changed=(a[f'{prefix}_{level}']!=a[f'representative_{level}']).sum()
            plotted=(a[f'{method}_{level}_category']!=a[f'primary_{level}_category']).sum()
            shifts=[]
            for label,_,n in GROUPS:
                g=a[a.label==label]
                p=g[f'primary_{level}_category'].value_counts();q=g[f'{method}_{level}_category'].value_counts()
                shifts.append(max(abs(q.get(c,0)-p.get(c,0))*100/n for c in categories[level]))
            lines.append(f'| {method} | {level} | {changed:,} ({100*changed/13880:.2f}%) | {plotted:,} ({100*plotted/13880:.2f}%) | {shifts[0]:.3f} / {shifts[1]:.3f} percentage points |')
    lines += ['', 'Exact alternative category counts, all-hit prevalence across all families/clans, and per-class ambiguity counts are saved in the statistics TSV. Per-protein representative and alternate assignments are saved separately.','',
              '## Interpretation and caption','',
              'The **clan candidate is the more interpretable of these two displays**, because it groups related domain families. However, neither top-seven display summarizes most of this biologically diverse cohort: residual categories dominate, particularly at family level. The figure should not imply that the named categories are the dominant biology of the entire dataset. A fuller family distribution would require a different visualization or many more categories.',
              'Suggested caption: **“Soluble and insoluble proteins span diverse Pfam families; broader clan groupings reveal their domain-family composition.”**',
              'Methods footnote: “Pfam 37.0; one coverage-based representative domain family per protein; seven categories selected by combined prevalence. Mixed-family membership uses all Pfam hits.”',
              'This is descriptive composition, not evidence for a causal solubility mechanism or ESM representation. Broad clan names do not establish precise function for every member.','',
              '## Reference metadata and reproduction','',
              'No Pfam searches were rerun. Existing Pfam-A.clans.tsv supplies family descriptions and clan membership. Full clan names were extracted into the small pfam37_clan_names.tsv reference from pinned [Pfam 37.0 Pfam-C.gz](https://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam37.0/Pfam-C.gz), #=GF AC/ID/DE fields (Latin-1 decoding). Download MD5 verified: `5646960c822e6bada8be93afc5b4aa09`; SHA-256: `a3bd7c5fa1d07929b7a600277281cb6088a9738d6b1220a4f556ad444e4815ce`. No reference download is needed to regenerate figures.',
              'From repository root:', '```sh','python Violet969-PLM_Sol-3589a95/curation/presentation_assets/family_distribution_plot.py','```',
              f'pandas {pd.__version__}; matplotlib {matplotlib.__version__}. PNGs have white backgrounds, 2640 × 2600 pixels, 400 dpi; SVG text remains editable. Bar proportions use exact integer ratios; only displayed labels round to one decimal. The TSV includes exact rational percentages.',
              'Checks include cohort counts; stable ID/hash/split/label/cluster alignment; completed annotation checksum; GA thresholds; valid coordinates; domain/membership consistency; exclusive/exhaustive partitions; order-invariant assignment; and unchanged input hashes. No independent-test reads, annotation reruns or model execution. The original DeepLoc slide is unavailable locally; the panels follow the requested aligned bars and navy/purple/muted style.','',
              'Input SHA-256:','']
    lines += [f'- `{name}`: `{digest}`' for name,digest in hashes.items()]
    (out/'family_distribution_summary.md').write_text('\n'.join(lines)+'\n')


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out',type=Path,default=OUT)
    main(parser.parse_args().out)
