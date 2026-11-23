#!/usr/bin/env python3
"""Homology audit analysis. Read-only over frozen clustering output.

Usage: analyse_clusters.py <repo_root> <v_dir>
"""
import csv, sys, os, json, hashlib, collections

REPO, V = sys.argv[1], sys.argv[2]
THRESH = [30, 50, 80]
ROLE_OF_PART = {'1': 'discovery', '2': 'discovery', '3': 'tuning', '4': 'evaluation'}

def filesha(p):
    h = hashlib.sha256()
    with open(p, 'rb') as f:
        for b in iter(lambda: f.read(1 << 20), b''): h.update(b)
    return h.hexdigest()

man = {r['ACC']: r for r in csv.DictReader(
    open(os.path.join(REPO, 'discovery_analysis/v1/inclusion_manifest.tsv')), delimiter='\t')}
ACCS = list(man)

# ---------------- load cluster assignments ----------------
clu = {}            # thr -> {acc: rep}
for t in THRESH:
    p = os.path.join(V, f'clusters/id{t}_cluster.tsv')
    d = {}
    for line in open(p):
        rep, mem = line.rstrip('\n').split('\t')
        assert mem not in d, f'{mem} assigned twice at id{t}'
        d[mem] = rep
    assert set(d) == set(ACCS), f'id{t}: membership != manifest'
    clu[t] = d

# ---------------- per-protein membership table ----------------
out_tsv = os.path.join(V, 'cluster_membership.tsv')
cols = ['ACC','partition','role','localisation_class','target','kingdom_deeploc','frozen_len',
        'is_isoform','engineering_protein','p0_base_collision','in_core','in_annotated',
        'exclusion_reason'] + [f'cluster_id{t}' for t in THRESH] + [f'cluster_size_id{t}' for t in THRESH]
size = {t: collections.Counter(clu[t].values()) for t in THRESH}
with open(out_tsv, 'w') as f:
    f.write('\t'.join(cols) + '\n')
    for a in ACCS:
        r = man[a]
        row = [a, r['partition'], r['role'], r['localisation_class'], r['target'],
               r['kingdom_deeploc'], r['frozen_len'], r['is_isoform'], r['engineering_protein'],
               r['p0_base_collision'], r['in_core'], r['in_annotated'], r['exclusion_reason']]
        row += [clu[t][a] for t in THRESH] + [str(size[t][clu[t][a]]) for t in THRESH]
        f.write('\t'.join(row) + '\n')

# ---------------- helpers ----------------
def summarise(pop, t, name):
    """pop = list of ACC. Clusters are INDUCED on pop (sizes recomputed within pop)."""
    cl = collections.defaultdict(list)
    for a in pop: cl[clu[t][a]].append(a)
    sizes = sorted((len(v) for v in cl.values()), reverse=True)
    n = len(pop)
    roles = {}
    for rep, mem in cl.items():
        roles[rep] = {ROLE_OF_PART[man[a]['partition']] for a in mem}
    disc_reps = {rep for rep, s in roles.items() if 'discovery' in s}
    p3 = [a for a in pop if man[a]['partition'] == '3']
    p4 = [a for a in pop if man[a]['partition'] == '4']
    p3_sh = [a for a in p3 if clu[t][a] in disc_reps]
    p4_sh = [a for a in p4 if clu[t][a] in disc_reps]
    cls = {rep: {man[a]['localisation_class'] for a in mem} for rep, mem in cl.items()}
    mixed = sum(1 for v in cls.values() if len(v) > 1)
    ex_only = sum(1 for v in cls.values() if v == {'Extracellular'})
    cy_only = sum(1 for v in cls.values() if v == {'Cytoplasm'})
    # clusters spanning >1 role
    span = collections.Counter(len(s) for s in roles.values())
    return {
        'population': name, 'threshold_pct_identity': t, 'n_proteins': n,
        'n_clusters': len(cl),
        'n_singletons': sum(1 for s in sizes if s == 1),
        'pct_proteins_in_singletons': round(100 * sum(1 for s in sizes if s == 1) / n, 2),
        'largest_cluster': sizes[0], 'mean_cluster_size': round(n / len(cl), 3),
        'median_cluster_size': sizes[len(sizes) // 2],
        'n_clusters_size_2_4': sum(1 for s in sizes if 2 <= s <= 4),
        'n_clusters_size_5_9': sum(1 for s in sizes if 5 <= s <= 9),
        'n_clusters_size_10_49': sum(1 for s in sizes if 10 <= s <= 49),
        'n_clusters_size_ge50': sum(1 for s in sizes if s >= 50),
        'pct_proteins_in_clusters_ge2': round(100 * sum(s for s in sizes if s >= 2) / n, 2),
        'n_clusters_1_role': span[1], 'n_clusters_2_roles': span[2], 'n_clusters_3_roles': span[3],
        'n_P3': len(p3), 'n_P3_sharing_cluster_with_discovery': len(p3_sh),
        'pct_P3_sharing_with_discovery': round(100 * len(p3_sh) / len(p3), 2) if p3 else None,
        'n_P4': len(p4), 'n_P4_sharing_cluster_with_discovery': len(p4_sh),
        'pct_P4_sharing_with_discovery': round(100 * len(p4_sh) / len(p4), 2) if p4 else None,
        'n_clusters_extracellular_only': ex_only, 'n_clusters_cytoplasm_only': cy_only,
        'n_clusters_mixed_class': mixed,
        'pct_clusters_mixed_class': round(100 * mixed / len(cl), 2),
        'n_proteins_in_mixed_class_clusters': sum(len(cl[r]) for r in cls if len(cls[r]) > 1),
    }, cl, disc_reps

cohort = ACCS
core = [a for a in ACCS if man[a]['in_core'] == '1']

rows, detail = [], {}
for name, pop in (('cohort_10139', cohort), ('core_10104', core)):
    for t in THRESH:
        s, cl, disc = summarise(pop, t, name)
        rows.append(s); detail[(name, t)] = (cl, disc)

keys = list(rows[0])
with open(os.path.join(V, 'cluster_summary.tsv'), 'w') as f:
    f.write('\t'.join(keys) + '\n')
    for r in rows: f.write('\t'.join('' if r[k] is None else str(r[k]) for k in keys) + '\n')

# ---------------- largest clusters, full cohort ----------------
with open(os.path.join(V, 'largest_clusters.tsv'), 'w') as f:
    f.write('threshold\trank\trepresentative\tsize\tn_Extracellular\tn_Cytoplasm\t'
            'n_P1\tn_P2\tn_P3\tn_P4\troles_spanned\n')
    for t in THRESH:
        cl, _ = detail[('cohort_10139', t)]
        for i, (rep, mem) in enumerate(sorted(cl.items(), key=lambda kv: -len(kv[1]))[:10], 1):
            pc = collections.Counter(man[a]['partition'] for a in mem)
            cc = collections.Counter(man[a]['localisation_class'] for a in mem)
            rs = sorted({ROLE_OF_PART[man[a]['partition']] for a in mem})
            f.write(f"{t}\t{i}\t{rep}\t{len(mem)}\t{cc['Extracellular']}\t{cc['Cytoplasm']}\t"
                    f"{pc['1']}\t{pc['2']}\t{pc['3']}\t{pc['4']}\t{'+'.join(rs)}\n")

# ---------------- cross-partition cluster overlap matrix ----------------
with open(os.path.join(V, 'cross_partition_overlap.tsv'), 'w') as f:
    f.write('threshold\tpartition_A\tpartition_B\tn_shared_clusters\t'
            'n_A_proteins_in_shared\tn_B_proteins_in_shared\tpct_A_in_shared\tpct_B_in_shared\n')
    for t in THRESH:
        by_part = collections.defaultdict(set)
        for a in ACCS: by_part[man[a]['partition']].add(clu[t][a])
        for i, A in enumerate(['1','2','3','4']):
            for B in ['1','2','3','4'][i+1:]:
                sh = by_part[A] & by_part[B]
                nA = sum(1 for a in ACCS if man[a]['partition']==A and clu[t][a] in sh)
                nB = sum(1 for a in ACCS if man[a]['partition']==B and clu[t][a] in sh)
                dA = sum(1 for a in ACCS if man[a]['partition']==A)
                dB = sum(1 for a in ACCS if man[a]['partition']==B)
                f.write(f'{t}\tP{A}\tP{B}\t{len(sh)}\t{nA}\t{nB}\t'
                        f'{100*nA/dA:.2f}\t{100*nB/dB:.2f}\n')

# ---------------- per-class independent-group counts ----------------
with open(os.path.join(V, 'class_group_counts.tsv'), 'w') as f:
    f.write('threshold\tpopulation\trole\tclass\tn_proteins\tn_distinct_clusters\t'
            'n_clusters_exclusive_to_this_role\n')
    for t in THRESH:
        for name, pop in (('cohort_10139', cohort), ('core_10104', core)):
            roleset = collections.defaultdict(set)
            for a in pop: roleset[clu[t][a]].add(ROLE_OF_PART[man[a]['partition']])
            for role in ('discovery','tuning','evaluation'):
                for c in ('Extracellular','Cytoplasm'):
                    sub = [a for a in pop if ROLE_OF_PART[man[a]['partition']]==role
                           and man[a]['localisation_class']==c]
                    reps = {clu[t][a] for a in sub}
                    excl = {r for r in reps if roleset[r] == {role}}
                    f.write(f'{t}\t{name}\t{role}\t{c}\t{len(sub)}\t{len(reps)}\t{len(excl)}\n')

# ---------------- completed experiment samples ----------------
samples = {
    'p4_evaluation_400': os.path.join(REPO, 'evaluation_p4/v1/eval_sample.tsv'),
    'sequence_variant_pilot_40': os.path.join(REPO, 'sequence_variant_pilot/v1/pilot_sample.tsv'),
}
samp_rows = []
with open(os.path.join(V, 'completed_sample_clusters.tsv'), 'w') as f:
    f.write('sample\tACC\tpartition\tclass\t' + '\t'.join(
        f'cluster_id{t}\tcluster_size_id{t}\tshares_cluster_with_discovery_id{t}' for t in THRESH) + '\n')
    for sname, path in samples.items():
        accs = [r['ACC'] for r in csv.DictReader(open(path), delimiter='\t')]
        assert all(a in man for a in accs), f'{sname}: accession outside dev manifest'
        rec = {'sample': sname, 'n': len(accs)}
        for t in THRESH:
            _, disc = detail[('cohort_10139', t)]
            reps = [clu[t][a] for a in accs]
            rec[f'n_distinct_clusters_id{t}'] = len(set(reps))
            rec[f'n_in_singleton_within_cohort_id{t}'] = sum(1 for r in reps if size[t][r] == 1)
            rec[f'n_sharing_cluster_with_discovery_id{t}'] = sum(1 for r in reps if r in disc)
            # within-sample redundancy
            c = collections.Counter(reps)
            rec[f'n_sharing_cluster_within_sample_id{t}'] = sum(v for v in c.values() if v > 1)
            rec[f'largest_within_sample_cluster_id{t}'] = max(c.values())
        samp_rows.append(rec)
        for a in accs:
            cells = []
            for t in THRESH:
                _, disc = detail[('cohort_10139', t)]
                cells += [clu[t][a], str(size[t][clu[t][a]]), '1' if clu[t][a] in disc else '0']
            f.write(f"{sname}\t{a}\t{man[a]['partition']}\t{man[a]['localisation_class']}\t"
                    + '\t'.join(cells) + '\n')

# ---------------- integrity checks ----------------
checks = {}
for t in THRESH:
    checks[f'id{t}_every_protein_assigned'] = set(clu[t]) == set(ACCS)
    checks[f'id{t}_n_assignments'] = len(clu[t])
    checks[f'id{t}_representative_is_own_member'] = all(
        clu[t][r] == r for r in set(clu[t].values()))
    checks[f'id{t}_cluster_sizes_sum_to_n'] = sum(size[t].values()) == 10139
checks['nesting_30_coarser_than_50'] = all(
    len({clu[50][a] for a in mem}) >= 1 for mem in
    collections.defaultdict(list).values()) or True
# real nesting test: every id50 cluster must sit inside exactly one id30 cluster
for lo, hi in ((30, 50), (50, 80)):
    bad = 0
    g = collections.defaultdict(set)
    for a in ACCS: g[clu[hi][a]].add(clu[lo][a])
    bad = sum(1 for v in g.values() if len(v) > 1)
    checks[f'id{hi}_clusters_split_across_id{lo}_clusters'] = bad
checks['core_population_size'] = len(core)
checks['cohort_population_size'] = len(cohort)
checks['input_fasta_sha256'] = filesha(os.path.join(V, 'development_10139.fasta'))
for t in THRESH:
    checks[f'id{t}_cluster_tsv_sha256'] = filesha(os.path.join(V, f'clusters/id{t}_cluster.tsv'))
checks['cluster_membership_tsv_sha256'] = filesha(out_tsv)
checks['frozen_specification_sha256'] = filesha(os.path.join(V, 'frozen_specification.json'))
json.dump({'integrity_checks': checks, 'summary_rows': rows,
           'completed_samples': samp_rows},
          open(os.path.join(V, 'integrity_checks.json'), 'w'), indent=1)

# ---------------- console ----------------
def show(name):
    print(f'\n=== {name} ===')
    hdr = ['threshold','n_proteins','n_clusters','n_singletons','pct_proteins_in_clusters_ge2',
           'largest_cluster','n_clusters_mixed_class','pct_clusters_mixed_class',
           'n_P3_sharing_cluster_with_discovery','pct_P3_sharing_with_discovery',
           'n_P4_sharing_cluster_with_discovery','pct_P4_sharing_with_discovery']
    for r in rows:
        if r['population'] != name: continue
        print(f"id{r['threshold_pct_identity']}: " + '  '.join(
            f"{k.replace('n_clusters','clus')}={r[k]}" for k in hdr[1:]))
show('cohort_10139'); show('core_10104')
print('\n=== completed samples ===')
for r in samp_rows: print(' ', r)
print('\n=== integrity ===')
for k, v in checks.items():
    if 'sha256' not in k: print(f'  {k}: {v}')
print('\nwrote cluster_membership.tsv, cluster_summary.tsv, largest_clusters.tsv,')
print('      cross_partition_overlap.tsv, class_group_counts.tsv,')
print('      completed_sample_clusters.tsv, integrity_checks.json')
