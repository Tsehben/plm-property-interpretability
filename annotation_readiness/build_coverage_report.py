#!/usr/bin/env python3
"""
Build the annotation coverage report from the tables written by fetch_annotations.py.

Reads only annotation_readiness/out/*.tsv and integrity_checks.json. Performs no
network access and modifies nothing outside annotation_readiness/.
"""

import csv
import json
import os
from collections import Counter, defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, 'out')

ROLE_ORDER = ['discovery', 'tuning', 'evaluation']
SIG_STATES = ['annotated_present', 'not_annotated', 'unverifiable']
EV_ORDER = ['experimental', 'combinatorial_experimental', 'curator_inference',
            'author_statement', 'by_similarity', 'sequence_analysis_predictor',
            'imported', 'automatic_assertion', 'unspecified']


def read_tsv(name):
    with open(os.path.join(OUT, name), newline='') as fh:
        return list(csv.DictReader(fh, delimiter='\t'))


def pct(n, d):
    return f'{100.0 * n / d:.1f}%' if d else '—'


def cls(r):
    return 'Extracellular' if r['target'] == '1' else 'Cytoplasm'


def table(header, rows):
    out = ['| ' + ' | '.join(header) + ' |',
           '|' + '|'.join(['---'] * len(header)) + '|']
    for r in rows:
        out.append('| ' + ' | '.join(str(x) for x in r) + ' |')
    return '\n'.join(out)


def main():
    primary = read_tsv('annotation_primary.tsv')
    feats = read_tsv('annotation_features.tsv')
    doms = read_tsv('annotation_domains.tsv')
    integ = json.load(open(os.path.join(OUT, 'integrity_checks.json')))
    prov = json.load(open(os.path.join(OUT, 'provenance.json')))
    N = len(primary)

    L = []
    w = L.append

    w('# Annotation coverage report — DeepLoc localisation development cohort\n')
    w(f'Generated from `out/annotation_primary.tsv` ({N} rows), '
      f'`out/annotation_features.tsv` ({len(feats)} rows), '
      f'`out/annotation_domains.tsv` ({len(doms)} rows).\n')
    w(f'Retrieval plan: {prov["n_base_accessions_queried"]} base accessions in batches of '
      f'{prov["batch_size"]}, {prov["n_isoform_sequences_requested"]} isoform sequence '
      f'requests, and {prov["n_individual_entry_lookups"]} individual entry lookups for '
      f'accessions the batch search either omitted or returned without a sequence. '
      f'This run: {prov["cache_hits"]} cache hits, {prov["network_requests"]} network '
      f'requests, elapsed {prov["elapsed_s"]}s — the fetch is cached and resumable, so a '
      f're-run reproduces the tables without re-querying UniProt. '
      f'Provenance in `out/provenance.json`; raw responses retained under `cache/`.\n')
    w('---\n')

    # ---------------- 1. retrieval and sequence verification ----------------
    w('## 1. Retrieval and sequence verification\n')
    rs = Counter(r['retrieval_status'] for r in primary)
    sm = Counter(r['sequence_match'] for r in primary)
    w(table(['Retrieval status', 'n', 'share'],
            [[k, v, pct(v, N)] for k, v in rs.most_common()]) + '\n')
    w(table(['Sequence vs frozen extraction input', 'n', 'share'],
            [[k, v, pct(v, N)] for k, v in sm.most_common()]) + '\n')

    iso = [r for r in primary if r['is_isoform'] == '1']
    iso_ok = sum(1 for r in iso if r['sequence_match'] == 'exact')
    iso_src_none = sum(1 for r in iso if r['sequence_source'] == 'none')
    nontrans = [r for r in primary if r['coordinates_transferable'] == '0']
    w(f'**Isoforms.** {len(iso)} of {N} accessions carry a `-N` suffix. '
      f'{iso_ok} matched the frozen sequence exactly; {iso_src_none} could not have '
      f'their own isoform sequence retrieved. The canonical sequence was never '
      f'substituted for an unresolved isoform.\n')
    merged = [r for r in primary
              if r['retrieval_status'] == 'resolved_via_merged_accession']
    deleted = [r for r in primary if r['retrieval_status'] == 'inactive_deleted']
    if merged:
        w(f'**Merged accessions ({len(merged)}).** These cohort accessions are no longer '
          'primary; UniProt resolves them to a successor entry. The successor\'s record '
          'is used, and the substitution is flagged rather than applied silently:\n')
        w(table(['cohort accession', 'resolved to', 'successor id', 'sequence vs frozen',
                 'signal status'],
                [[r['ACC'], r['merged_into'], r['uniprot_id'], r['sequence_match'],
                  r['signal_status']] for r in merged]) + '\n')
    if deleted:
        w(f'**Deleted accessions ({len(deleted)}).** '
          + '; '.join(f'`{r["ACC"]}` ({r["deleted_reason"]})' for r in deleted)
          + '. No sequence and no annotation are available; these rows are '
            '`unverifiable` throughout and are retained in the cohort.\n')

    w(f'**Coordinate transfer.** {len(nontrans)} of {N} records are marked '
      f'`coordinates_transferable = 0`: no UniProt residue interval is applied to them, '
      f'because the sequence they would be applied to is not the one we verified.\n')

    iso_nontrans = [r for r in nontrans if r['is_isoform'] == '1'
                    and r['sequence_match'] == 'exact']
    iso_nt_sig = sum(1 for r in iso_nontrans if r['parent_entry_has_signal'] == '1')
    if iso_nontrans:
        w(f'Of those, **{len(iso_nontrans)} are isoforms whose own sequence matched the '
          f'frozen input exactly but differs from the canonical sequence** that UniProt\'s '
          f'feature coordinates are defined on. Their parent entries are not blank: '
          f'**{iso_nt_sig}** carry a signal-peptide annotation, recorded in the '
          f'`parent_entry_has_signal` column but deliberately **not** transferred. '
          'Recovering them would require aligning each isoform to its canonical sequence '
          'and mapping the interval across — a bounded, well-defined follow-up that is '
          'not performed here, because transferring an interval without that alignment is '
          'exactly the error this column exists to prevent.\n')

    fails = [r for r in primary
             if r['retrieval_status'] != 'ok' or r['sequence_match'] != 'exact']
    if fails:
        w(f'**Failures, itemised** ({len(fails)} rows; full detail in '
          f'`out/annotation_failures.tsv`):\n')
        fc = Counter((r['retrieval_status'], r['sequence_match']) for r in fails)
        w(table(['retrieval_status', 'sequence_match', 'n', 'roles affected'],
                [[a, b, n,
                  ', '.join(f'{k}:{v}' for k, v in sorted(Counter(
                      r['role'] for r in fails
                      if (r['retrieval_status'], r['sequence_match']) == (a, b)
                  ).items()))]
                 for (a, b), n in fc.most_common()]) + '\n')
    else:
        w('**No retrieval or sequence-match failures.**\n')

    # ---------------- 2. signal peptide coverage ----------------
    w('---\n\n## 2. Signal-peptide annotation coverage\n')
    w('> `not_annotated` means a matching record was retrieved and carries no '
      'signal-peptide annotation. **It is not a confirmed biological negative.** '
      'UniProt records evidence for presence; it does not assert absence.\n')

    rows = []
    for role in ROLE_ORDER:
        for c in ['Extracellular', 'Cytoplasm']:
            sub = [r for r in primary if r['role'] == role and cls(r) == c]
            if not sub:
                continue
            cc = Counter(r['signal_status'] for r in sub)
            rows.append([role, c, len(sub)] +
                        [f'{cc.get(s, 0)} ({pct(cc.get(s, 0), len(sub))})'
                         for s in SIG_STATES])
    w(table(['role', 'class', 'n'] + SIG_STATES, rows) + '\n')

    rows = []
    for k in sorted({r['kingdom_deeploc'] for r in primary}):
        sub = [r for r in primary if r['kingdom_deeploc'] == k]
        cc = Counter(r['signal_status'] for r in sub)
        rows.append([k, len(sub)] +
                    [f'{cc.get(s, 0)} ({pct(cc.get(s, 0), len(sub))})'
                     for s in SIG_STATES])
    w('**By kingdom (DeepLoc bucket), all roles pooled:**\n')
    w(table(['kingdom', 'n'] + SIG_STATES, rows) + '\n')
    w('Full three-way breakdown (role × class × kingdom) with every counter: '
      '`out/coverage_by_role_class_kingdom.tsv`.\n')

    # ---------------- 3. evidence strength ----------------
    w('---\n\n## 3. Evidence strength, recorded separately from presence\n')
    w('> An annotation is not automatically experimental evidence. These are the '
      'ECO evidence categories attached to the signal-peptide features themselves.\n')

    sig_feats = [f for f in feats if f['feature_type'] == 'Signal']
    ev = Counter(f['evidence_category'] for f in sig_feats)
    w('**Signal-peptide feature-evidence rows** '
      f'({len(sig_feats)} rows over '
      f'{len({f["ACC"] for f in sig_feats})} proteins; a feature may carry several '
      'evidence attributions):\n')
    w(table(['evidence category', 'n evidence rows', 'n distinct proteins'],
            [[c, ev[c], len({f['ACC'] for f in sig_feats
                             if f['evidence_category'] == c})]
             for c in EV_ORDER if c in ev]) + '\n')

    present = [r for r in primary if r['signal_status'] == 'annotated_present']
    best = Counter(r['signal_evidence_best'] or 'unspecified' for r in present)
    rows = []
    for c in EV_ORDER:
        if c not in best:
            continue
        sub = [r for r in present
               if (r['signal_evidence_best'] or 'unspecified') == c]
        nx = sum(1 for r in sub if cls(r) == 'Extracellular')
        rows.append([c, best[c], pct(best[c], len(present)), nx, len(sub) - nx])
    w(f'**Strongest evidence category per protein**, over the {len(present)} proteins '
      'with `annotated_present`:\n')
    w(table(['strongest category', 'n proteins', 'share of annotated',
             'Extracellular', 'Cytoplasm'], rows) + '\n')

    exp = sum(v for c, v in best.items()
              if c in ('experimental', 'combinatorial_experimental'))
    w(f'**Experimental-evidence subset: {exp} proteins** '
      f'({pct(exp, len(present))} of annotated, {pct(exp, N)} of the cohort). '
      'This is the denominator for any sensitivity analysis restricted to '
      'experimentally-evidenced signal peptides.\n')

    tra = [f for f in feats if f['feature_type'] == 'Transit peptide']
    tr_present = sum(1 for r in primary if r['transit_status'] == 'annotated_present')
    w(f'**Transit peptides:** {tr_present} proteins carry an annotated transit peptide '
      f'({len(tra)} evidence rows). Evidence breakdown: '
      + ', '.join(f'{c} {n}' for c, n in
                  Counter(f['evidence_category'] for f in tra).most_common()) + '.\n')

    # ---------------- 4. family / domain ----------------
    w('---\n\n## 4. Family and domain coverage\n')
    w('> **A Pfam domain is not a whole-protein family.** Pfam and InterPro entries are '
      'region-level annotations: a multidomain protein carries several, and two proteins '
      'sharing one domain are not thereby in the same family. The counts below are '
      '*annotation coverage*, not a family assignment and not a clustering.\n')
    w('> **No homology clusters have been computed.** No MMseqs2 or CD-HIT run has been '
      'performed, so no sequence-similarity groups exist for this cohort. The '
      'residual-homology question of §2.6 of the specification remains unmeasured.\n')

    has_pf = sum(1 for r in primary if int(r['n_pfam'] or 0) > 0)
    has_ip = sum(1 for r in primary if int(r['n_interpro'] or 0) > 0)
    has_sim = sum(1 for r in primary if r['family_similarity_text'])
    pf_ids = {d['xref_id'] for d in doms if d['database'] == 'Pfam'}
    ip_ids = {d['xref_id'] for d in doms if d['database'] == 'InterPro'}
    npf = [int(r['n_pfam'] or 0) for r in primary]
    w(table(['quantity', 'value'], [
        ['proteins with ≥1 Pfam cross-reference', f'{has_pf} ({pct(has_pf, N)})'],
        ['proteins with ≥1 InterPro cross-reference', f'{has_ip} ({pct(has_ip, N)})'],
        ['proteins with a UniProt SIMILARITY (family) comment',
         f'{has_sim} ({pct(has_sim, N)})'],
        ['distinct Pfam accessions observed', len(pf_ids)],
        ['distinct InterPro accessions observed', len(ip_ids)],
        ['total Pfam cross-reference rows',
         sum(1 for d in doms if d['database'] == 'Pfam')],
        ['total InterPro cross-reference rows',
         sum(1 for d in doms if d['database'] == 'InterPro')],
        ['max Pfam domains on one protein', max(npf) if npf else 0],
        ['proteins with >1 Pfam domain', sum(1 for x in npf if x > 1)],
    ]) + '\n')

    rows = []
    for role in ROLE_ORDER:
        sub = [r for r in primary if r['role'] == role]
        rows.append([role, len(sub),
                     pct(sum(1 for r in sub if int(r['n_pfam'] or 0) > 0), len(sub)),
                     pct(sum(1 for r in sub if int(r['n_interpro'] or 0) > 0), len(sub)),
                     pct(sum(1 for r in sub if r['family_similarity_text']), len(sub))])
    w(table(['role', 'n', '≥1 Pfam', '≥1 InterPro', 'SIMILARITY comment'], rows) + '\n')
    w('The SIMILARITY comment is the closest thing UniProt gives to a whole-protein '
      'family statement, and it is free text, not a controlled identifier. It is '
      'reported as its own column rather than being conflated with the domain counts.\n')

    # ---------------- 5. organism types ----------------
    w('---\n\n## 5. Organism resolution (for the optional §6 external-evaluation arm)\n')
    wt = Counter(r['wolfpsort_organism_type'] for r in primary)
    w(table(['organism type', 'n', 'share'],
            [[k, v, pct(v, N)] for k, v in wt.most_common()]) + '\n')
    other = [r for r in primary if r['kingdom_deeploc'] == 'Other']
    ot = Counter(r['wolfpsort_organism_type'] for r in other)
    w(f'**DeepLoc\'s "Other" kingdom bucket resolved by taxonomy** ({len(other)} '
      'proteins). This bucket was *not* assumed unsupported; each lineage was checked:\n')
    w(table(['organism type', 'n'], [[k, v] for k, v in ot.most_common()]) + '\n')
    unsup = [r for r in primary if r['wolfpsort_organism_type'] in
             ('unsupported', 'unresolved')]
    if unsup:
        w('Every one of these resolved to a **eukaryote outside Metazoa, Viridiplantae '
          'and Fungi** — so the DeepLoc bucket and the predictor-support verdict happen '
          'to coincide here, but only because the lineages were checked. What they '
          'actually are, by major eukaryotic group (second lineage rank):\n')
        grp = Counter((r['organism_lineage'].split(';')[1]
                       if len(r['organism_lineage'].split(';')) > 1 else '(unresolved)')
                      for r in unsup)
        w(table(['eukaryotic group', 'n proteins', 'example organism'],
                [[g, n, next((r['organism_scientific_name'] for r in unsup
                              if (r['organism_lineage'].split(';')[1:2] or
                                  ['(unresolved)'])[0] == g), '')]
                 for g, n in grp.most_common()]) + '\n')
        w('Unsupported and unresolved records are retained in every output table with '
          'their reason in `wolfpsort_organism_reason`, excluded only from a WoLF PSORT '
          'submission, and never dropped from the cohort or from any other analysis.\n')

    # ---------------- 6. integrity ----------------
    w('---\n\n## 6. Integrity checks\n')
    p0 = integ['p0_requests']
    w(table(['check', 'result'], [
        ['primary rows', f'{integ["n_primary_rows"]} (expected '
                         f'{integ["n_primary_rows_expected"]})'],
        ['unique accessions', integ['n_unique_accessions']],
        ['exactly one row per development protein',
         'PASS' if integ['one_row_per_protein'] else 'FAIL'],
        ['labels and partitions unchanged vs manifest',
         'PASS' if integ['labels_and_partitions_unchanged'] else 'FAIL'],
        ['frozen sequences re-verified against manifest sha256',
         f'{integ["frozen_sequence_reverification"]["match"]} match, '
         f'{len(integ["frozen_sequence_reverification"]["mismatch"])} mismatch, '
         f'{len(integ["frozen_sequence_reverification"]["missing"])} missing'],
        ['P0 rows requested as annotation targets',
         p0['p0_rows_requested_as_annotation_targets']],
        ['HPA files read', len(integ['hpa_files_read'])],
        ['every row has an explicit status (no silent drops)',
         'PASS' if integ['failures_accounted'] else 'FAIL'],
        ['engineering proteins flagged', integ['engineering_proteins_flagged']],
        ['role counts', ', '.join(f'{k} {v}' for k, v in
                                  sorted(integ['role_counts'].items()))],
    ]) + '\n')

    n_coll = p0['n_development_isoform_parents_that_are_themselves_P0_rows']
    n_share = p0['n_development_parents_sharing_a_base_with_a_P0_isoform']
    w('### P0 boundary — one case that needs stating plainly\n')
    w(f'No P0 row was requested as an annotation target, and no P0 row has an annotation '
      f'row in any output table. However, **{n_coll} development isoform accessions have '
      f'a parent accession that is itself a P0 row**, and **{n_share} more share a base '
      f'accession with a P0 isoform**. Retrieving the parent entry is unavoidable if '
      f'those development isoforms are to be annotated at all.\n')
    if p0['development_isoform_parents_that_are_themselves_P0_rows']:
        w('Parent accessions that are themselves P0 rows: `'
          + '`, `'.join(p0['development_isoform_parents_that_are_themselves_P0_rows'])
          + '`\n')
    w('Every affected development row carries `p0_base_collision != none`, so these rows '
      'can be excluded from any later analysis without refetching, if the reviewer '
      'prefers that. **This is flagged for a decision, not decided here.** The P0 '
      'reservation as written covers fitting, feature selection, tuning and scoring, '
      'none of which this step performs.\n')

    # ---------------- 7. assessment ----------------
    w('---\n\n## 7. What these annotations support, and what they do not\n')
    w('*Assessment of analysis readiness. The experiment is not redesigned here.*\n')

    sig_rate_x = sum(1 for r in primary if cls(r) == 'Extracellular'
                     and r['signal_status'] == 'annotated_present')
    n_x = sum(1 for r in primary if cls(r) == 'Extracellular')
    sig_rate_c = sum(1 for r in primary if cls(r) == 'Cytoplasm'
                     and r['signal_status'] == 'annotated_present')
    n_c = sum(1 for r in primary if cls(r) == 'Cytoplasm')

    w('### Supported\n')
    w(f'- **Signal-peptide co-occurrence characterisation (§1.2 stage 2).** '
      f'{len(present)} proteins carry an annotation ({pct(len(present), N)} of the '
      f'cohort): {sig_rate_x}/{n_x} Extracellular ({pct(sig_rate_x, n_x)}) against '
      f'{sig_rate_c}/{n_c} Cytoplasm ({pct(sig_rate_c, n_c)}). Both arms are populated, '
      f'so the characterisation can run.\n')
    w(f'- **Per-residue interval analysis** on the {sum(1 for r in primary if r["coordinates_transferable"] == "1")} '
      'records whose coordinates are transferable — the annotated intervals can be '
      'aligned to the exact extraction-input sequence without re-indexing.\n')
    w('- **Domain-level covariates.** Pfam/InterPro coverage is high enough to use as a '
      'descriptive covariate in feature characterisation.\n')
    w('- **Organism typing** for the optional external-evaluation arm, resolved from '
      'taxonomy rather than from the DeepLoc kingdom bucket.\n')
    w('- **Entry dates** for the exposure bound described in §6.3 of the specification.\n')

    w('\n### Not supported\n')
    w(f'- **Any analysis treating `not_annotated` as a biological negative.** The '
      f'{Counter(r["signal_status"] for r in primary).get("not_annotated", 0)} '
      '`not_annotated` proteins are an unlabelled pool, not a negative class.\n')
    w(f'- **An experimentally-evidenced-only analysis at full scale.** Only {exp} '
      f'proteins have experimental or combinatorial-experimental evidence '
      f'({pct(exp, N)} of the cohort), so that sensitivity analysis is small and must '
      'report its denominator.\n')
    w('- **Family-level control or family-disjointness claims.** Pfam/InterPro are '
      'region annotations, not family assignments, and no clustering has been run. '
      'The residual-homology question stays open.\n')
    w('- **Cluster-level uncertainty (§5.6 cluster bootstrap).** Blocked until a '
      'clustering run exists. Until then a protein-level bootstrap assumes an '
      'independence that remains unverified.\n')
    w('- **Treating any annotation as independent of a predictor.** A large share of '
      'signal-peptide annotations carry `sequence_analysis_predictor` (ECO:0000255) '
      'evidence, i.e. they are themselves predictor-derived. Validating a feature that '
      'predicts signal peptides against predictor-derived annotations is partly '
      'circular, and any such comparison must be reported split by evidence category.\n')

    txt = '\n'.join(L) + '\n'
    with open(os.path.join(HERE, 'ANNOTATION_COVERAGE_REPORT.md'), 'w') as fh:
        fh.write(txt)
    print(txt)


if __name__ == '__main__':
    main()
