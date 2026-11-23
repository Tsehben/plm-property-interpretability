#!/usr/bin/env python3
"""
Annotation-readiness collection for the DeepLoc 2.0 localisation development cohort.

Retrieves UniProt metadata for the 10,139 frozen development proteins (partitions
P1-P4) and emits an accession-level table, linked feature/domain tables, coverage
counts, and integrity checks.

Scope limits enforced by this script, not merely documented:
  * Only development_P1_4 rows are used as request sources. P0 rows are never a
    request source. Where a development ISOFORM's parent accession happens to also
    be a P0 row, that is detected, flagged per row, counted, and reported -- not
    silently performed and not silently skipped.
  * HPA is never read.
  * The primary cohort, partitions, labels and frozen sequences are read-only.
  * No accession is dropped. Every one of the 10,139 appears exactly once in the
    primary table with an explicit status.

Network: UniProt REST (public, no authentication required or used). Responses are
cached on disk, so re-running is resumable and performs no redundant requests.

Usage:
    python3 fetch_annotations.py [--limit N] [--no-network]
"""

import argparse
import csv
import gzip
import hashlib
import json
import os
import random
import ssl
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from collections import Counter, defaultdict
from datetime import datetime, timezone

# --------------------------------------------------------------------------
# Paths
# --------------------------------------------------------------------------

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
MANIFEST = os.path.join(ROOT, 'DeepLoc2', 'baseline', 'deeploc_binary_manifest.tsv')
RAW_CSV = os.path.join(ROOT, 'DeepLoc2', 'raw', 'Swissprot_Train_Validation_dataset.csv')
CACHE = os.path.join(HERE, 'cache')
OUT = os.path.join(HERE, 'out')

# The eight proteins used in the intervention validation. Flagged, never dropped.
ENGINEERING = {
    'H2A0K6', 'K7ZGS2', 'O13747', 'O42873',
    'O43037', 'O43447', 'O60462-6', 'O61866',
}

# Partition -> experiment role, per LOCALISATION_EXPERIMENT_SPEC.md section 2.3.
ROLE_BY_PARTITION = {'1': 'discovery', '2': 'discovery', '3': 'tuning', '4': 'evaluation'}

# --------------------------------------------------------------------------
# UniProt REST
# --------------------------------------------------------------------------

BASE = 'https://rest.uniprot.org/uniprotkb'
# Generic client identifier. No personal or project-identifying information.
USER_AGENT = 'annotation-readiness-fetch/1.0 (python-urllib)'

FIELDS = ','.join([
    'accession', 'id', 'reviewed', 'protein_name',
    'organism_name', 'organism_id', 'lineage',
    'length', 'sequence',
    'ft_signal', 'ft_transit',
    'xref_pfam', 'xref_interpro',
    'cc_subcellular_location', 'cc_alternative_products', 'protein_families',
    'date_created', 'date_modified', 'date_sequence_modified', 'version',
])

BATCH = 100          # accessions per search request
MAX_RETRIES = 5
BASE_BACKOFF = 1.5   # seconds; exponential with jitter
POLITE_DELAY = 0.15  # seconds between network requests

# ECO evidence code -> category. Evidence STRENGTH is recorded separately from
# annotation PRESENCE: a UniProt SIGNAL feature is not automatically experimental.
EVIDENCE_CATEGORY = {
    'ECO:0000269': 'experimental',
    'ECO:0007744': 'combinatorial_experimental',
    'ECO:0000250': 'by_similarity',
    'ECO:0000255': 'sequence_analysis_predictor',
    'ECO:0000256': 'automatic_assertion',
    'ECO:0000259': 'automatic_assertion',
    'ECO:0000303': 'author_statement',
    'ECO:0000305': 'curator_inference',
    'ECO:0000312': 'imported',
    'ECO:0000313': 'imported',
}
# Ranked strongest-first, for reporting a single "best" category per protein.
EVIDENCE_RANK = [
    'experimental', 'combinatorial_experimental', 'curator_inference',
    'author_statement', 'by_similarity', 'sequence_analysis_predictor',
    'imported', 'automatic_assertion', 'unspecified',
]

_SSL = ssl.create_default_context()
NETWORK_LOG = []


def sha256_text(s):
    return hashlib.sha256(s.encode('utf-8')).hexdigest()


def cache_path(kind, key):
    d = os.path.join(CACHE, kind)
    os.makedirs(d, exist_ok=True)
    return os.path.join(d, hashlib.sha256(key.encode()).hexdigest()[:40] + '.json.gz')


def http_get(url, kind, cache_key, allow_network=True):
    """GET with on-disk caching, retry and exponential backoff.

    Returns (payload_dict, from_cache). payload holds status, body and the
    retrieval timestamp, so provenance survives in the cache itself.
    """
    cp = cache_path(kind, cache_key)
    if os.path.exists(cp):
        with gzip.open(cp, 'rt', encoding='utf-8') as fh:
            return json.load(fh), True
    if not allow_network:
        raise RuntimeError('cache miss with --no-network: ' + cache_key)

    last_err = None
    for attempt in range(MAX_RETRIES):
        try:
            req = urllib.request.Request(url, headers={
                'User-Agent': USER_AGENT,
                'Accept': 'application/json',
            })
            with urllib.request.urlopen(req, timeout=60, context=_SSL) as resp:
                body = resp.read().decode('utf-8')
                payload = {
                    'url': url,
                    'status': resp.status,
                    'body': body,
                    'retrieved_at': datetime.now(timezone.utc).isoformat(),
                }
            with gzip.open(cp, 'wt', encoding='utf-8') as fh:
                json.dump(payload, fh)
            NETWORK_LOG.append({'url': url, 'status': payload['status'],
                                'attempt': attempt + 1})
            time.sleep(POLITE_DELAY)
            return payload, False
        except urllib.error.HTTPError as e:
            if e.code == 404:
                payload = {'url': url, 'status': 404, 'body': '',
                           'retrieved_at': datetime.now(timezone.utc).isoformat()}
                with gzip.open(cp, 'wt', encoding='utf-8') as fh:
                    json.dump(payload, fh)
                NETWORK_LOG.append({'url': url, 'status': 404, 'attempt': attempt + 1})
                time.sleep(POLITE_DELAY)
                return payload, False
            last_err = f'HTTP {e.code}'
            if e.code in (429, 500, 502, 503, 504):
                wait = e.headers.get('Retry-After')
                delay = float(wait) if (wait and wait.isdigit()) else \
                    BASE_BACKOFF * (2 ** attempt) + random.uniform(0, 0.5)
                time.sleep(delay)
                continue
            raise
        except (urllib.error.URLError, TimeoutError, ssl.SSLError, OSError) as e:
            last_err = repr(e)
            time.sleep(BASE_BACKOFF * (2 ** attempt) + random.uniform(0, 0.5))
    raise RuntimeError(f'exhausted {MAX_RETRIES} retries for {url}: {last_err}')


def fetch_search_batch(accs, allow_network=True):
    query = ' OR '.join('accession:' + a for a in accs)
    params = urllib.parse.urlencode({
        'query': query, 'fields': FIELDS, 'format': 'json', 'size': str(len(accs)),
    })
    url = f'{BASE}/search?{params}'
    payload, cached = http_get(url, 'search', 'search|' + '|'.join(accs), allow_network)
    if payload['status'] != 200:
        return [], payload, cached
    return json.loads(payload['body']).get('results', []), payload, cached


def fetch_single_entry(acc, allow_network=True):
    """Individual lookup. Resolves inactive (deleted/merged) entries, which the
    search endpoint omits entirely."""
    url = f'{BASE}/{urllib.parse.quote(acc)}.json'
    payload, cached = http_get(url, 'entry', 'entry|' + acc, allow_network)
    if payload['status'] != 200 or not payload['body']:
        return None, payload, cached
    try:
        return json.loads(payload['body']), payload, cached
    except json.JSONDecodeError:
        return None, payload, cached


def fetch_isoform_fasta(acc, allow_network=True):
    """Isoform sequence by its own accession. The canonical sequence is NEVER
    substituted for a missing isoform."""
    url = f'{BASE}/{urllib.parse.quote(acc)}.fasta'
    payload, cached = http_get(url, 'isoform', 'isoform|' + acc, allow_network)
    if payload['status'] != 200 or not payload['body'].startswith('>'):
        return None, payload, cached
    lines = payload['body'].splitlines()
    seq = ''.join(l.strip() for l in lines[1:] if l and not l.startswith('>'))
    return (seq or None), payload, cached


# --------------------------------------------------------------------------
# Cohort
# --------------------------------------------------------------------------

def load_cohort():
    with open(MANIFEST, newline='') as fh:
        rows = list(csv.DictReader(fh, delimiter='\t'))
    dev = [r for r in rows if r['split_role'] == 'development_P1_4']
    p0 = [r for r in rows if r['partition'] == '0']
    assert all(r['partition'] != '0' for r in dev), 'P0 row present in development set'
    assert len({r['ACC'] for r in dev}) == len(dev), 'duplicate ACC in development set'
    return rows, dev, {r['ACC'] for r in p0}


def load_frozen_sequences(dev):
    """The exact extraction-input sequences, re-verified against the frozen
    manifest hashes before any comparison with UniProt."""
    want = {r['ACC'] for r in dev}
    seqs = {}
    with open(RAW_CSV, newline='') as fh:
        for r in csv.DictReader(fh):
            if r['ACC'] in want:
                seqs[r['ACC']] = r['Sequence'].strip()
    verified = {'match': 0, 'mismatch': [], 'missing': []}
    for r in dev:
        s = seqs.get(r['ACC'])
        if s is None:
            verified['missing'].append(r['ACC'])
        elif sha256_text(s) == r['sha256']:
            verified['match'] += 1
        else:
            verified['mismatch'].append(r['ACC'])
    return seqs, verified


def split_accession(acc):
    if '-' in acc:
        base, _, suf = acc.partition('-')
        return base, suf
    return acc, ''


def wolfpsort_organism_type(lineage):
    """Map an NCBI lineage onto WoLF PSORT's three organism types.

    Returns (type, reason). 'unsupported' and 'unresolved' are explicit flagged
    states, not grounds for dropping a record. DeepLoc's 'Other' kingdom bucket is
    NOT assumed unsupported -- the lineage decides.
    """
    if not lineage:
        return 'unresolved', 'no lineage retrieved'
    s = set(lineage)
    hits = [n for n in ('Metazoa', 'Viridiplantae', 'Fungi') if n in s]
    if len(hits) == 1:
        return {'Metazoa': 'animal', 'Viridiplantae': 'plant',
                'Fungi': 'fungi'}[hits[0]], f'lineage contains {hits[0]}'
    if len(hits) > 1:
        return 'unresolved', 'lineage matches multiple organism types: ' + ','.join(hits)
    top = lineage[0] if lineage else '?'
    return 'unsupported', f'lineage outside Metazoa/Viridiplantae/Fungi (root: {top})'


# --------------------------------------------------------------------------
# Record parsing
# --------------------------------------------------------------------------

def parse_features(entry, ftype):
    out = []
    for f in entry.get('features', []) or []:
        if f.get('type') != ftype:
            continue
        loc = f.get('location', {}) or {}
        st, en = loc.get('start', {}) or {}, loc.get('end', {}) or {}
        evs = f.get('evidences', []) or []
        out.append({
            'start': st.get('value'), 'start_modifier': st.get('modifier', ''),
            'end': en.get('value'), 'end_modifier': en.get('modifier', ''),
            'description': (f.get('description') or '').replace('\t', ' '),
            'evidences': [{
                'evidence_code': e.get('evidenceCode', ''),
                'evidence_source': e.get('source', ''),
                'evidence_id': e.get('id', ''),
                'evidence_category': EVIDENCE_CATEGORY.get(
                    e.get('evidenceCode', ''), 'unspecified'),
            } for e in evs],
        })
    return out


def parse_xrefs(entry):
    out = []
    for x in entry.get('uniProtKBCrossReferences', []) or []:
        db = x.get('database')
        if db not in ('Pfam', 'InterPro'):
            continue
        name = ''
        for p in x.get('properties', []) or []:
            if p.get('key') in ('EntryName', 'Description'):
                name = p.get('value', '')
        out.append({'database': db, 'xref_id': x.get('id', ''), 'entry_name': name})
    return out


def comment_text(entry, ctype):
    chunks = []
    for c in entry.get('comments', []) or []:
        if c.get('commentType') != ctype:
            continue
        for t in c.get('texts', []) or []:
            if t.get('value'):
                chunks.append(t['value'])
        for loc in c.get('subcellularLocations', []) or []:
            v = (loc.get('location') or {}).get('value')
            if v:
                chunks.append(v)
    return ' | '.join(chunks).replace('\t', ' ')


def best_evidence(cats):
    for c in EVIDENCE_RANK:
        if c in cats:
            return c
    return ''


# --------------------------------------------------------------------------
# Main
# --------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--limit', type=int, default=0,
                    help='process only the first N development proteins (smoke test)')
    ap.add_argument('--no-network', action='store_true',
                    help='fail on cache miss instead of making requests')
    args = ap.parse_args()

    os.makedirs(OUT, exist_ok=True)
    os.makedirs(CACHE, exist_ok=True)
    t0 = time.time()
    started = datetime.now(timezone.utc).isoformat()

    all_rows, dev, p0_accs = load_cohort()
    if args.limit:
        dev = dev[:args.limit]
    frozen_seqs, seq_verify = load_frozen_sequences(dev)
    print(f'[cohort] development rows: {len(dev)}  '
          f'frozen-sequence re-verification: match={seq_verify["match"]} '
          f'mismatch={len(seq_verify["mismatch"])} missing={len(seq_verify["missing"])}',
          flush=True)

    # ---- request plan: sources are development rows only -------------------
    base_of = {r['ACC']: split_accession(r['ACC'])[0] for r in dev}
    bases = sorted(set(base_of.values()))
    p0_base_accs = {split_accession(a)[0] for a in p0_accs}

    # A development isoform's parent accession may itself be a P0 row. Detect,
    # flag and count; never request a P0 row as its own annotation target.
    base_is_p0_row = {b for b in bases if b in p0_accs}
    base_shares_p0 = {b for b in bases if b not in p0_accs and b in p0_base_accs}
    print(f'[plan] unique base accessions to query: {len(bases)}  '
          f'(parents that are themselves P0 rows: {len(base_is_p0_row)}; '
          f'parents sharing a base with a P0 isoform: {len(base_shares_p0)})', flush=True)

    isoform_accs = sorted({r['ACC'] for r in dev if split_accession(r['ACC'])[1]})
    print(f'[plan] isoform accessions needing their own sequence: {len(isoform_accs)}',
          flush=True)

    # ---- stage 1: batched entry retrieval ----------------------------------
    entries, provenance = {}, {}
    n_cached = n_net = 0
    for i in range(0, len(bases), BATCH):
        chunk = bases[i:i + BATCH]
        results, payload, cached = fetch_search_batch(chunk, not args.no_network)
        n_cached += cached
        n_net += (not cached)
        for e in results:
            pa = e.get('primaryAccession')
            if pa:
                entries[pa] = e
                provenance[pa] = {'url': payload['url'],
                                  'retrieved_at': payload['retrieved_at'],
                                  'method': 'search_batch'}
            # A requested accession may be a SECONDARY accession of this entry.
            for sa in e.get('secondaryAccessions', []) or []:
                entries.setdefault(sa, e)
                provenance.setdefault(sa, {'url': payload['url'],
                                           'retrieved_at': payload['retrieved_at'],
                                           'method': 'search_batch_secondary'})
        if (i // BATCH) % 10 == 0:
            print(f'[search] {i + len(chunk)}/{len(bases)} base accessions '
                  f'(cache hits {n_cached}, network {n_net})', flush=True)

    # ---- stage 2: individual lookup for anything unresolved ----------------
    # Two cases need the individual endpoint:
    #   (a) accessions the batch search did not return at all;
    #   (b) accessions the search DID return, but as a stub carrying no sequence.
    # Case (b) is how an inactive (merged or deleted) entry surfaces through the
    # search endpoint. Without this check such an accession would be recorded as
    # 'ok' with an empty sequence instead of as merged/deleted.
    omitted = [b for b in bases if b not in entries]
    stubs = [b for b in bases if b in entries
             and not ((entries[b].get('sequence') or {}).get('value'))]
    for b in stubs:
        entries.pop(b, None)
    missing = sorted(set(omitted) | set(stubs))
    print(f'[resolve] base accessions needing individual lookup: {len(missing)} '
          f'(not returned by search: {len(omitted)}; '
          f'returned without a sequence: {len(stubs)})', flush=True)
    inactive = {}
    for j, acc in enumerate(missing):
        e, payload, cached = fetch_single_entry(acc, not args.no_network)
        n_cached += cached
        n_net += (not cached)
        if e is None:
            inactive[acc] = {'reason_type': 'NOT_FOUND', 'merged_into': '',
                             'http_status': payload['status']}
            provenance[acc] = {'url': payload['url'],
                               'retrieved_at': payload['retrieved_at'],
                               'method': 'entry_lookup'}
        elif e.get('entryType', '').lower().startswith('inactive') or 'inactiveReason' in e:
            ir = e.get('inactiveReason', {}) or {}
            inactive[acc] = {
                'reason_type': ir.get('inactiveReasonType', 'INACTIVE'),
                'merged_into': ','.join(ir.get('mergeDemergeTo', []) or []),
                'deleted_reason': ir.get('deletedReason', ''),
                'http_status': payload['status'],
            }
            provenance[acc] = {'url': payload['url'],
                               'retrieved_at': payload['retrieved_at'],
                               'method': 'entry_lookup_inactive'}
        else:
            entries[acc] = e
            provenance[acc] = {'url': payload['url'],
                               'retrieved_at': payload['retrieved_at'],
                               'method': 'entry_lookup'}
        if j and j % 25 == 0:
            print(f'[resolve] {j}/{len(missing)}', flush=True)

    # ---- stage 3: isoform sequences ----------------------------------------
    iso_seqs, iso_fail = {}, {}
    for j, acc in enumerate(isoform_accs):
        s, payload, cached = fetch_isoform_fasta(acc, not args.no_network)
        n_cached += cached
        n_net += (not cached)
        if s:
            iso_seqs[acc] = s
        else:
            iso_fail[acc] = payload['status']
        provenance[acc] = {'url': payload['url'],
                           'retrieved_at': payload['retrieved_at'],
                           'method': 'isoform_fasta'}
        if j and j % 50 == 0:
            print(f'[isoform] {j}/{len(isoform_accs)}', flush=True)
    print(f'[isoform] resolved {len(iso_seqs)}/{len(isoform_accs)}, '
          f'unresolved {len(iso_fail)}', flush=True)

    # ---- stage 4: assemble -------------------------------------------------
    primary, feat_rows, dom_rows = [], [], []

    for r in dev:
        acc = r['ACC']
        base, suf = split_accession(acc)
        is_iso = bool(suf)
        frozen = frozen_seqs.get(acc)
        frozen_sha = r['sha256']
        e = entries.get(base)

        row = {
            'ACC': acc, 'base_accession': base,
            'is_isoform': int(is_iso), 'isoform_suffix': suf,
            'partition': r['partition'],
            'role': ROLE_BY_PARTITION[r['partition']],
            'target': r['target'], 'kingdom_deeploc': r['kingdom'],
            'frozen_len': r['len'],
            'has_additional_localization': r['has_additional_localization'],
            'engineering_protein': int(acc in ENGINEERING),
            'frozen_sha256': frozen_sha,
            'p0_base_collision': ('base_is_p0_row' if base in base_is_p0_row else
                                  'base_shares_with_p0_isoform' if base in base_shares_p0
                                  else 'none'),
        }

        # -- retrieval status
        if e is None:
            inact = inactive.get(base)
            if inact and inact['reason_type'] == 'MERGED':
                row['retrieval_status'] = 'inactive_merged'
            elif inact and inact['reason_type'] == 'DELETED':
                row['retrieval_status'] = 'inactive_deleted'
            elif inact:
                row['retrieval_status'] = 'not_found'
            else:
                row['retrieval_status'] = 'not_found'
            row['inactive_reason'] = (inact or {}).get('reason_type', '')
            row['merged_into'] = (inact or {}).get('merged_into', '')
            row['deleted_reason'] = (inact or {}).get('deleted_reason', '')
        elif (e.get('primaryAccession') or base) != base:
            # The requested accession resolved to a DIFFERENT primary accession: it is
            # a secondary accession of, or was merged into, the entry we received.
            # Recording this as plain 'ok' would silently substitute another
            # accession's record for the one the cohort names.
            row['retrieval_status'] = 'resolved_via_merged_accession'
            row['inactive_reason'] = 'MERGED_OR_SECONDARY'
            row['merged_into'] = e.get('primaryAccession', '')
            row['deleted_reason'] = ''
        else:
            row['retrieval_status'] = 'ok'
            row['inactive_reason'] = ''
            row['merged_into'] = ''
            row['deleted_reason'] = ''

        if e is None:
            row.update({
                'uniprot_primary_accession': '', 'uniprot_id': '', 'entry_type': '',
                'organism_scientific_name': '', 'organism_taxon_id': '',
                'organism_lineage': '', 'wolfpsort_organism_type': 'unresolved',
                'wolfpsort_organism_reason': 'no record retrieved',
                'entry_version': '', 'first_public_date': '',
                'last_annotation_update': '', 'last_sequence_update': '',
                'uniprot_sequence_length': '', 'uniprot_sequence_sha256': '',
                'sequence_source': 'none', 'sequence_match': 'unverifiable',
                'sequence_match_detail': 'no record retrieved',
                'coordinates_transferable': 0,
                'signal_status': 'unverifiable',
                'signal_status_reason': 'no record retrieved',
                'n_signal_features': 0, 'signal_evidence_best': '',
                'signal_evidence_categories': '', 'parent_entry_has_signal': '',
                'n_transit_features': 0, 'transit_status': 'unverifiable',
                'n_pfam': 0, 'n_interpro': 0, 'pfam_ids': '', 'interpro_ids': '',
                'family_similarity_text': '', 'subcellular_location_text': '',
                'retrieved_at': (provenance.get(base) or {}).get('retrieved_at', ''),
                'query_url': (provenance.get(base) or {}).get('url', ''),
                'query_method': (provenance.get(base) or {}).get('method', ''),
            })
            primary.append(row)
            continue

        org = e.get('organism', {}) or {}
        lineage = org.get('lineage', []) or []
        otype, oreason = wolfpsort_organism_type(lineage)
        audit = e.get('entryAudit', {}) or {}
        canon = (e.get('sequence', {}) or {}).get('value', '') or ''

        # -- sequence resolution: isoforms are never served the canonical sequence
        if is_iso:
            up_seq = iso_seqs.get(acc)
            src = 'isoform_fasta' if up_seq else 'none'
        else:
            up_seq = canon or None
            src = 'canonical' if up_seq else 'none'

        if up_seq is None:
            seq_match = 'unverifiable'
            detail = ('isoform sequence could not be retrieved (HTTP '
                      f'{iso_fail.get(acc, "?")}); canonical NOT substituted'
                      if is_iso else 'entry carries no sequence')
            up_sha, up_len = '', ''
        else:
            up_sha, up_len = sha256_text(up_seq), len(up_seq)
            if up_sha == frozen_sha:
                seq_match, detail = 'exact', 'sha256 identical to frozen input'
            else:
                seq_match = 'mismatch'
                detail = (f'sha256 differs (frozen len {r["len"]}, retrieved len {up_len})')

        # Residue coordinates are only transferable onto a sequence we verified.
        # UniProt feature coordinates are defined on the CANONICAL sequence, so an
        # isoform whose sequence differs from canonical cannot inherit them.
        if seq_match != 'exact':
            transferable = False
            trans_reason = 'sequence does not match frozen input'
        elif is_iso and up_seq != canon:
            transferable = False
            trans_reason = 'isoform sequence differs from canonical; coordinates not transferable'
        else:
            transferable = True
            trans_reason = ''

        sig = parse_features(e, 'Signal')
        tra = parse_features(e, 'Transit peptide')
        xrefs = parse_xrefs(e)

        sig_cats = sorted({ev['evidence_category'] for f in sig for ev in f['evidences']}
                          or ({'unspecified'} if sig else set()))

        # -- signal-peptide status.
        # Gated on SEQUENCE verification, not on accession currency: if a merged
        # accession's successor entry carries exactly our frozen sequence, its
        # annotation does apply to that sequence. The merge is flagged separately in
        # retrieval_status so the substitution is never invisible.
        if seq_match == 'mismatch':
            s_status, s_reason = 'unverifiable', 'retrieved sequence does not match frozen input'
        elif seq_match == 'unverifiable':
            s_status, s_reason = 'unverifiable', detail
        elif not transferable:
            s_status, s_reason = 'unverifiable', trans_reason
        elif sig:
            s_status, s_reason = 'annotated_present', 'signal feature on verified matching sequence'
        else:
            s_status, s_reason = 'not_annotated', (
                'matching record retrieved, no signal feature; '
                'NOT a confirmed biological negative')

        if not transferable or seq_match != 'exact':
            t_status = 'unverifiable'
        elif tra:
            t_status = 'annotated_present'
        else:
            t_status = 'not_annotated'

        pfam = [x for x in xrefs if x['database'] == 'Pfam']
        ipro = [x for x in xrefs if x['database'] == 'InterPro']

        row.update({
            'uniprot_primary_accession': e.get('primaryAccession', ''),
            'uniprot_id': e.get('uniProtkbId', ''),
            'entry_type': e.get('entryType', ''),
            'organism_scientific_name': org.get('scientificName', ''),
            'organism_taxon_id': org.get('taxonId', ''),
            'organism_lineage': ';'.join(lineage),
            'wolfpsort_organism_type': otype,
            'wolfpsort_organism_reason': oreason,
            'entry_version': audit.get('entryVersion', ''),
            'first_public_date': audit.get('firstPublicDate', ''),
            'last_annotation_update': audit.get('lastAnnotationUpdateDate', ''),
            'last_sequence_update': audit.get('lastSequenceUpdateDate', ''),
            'uniprot_sequence_length': up_len,
            'uniprot_sequence_sha256': up_sha,
            'sequence_source': src,
            'sequence_match': seq_match,
            'sequence_match_detail': detail if seq_match != 'exact' else '',
            'coordinates_transferable': int(transferable),
            'signal_status': s_status,
            'signal_status_reason': s_reason,
            'n_signal_features': len(sig),
            'signal_evidence_best': best_evidence(sig_cats),
            'signal_evidence_categories': ';'.join(sig_cats),
            'parent_entry_has_signal': int(bool(sig)),
            'n_transit_features': len(tra),
            'transit_status': t_status,
            'n_pfam': len(pfam), 'n_interpro': len(ipro),
            'pfam_ids': ';'.join(x['xref_id'] for x in pfam),
            'interpro_ids': ';'.join(x['xref_id'] for x in ipro),
            'family_similarity_text': comment_text(e, 'SIMILARITY'),
            'subcellular_location_text': comment_text(e, 'SUBCELLULAR LOCATION'),
            'retrieved_at': (provenance.get(base) or {}).get('retrieved_at', ''),
            'query_url': (provenance.get(base) or {}).get('url', ''),
            'query_method': (provenance.get(base) or {}).get('method', ''),
        })
        primary.append(row)

        for ftype, flist in (('Signal', sig), ('Transit peptide', tra)):
            for f in flist:
                evs = f['evidences'] or [{'evidence_code': '', 'evidence_source': '',
                                          'evidence_id': '',
                                          'evidence_category': 'unspecified'}]
                for ev in evs:
                    feat_rows.append({
                        'ACC': acc, 'base_accession': base, 'feature_type': ftype,
                        'start': f['start'], 'start_modifier': f['start_modifier'],
                        'end': f['end'], 'end_modifier': f['end_modifier'],
                        'description': f['description'],
                        'evidence_code': ev['evidence_code'],
                        'evidence_source': ev['evidence_source'],
                        'evidence_id': ev['evidence_id'],
                        'evidence_category': ev['evidence_category'],
                        'coordinates_apply_to_frozen_sequence': int(transferable),
                    })
        for x in xrefs:
            dom_rows.append({'ACC': acc, 'base_accession': base,
                             'database': x['database'], 'xref_id': x['xref_id'],
                             'entry_name': x['entry_name']})

    # ---- write tables ------------------------------------------------------
    def write_tsv(path, rows, cols):
        with open(path, 'w', newline='') as fh:
            w = csv.DictWriter(fh, fieldnames=cols, delimiter='\t',
                               extrasaction='ignore', lineterminator='\n')
            w.writeheader()
            w.writerows(rows)

    pcols = list(primary[0].keys())
    write_tsv(os.path.join(OUT, 'annotation_primary.tsv'), primary, pcols)
    write_tsv(os.path.join(OUT, 'annotation_features.tsv'), feat_rows, [
        'ACC', 'base_accession', 'feature_type', 'start', 'start_modifier', 'end',
        'end_modifier', 'description', 'evidence_code', 'evidence_source',
        'evidence_id', 'evidence_category', 'coordinates_apply_to_frozen_sequence'])
    write_tsv(os.path.join(OUT, 'annotation_domains.tsv'), dom_rows, [
        'ACC', 'base_accession', 'database', 'xref_id', 'entry_name'])

    failures = [r for r in primary
                if r['retrieval_status'] != 'ok' or r['sequence_match'] != 'exact']
    write_tsv(os.path.join(OUT, 'annotation_failures.tsv'), failures, [
        'ACC', 'base_accession', 'is_isoform', 'partition', 'role', 'target',
        'kingdom_deeploc', 'retrieval_status', 'inactive_reason', 'merged_into',
        'deleted_reason', 'sequence_source', 'sequence_match', 'sequence_match_detail',
        'signal_status', 'signal_status_reason', 'query_url'])

    # ---- coverage ----------------------------------------------------------
    cov = defaultdict(Counter)
    for r in primary:
        key = (r['role'], 'Extracellular' if r['target'] == '1' else 'Cytoplasm',
               r['kingdom_deeploc'])
        c = cov[key]
        c['n'] += 1
        c['signal_' + r['signal_status']] += 1
        c['seqmatch_' + r['sequence_match']] += 1
        c['retrieval_' + r['retrieval_status']] += 1
        if r['signal_status'] == 'annotated_present':
            c['sigev_' + (r['signal_evidence_best'] or 'unspecified')] += 1
        if r['transit_status'] == 'annotated_present':
            c['transit_present'] += 1
        if int(r['n_pfam'] or 0) > 0:
            c['has_pfam'] += 1
        if int(r['n_interpro'] or 0) > 0:
            c['has_interpro'] += 1
        if r['family_similarity_text']:
            c['has_similarity_comment'] += 1
        c['wolf_' + r['wolfpsort_organism_type']] += 1

    allkeys = sorted({k for c in cov.values() for k in c})
    with open(os.path.join(OUT, 'coverage_by_role_class_kingdom.tsv'), 'w',
              newline='') as fh:
        w = csv.writer(fh, delimiter='\t', lineterminator='\n')
        w.writerow(['role', 'localisation_class', 'kingdom_deeploc'] + allkeys)
        for key in sorted(cov):
            w.writerow(list(key) + [cov[key].get(k, 0) for k in allkeys])

    # ---- integrity ---------------------------------------------------------
    manifest_fp = hashlib.sha256(
        '\n'.join(f"{r['ACC']}\t{r['target']}\t{r['partition']}" for r in dev)
        .encode()).hexdigest()
    out_fp = hashlib.sha256(
        '\n'.join(f"{r['ACC']}\t{r['target']}\t{r['partition']}" for r in primary)
        .encode()).hexdigest()

    status_counts = Counter(r['retrieval_status'] for r in primary)
    seq_counts = Counter(r['sequence_match'] for r in primary)
    sig_counts = Counter(r['signal_status'] for r in primary)

    integrity = {
        'n_primary_rows': len(primary),
        'n_primary_rows_expected': len(dev),
        'n_unique_accessions': len({r['ACC'] for r in primary}),
        'one_row_per_protein': len(primary) == len({r['ACC'] for r in primary}) == len(dev),
        'labels_and_partitions_unchanged': manifest_fp == out_fp,
        'label_partition_fingerprint': out_fp,
        'manifest_sha256': hashlib.sha256(open(MANIFEST, 'rb').read()).hexdigest(),
        'frozen_sequence_reverification': {
            'match': seq_verify['match'],
            'mismatch': seq_verify['mismatch'],
            'missing': seq_verify['missing'],
        },
        'p0_requests': {
            'p0_rows_requested_as_annotation_targets': 0,
            'rule': ('request sources are development_P1_4 rows only; the queried set is '
                     'the set of their base accessions'),
            'development_isoform_parents_that_are_themselves_P0_rows':
                sorted(base_is_p0_row),
            'n_development_isoform_parents_that_are_themselves_P0_rows':
                len(base_is_p0_row),
            'development_parents_sharing_a_base_with_a_P0_isoform':
                sorted(base_shares_p0),
            'n_development_parents_sharing_a_base_with_a_P0_isoform':
                len(base_shares_p0),
            'note': ('These parent records were retrieved because a DEVELOPMENT isoform '
                     'requires its parent entry. No P0 row received an annotation row. '
                     'Affected development rows carry p0_base_collision != none.'),
        },
        'retrieval_status_counts': dict(status_counts),
        'sequence_match_counts': dict(seq_counts),
        'signal_status_counts': dict(sig_counts),
        'failures_accounted': (
            sum(status_counts.values()) == len(dev)
            and sum(seq_counts.values()) == len(dev)
            and sum(sig_counts.values()) == len(dev)
        ),
        'n_failure_rows': len(failures),
        'engineering_proteins_flagged': sum(r['engineering_protein'] for r in primary),
        'role_counts': dict(Counter(r['role'] for r in primary)),
        'hpa_files_read': [],
    }
    with open(os.path.join(OUT, 'integrity_checks.json'), 'w') as fh:
        json.dump(integrity, fh, indent=2)

    provenance_doc = {
        'started_at': started,
        'finished_at': datetime.now(timezone.utc).isoformat(),
        'elapsed_s': round(time.time() - t0, 1),
        'script': os.path.basename(__file__),
        'script_sha256': hashlib.sha256(open(__file__, 'rb').read()).hexdigest(),
        'endpoint_base': BASE,
        'fields_requested': FIELDS,
        'user_agent': USER_AGENT,
        'authentication': 'none required or used; UniProt REST is public',
        'batch_size': BATCH, 'max_retries': MAX_RETRIES,
        'backoff': 'exponential from %.1fs with jitter; Retry-After honoured' % BASE_BACKOFF,
        'cache_dir': os.path.relpath(CACHE, ROOT),
        'cache_hits': n_cached, 'network_requests': n_net,
        'n_base_accessions_queried': len(bases),
        'n_isoform_sequences_requested': len(isoform_accs),
        'n_individual_entry_lookups': len(missing),
        'inputs': {
            'manifest': os.path.relpath(MANIFEST, ROOT),
            'manifest_sha256': integrity['manifest_sha256'],
            'sequences': os.path.relpath(RAW_CSV, ROOT),
            'sequences_sha256_note': 'per-protein sha256 re-verified against the manifest',
        },
        'network_log_sample': NETWORK_LOG[:20],
        'n_network_log_entries': len(NETWORK_LOG),
    }
    with open(os.path.join(OUT, 'provenance.json'), 'w') as fh:
        json.dump(provenance_doc, fh, indent=2)

    print('\n=== integrity ===')
    print(json.dumps({k: v for k, v in integrity.items()
                      if k not in ('p0_requests', 'frozen_sequence_reverification')},
                     indent=2))
    print('p0_requests.p0_rows_requested_as_annotation_targets =',
          integrity['p0_requests']['p0_rows_requested_as_annotation_targets'])
    print('development isoform parents that are themselves P0 rows =',
          integrity['p0_requests']
          ['n_development_isoform_parents_that_are_themselves_P0_rows'])
    print(f'\nelapsed {provenance_doc["elapsed_s"]}s  '
          f'cache hits {n_cached}  network {n_net}')


if __name__ == '__main__':
    main()
