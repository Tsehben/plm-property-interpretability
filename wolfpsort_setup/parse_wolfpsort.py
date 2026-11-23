#!/usr/bin/env python3
"""
Parser for runWolfPsortSummary output, plus the extracellular/cytoplasmic
readout defined in WOLFPSORT_READINESS.md.

What the numbers are
--------------------
Per the packaged documentation (bin/runWolfPsortSummary, =head1 DESCRIPTION),
each number "roughly corresponds to the number of nearest neighbours" with that
localisation site, computed through a class-similarity matrix (data/*.wolfu):

    score(c) = sum over the k neighbours n of  U[c, class(n)]

so for a compound class the documentation gives, e.g.

    score(extr_plas) = #extr_plas + 0.5*#extr + 0.5*#plas

They are similarity-weighted neighbour vote counts. They are NOT probabilities,
they are not normalised, and across a line they do not sum to k. Nothing in this
module labels them as probabilities, and callers should not either.

Because score("extr") already absorbs half-credit from every extr_* compound
class, adding a compound score to the single-site score would double count the
same neighbours. The readout below therefore uses the single-site scores only.

Why identifier accounting matters
---------------------------------
checkFastaInput.pl silently removes records: it writes "#Warning:" lines to
stdout, which the downstream feature computer then consumes, so they never reach
the caller, and the exit status is 0 regardless. Records are dropped when the
cleaned sequence is shorter than 30 or longer than 10000 residues, or contains
the residue 'O'. Non-letter characters are deleted before that length check.
A caller that does not reconcile output identifiers against input identifiers
will silently lose proteins.
"""

import json
import re
import sys

SINGLE_SITES = ("extr", "cysk", "E.R.", "golg", "mito", "nucl", "plas",
                "pero", "vacu", "chlo", "lyso", "cyto")
K_LINE = re.compile(r"^#\s*k used for kNN is:\s*([0-9.]+)\s*$")
PAIR = re.compile(r"^(\S+?):?\s+(-?[0-9]+(?:\.[0-9]+)?)$")

SKIP_SHORT = re.compile(r"^#Warning: Sequence (\S+) skipped because length (\d+) less than")
SKIP_LONG = re.compile(r"^#Warning: Sequence (\S+) skipped because length (\d+) greater than")
SKIP_O = re.compile(r"^#Warning: Sequence (\S+) skipped because it contains a non-standard")


def parse_fasta_ids(path):
    """Identifiers exactly as WoLF PSORT derives them: first whitespace-delimited
    token after '>'. Bare '>' lines become queryProtein1, queryProtein2, ...;
    a headerless input becomes 'queryProtein'."""
    ids, unlabeled, seen_header = [], 0, False
    with open(path) as fh:
        for line in fh:
            s = line.lstrip()
            if not s.startswith(">"):
                if not seen_header and s.strip() and not ids:
                    ids.append("queryProtein")
                    seen_header = True
                continue
            seen_header = True
            rest = s[1:].strip()
            if rest:
                ids.append(rest.split()[0])
            else:
                unlabeled += 1
                ids.append(f"queryProtein{unlabeled}")
    return ids


def parse_summary(path):
    """-> (k, [record...]) where record = {id, sites: {site: score}, order: [...], line}."""
    k, recs = None, []
    with open(path) as fh:
        for raw in fh:
            line = raw.rstrip("\n")
            if not line.strip():
                continue
            m = K_LINE.match(line)
            if m:
                k = float(m.group(1))
                continue
            if line.startswith("#"):
                continue
            head, _, rest = line.partition(" ")
            sites, order = {}, []
            ok = True
            for chunk in rest.split(","):
                chunk = chunk.strip()
                if not chunk:
                    continue
                pm = PAIR.match(chunk)
                if not pm:
                    ok = False
                    break
                site, score = pm.group(1), float(pm.group(2))
                # a site repeated on one line would make the readout ambiguous
                if site in sites:
                    ok = False
                    break
                sites[site] = score
                order.append(site)
            if not ok or not sites:
                recs.append(dict(id=head, sites={}, order=[], line=line, parse_ok=False))
            else:
                recs.append(dict(id=head, sites=sites, order=order, line=line, parse_ok=True))
    return k, recs


def parse_preflight(path):
    """-> {id: reason} for records the input checker removed."""
    out = {}
    try:
        fh = open(path)
    except OSError:
        return out
    with fh:
        for line in fh:
            for rx, why in ((SKIP_SHORT, "skipped: shorter than 30 residues"),
                            (SKIP_LONG, "skipped: longer than 10000 residues"),
                            (SKIP_O, "skipped: contains residue 'O'")):
                m = rx.match(line.rstrip("\n"))
                if m:
                    out[m.group(1)] = why
    return out


def reconcile(input_ids, recs, skipped):
    """One row per INPUT record, in input order. Output identifiers are matched
    positionally within each identifier, so duplicate identifiers are preserved
    rather than collapsed."""
    from collections import defaultdict
    pool = defaultdict(list)
    for r in recs:
        pool[r["id"]].append(r)
    cursor = defaultdict(int)
    rows = []
    for i, iid in enumerate(input_ids):
        j = cursor[iid]
        cursor[iid] += 1
        r = pool[iid][j] if j < len(pool[iid]) else None
        if r is None:
            status = skipped.get(iid, "missing: no output line and no recorded reason")
            rows.append(dict(input_index=i, id=iid, status=status, sites={},
                             top_site=None, top_score=None, line=""))
        elif not r["parse_ok"]:
            rows.append(dict(input_index=i, id=iid, status="unparsable output line",
                             sites={}, top_site=None, top_score=None, line=r["line"]))
        else:
            top = r["order"][0]           # the program prints sites in rank order
            rows.append(dict(input_index=i, id=iid, status="predicted", sites=r["sites"],
                             top_site=top, top_score=r["sites"][top], line=r["line"]))
    extra = [r["id"] for r in recs if cursor[r["id"]] == 0]
    return rows, sorted(set(extra))


def readout(row):
    """Extracellular-minus-cytoplasmic readout for one reconciled row.

    E = score of the single-site class 'extr' (0 if the class is absent, which is
        a genuine zero: no neighbour carried that site)
    C = score of the single-site class 'cyto', likewise
    delta = E - C, in similarity-weighted neighbour votes

    Returns None when the record has no prediction, so the caller can drop the
    pair and count it rather than impute a zero.

    Compound classes are NOT added in: score('extr') already contains half-credit
    for every extr_* neighbour, so adding extr_plas would double count. They are
    reported as an ambiguity flag instead.
    """
    if row["status"] != "predicted":
        return None
    s = row["sites"]
    e, c = float(s.get("extr", 0.0)), float(s.get("cyto", 0.0))
    top = row["top_site"]
    parts = top.split("_")
    return dict(
        extr_score=e, cyto_score=c, delta=e - c,
        top_site=top,
        top_is_compound=len(parts) > 1,
        top_mentions_extr="extr" in parts,
        top_mentions_cyto="cyto" in parts,
        top_ambiguous_for_readout=("extr" in parts and "cyto" in parts),
        extr_class_absent="extr" not in s,
        cyto_class_absent="cyto" not in s)


def main():
    if len(sys.argv) not in (3, 4):
        print("usage: parse_wolfpsort.py <input.fasta> <summary.txt> [preflight.txt]",
              file=sys.stderr)
        return 2
    ids = parse_fasta_ids(sys.argv[1])
    k, recs = parse_summary(sys.argv[2])
    skipped = parse_preflight(sys.argv[3]) if len(sys.argv) == 4 else {}
    rows, extra = reconcile(ids, recs, skipped)
    for r in rows:
        r["readout"] = readout(r)
    print(json.dumps(dict(k_used_for_knn=k, n_input_records=len(ids),
                          n_predicted=sum(1 for r in rows if r["status"] == "predicted"),
                          unexpected_output_ids=extra, records=rows), indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())
