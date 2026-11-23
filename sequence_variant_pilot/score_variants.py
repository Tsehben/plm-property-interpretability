#!/usr/bin/env python3
"""
Score the frozen P3 pilot's originals and variants with the LOCAL, network-disabled
WoLF PSORT container, animal mode. No sequence leaves this machine.

Readout (frozen in generate_variants.py CONFIG["readout"]):
    margin(s) = single-site score("extr") - single-site score("cyto")
    D         = margin(variant) - margin(original)

Extraction -- CORRECTED, and the correction is verified here rather than assumed.
The task asked to verify from the implementation that an omitted class truly means
zero. IT DOES NOT. runWolfPsortSummary prints only a selected subset of the classes
with non-zero score; on packaged animal fixtures 6-12 non-zero classes per record
were omitted from the printed line. Treating an absent class as zero would
therefore be wrong.

So scores are recomputed for EVERY class from the k nearest neighbours that
--print-neighbors reports, through the shipped class-similarity matrix
data/animal.wolfu, using the formula the package documents:

    score(c) = sum over the k neighbours n of  U[c, class(n)]

and it is ASSERTED, for every record, that this reproduces every score the
program itself printed, exactly. The program's own printed numbers are the
check; the recomputation only supplies the classes it chose not to print.

These are predictor scores: similarity-weighted neighbour vote counts. They are
not probabilities and not biological measurements.

Usage:  python3 score_variants.py
"""

import glob
import json
import os
import subprocess
import sys
import time
from collections import Counter

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
OUT = os.path.join(HERE, "v1")
PER = os.path.join(OUT, "per_protein")
WOLF = os.path.join(ROOT, "wolfpsort_setup")
RAW = os.path.join(OUT, "wolfpsort_raw")
sys.path.insert(0, os.path.join(ROOT, "candidate_discovery"))
sys.path.insert(0, WOLF)
import run_discovery as RD  # noqa: E402
import parse_wolfpsort as P  # noqa: E402

ORGANISM = "animal"
WOLFU = os.path.join(WOLF, "src_download", "wolfpsort_repo", "data", f"{ORGANISM}.wolfu")
TOL = 1e-3          # printed scores carry up to 5 decimals (e.g. 12.8333)


def load_similarity():
    with open(WOLFU) as fh:
        classes = fh.readline().split()
        rows = [ln.split() for ln in fh if ln.strip()]
    assert len(rows) == len(classes), (len(rows), len(classes))
    U = {}
    for i, r in enumerate(rows):
        assert len(r) == len(classes)
        for j, v in enumerate(r):
            U[(classes[i], classes[j])] = float(v)
    for c in classes:
        assert U[(c, c)] == 1.0, c
    return classes, U


def parse_with_neighbors(path):
    """-> (k, [{id, printed:{site:score}, neighbours:[(id,class,score)], line}])"""
    k, recs, cur = None, [], None
    with open(path) as fh:
        for raw in fh:
            line = raw.rstrip("\n")
            if not line.strip():
                continue
            m = P.K_LINE.match(line)
            if m:
                k = float(m.group(1))
                continue
            if line.startswith("#"):
                continue
            if line.startswith(" ") or line.startswith("\t"):
                f = line.split()
                if cur is not None and len(f) >= 2:
                    cur["neighbours"].append(
                        (f[0], f[1], float(f[2]) if len(f) > 2 else float("nan")))
                continue
            head, _, rest = line.partition(" ")
            printed, ok = {}, True
            for chunk in rest.split(","):
                chunk = chunk.strip()
                if not chunk:
                    continue
                pm = P.PAIR.match(chunk)
                if not pm or pm.group(1) in printed:
                    ok = False
                    break
                printed[pm.group(1)] = float(pm.group(2))
            cur = dict(id=head, printed=printed, parse_ok=ok and bool(printed),
                       neighbours=[], line=line)
            recs.append(cur)
    return k, recs


def main():
    t0 = time.time()
    os.makedirs(RAW, exist_ok=True)
    classes, U = load_similarity()
    frozen = json.load(open(os.path.join(OUT, "frozen_protocol.json")))

    # ------------------------------------------------------- build the query FASTA
    entries, per = [], {}
    for fp in sorted(glob.glob(os.path.join(PER, "*.json"))):
        d = json.load(open(fp))
        per[d["protein"]] = d
        entries.append((f"{d['protein']}__orig", d["original_sequence"]))
        for v in d["variants"]:
            entries.append((v["variant_id"], v["sequence"]))
    assert len(per) == frozen["sample"]["n"], (len(per), frozen["sample"]["n"])
    expected = frozen["sample"]["n"] * (1 + frozen["config"]["n_replicates"]
                                        * len(frozen["config"]["arms"]))
    assert len(entries) == expected, (len(entries), expected)
    assert len({e[0] for e in entries}) == len(entries), "duplicate query identifiers"
    for _, s in entries:
        assert 30 <= len(s) <= 10000 and "O" not in s, "would be silently dropped"

    fa = os.path.join(OUT, "wolfpsort_queries.fasta")
    with open(fa, "w") as fh:
        for ident, s in entries:
            fh.write(f">{ident}\n")
            for i in range(0, len(s), 60):
                fh.write(s[i:i + 60] + "\n")

    # ------------------------------------------------------------------- run it
    if not os.path.exists(os.path.join(RAW, "summary.txt")):
        env = dict(os.environ, WOLFPSORT_FLAGS="--print-neighbors")
        r = subprocess.run([os.path.join(WOLF, "wolfpsort.sh"), ORGANISM, fa, RAW],
                           capture_output=True, text=True, env=env)
        if r.returncode != 0:
            raise RuntimeError(f"wrapper failed: {r.stderr[:1000]}")
    wolf_seconds = time.time() - t0

    # ---------------------------------------------------- parse and recompute
    k, recs = parse_with_neighbors(os.path.join(RAW, "summary.txt"))
    skipped = P.parse_preflight(os.path.join(RAW, "preflight.txt"))
    by_id = {}
    for r in recs:
        assert r["id"] not in by_id, f"duplicate output id {r['id']}"
        by_id[r["id"]] = r

    integrity = Counter()
    omitted_nonzero = []
    rows = []
    for ident, seq in entries:
        r = by_id.get(ident)
        if r is None:
            rows.append(dict(query_id=ident, status=skipped.get(
                ident, "missing: no output line and no recorded reason")))
            integrity["missing_record"] += 1
            continue
        if not r["parse_ok"]:
            rows.append(dict(query_id=ident, status="unparsable output line"))
            integrity["unparsable"] += 1
            continue
        nb = Counter(c for _, c, _ in r["neighbours"])
        integrity["records_parsed"] += 1
        if sum(nb.values()) != int(k):
            integrity["neighbour_count_mismatch"] += 1
        full = {c: sum(U[(c, cls)] * n for cls, n in nb.items()) for c in classes}
        # the program's own printed numbers verify the recomputation
        for c, v in r["printed"].items():
            assert c in full, (ident, c)
            assert abs(full[c] - v) <= TOL, (ident, c, v, full[c])
            integrity["printed_scores_verified"] += 1
        nz = {c for c, v in full.items() if v > 0}
        om = sorted(nz - set(r["printed"]))
        if om:
            integrity["records_with_nonzero_class_omitted"] += 1
            omitted_nonzero.append(len(om))
        if "extr" in nz and "extr" not in r["printed"]:
            integrity["extr_nonzero_but_omitted"] += 1
        if "cyto" in nz and "cyto" not in r["printed"]:
            integrity["cyto_nonzero_but_omitted"] += 1
        top = max(full.items(), key=lambda t: (t[1], t[0]))[0]
        parts = top.split("_")
        rows.append(dict(
            query_id=ident, status="predicted", k=k,
            n_neighbours=sum(nb.values()),
            extr_score=f"{full['extr']:.17g}", cyto_score=f"{full['cyto']:.17g}",
            margin=f"{full['extr'] - full['cyto']:.17g}",
            n_classes_nonzero=len(nz), n_classes_printed=len(r["printed"]),
            n_nonzero_classes_omitted=len(om),
            top_class_recomputed=top,
            top_class_printed=r["line"].split(" ", 1)[1].split(",")[0].split()[0],
            top_is_compound=int(len(parts) > 1),
            top_mentions_extr=int("extr" in parts), top_mentions_cyto=int("cyto" in parts),
            printed_line=r["line"].split(" ", 1)[1],
            all_scores=json.dumps({c: round(full[c], 6) for c in classes if full[c] > 0},
                                  sort_keys=True)))

    fields = []
    for r in rows:
        for kk in r:
            if kk not in fields:
                fields.append(kk)
    RD.write_tsv(os.path.join(OUT, "wolfpsort_scores.tsv"), fields, rows)

    verdict = dict(
        question="does a class omitted from the printed summary line truly mean zero?",
        answer="NO",
        evidence=dict(
            records_scored=integrity["records_parsed"],
            records_with_at_least_one_nonzero_class_omitted=integrity[
                "records_with_nonzero_class_omitted"],
            max_nonzero_classes_omitted_in_one_record=max(omitted_nonzero, default=0),
            mean_nonzero_classes_omitted=(sum(omitted_nonzero) / len(omitted_nonzero))
            if omitted_nonzero else 0.0,
            records_where_extr_was_nonzero_but_omitted=integrity["extr_nonzero_but_omitted"],
            records_where_cyto_was_nonzero_but_omitted=integrity["cyto_nonzero_but_omitted"],
            printed_scores_verified_against_recomputation=integrity["printed_scores_verified"],
            printed_score_mismatches=0,
        ),
        consequence="margin is computed from the full recomputed class-score vector, not "
                    "from the printed line. Every score the program printed is asserted to "
                    "match the recomputation exactly, so the recomputation is verified "
                    "against the implementation rather than substituted for it.",
        formula="score(c) = sum over the k neighbours n of U[c, class(n)], U = data/animal.wolfu",
    )
    summary = dict(
        organism_model=ORGANISM, k=k,
        n_queries=len(entries), n_predicted=integrity["records_parsed"],
        n_missing=integrity["missing_record"], n_unparsable=integrity["unparsable"],
        neighbour_count_mismatches=integrity["neighbour_count_mismatch"],
        omitted_class_verdict=verdict,
        wolfpsort_seconds=round(wolf_seconds, 1),
        scores_naming="predictor scores (similarity-weighted neighbour vote counts); "
                      "not probabilities, not biological measurements",
        raw_preserved=dict(summary=RD.sha256_file(os.path.join(RAW, "summary.txt")),
                           preflight=RD.sha256_file(os.path.join(RAW, "preflight.txt")),
                           queries=RD.sha256_file(fa)),
    )
    with open(os.path.join(OUT, "wolfpsort_summary.json"), "w") as fh:
        json.dump(summary, fh, indent=2)
    print(json.dumps(summary, indent=2))
    assert integrity["missing_record"] == 0 and integrity["unparsable"] == 0, \
        "records lost; stop rather than proceed"
    return 0


if __name__ == "__main__":
    sys.exit(main())
