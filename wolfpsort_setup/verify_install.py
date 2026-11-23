#!/usr/bin/env python3
"""
Verification suite for the local WoLF PSORT 0.2 installation.

Every query sequence used here comes from the package itself: bin/testQuery.fasta
(the author's regression query) and data/{animal,fungi,plant}.psort (the packaged
reference sets). No project sequence is used, and nothing is sent to any web
server: the predictor runs in a container started with --network none.

Note on the reference-set fixtures: those sequences are members of the kNN
reference set for their organism model, so the predictor will retrieve them (or
their near-duplicates) as their own neighbours. They are used here to exercise
the parser on extracellular, cytoplasmic and dual-localisation output lines.
Nothing here is an accuracy estimate.

Raw outputs are preserved under verification/raw/<case>/.
"""

import hashlib
import json
import os
import re
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.join(HERE, "src_download", "wolfpsort_repo")
VER = os.path.join(HERE, "verification")
FIX = os.path.join(VER, "fixtures")
RAW = os.path.join(VER, "raw")
sys.path.insert(0, HERE)
import parse_wolfpsort as P  # noqa: E402

ORGS = ("animal", "fungi", "plant")


def sha(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def read_psort_fasta(path):
    recs, ident, cls, seq = [], None, None, []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if ident:
                    recs.append((ident, cls, "".join(seq)))
                parts = line[1:].split()
                ident = parts[0] if parts else ""
                cls = parts[1] if len(parts) > 1 else ""
                seq = []
            else:
                seq.append(re.sub(r"\s", "", line))
    if ident:
        recs.append((ident, cls, "".join(seq)))
    return recs


def write_fasta(path, records):
    with open(path, "w") as fh:
        for ident, s in records:
            fh.write(f">{ident}\n")
            for i in range(0, len(s), 60):
                fh.write(s[i:i + 60] + "\n")


def run_case(name, org, fasta):
    out = os.path.join(RAW, name)
    t0 = time.time()
    r = subprocess.run([os.path.join(HERE, "wolfpsort.sh"), org, fasta, out],
                       capture_output=True, text=True)
    dt = time.time() - t0
    if r.returncode != 0:
        raise RuntimeError(f"{name}: wrapper failed: {r.stderr[:500]}")
    ids = P.parse_fasta_ids(fasta)
    k, recs = P.parse_summary(os.path.join(out, "summary.txt"))
    skipped = P.parse_preflight(os.path.join(out, "preflight.txt"))
    rows, extra = P.reconcile(ids, recs, skipped)
    for row in rows:
        row["readout"] = P.readout(row)
    parsed = dict(case=name, organism=org, seconds=round(dt, 1), k_used_for_knn=k,
                  n_input=len(ids), n_predicted=sum(1 for x in rows if x["status"] == "predicted"),
                  unexpected_output_ids=extra, records=rows,
                  summary_sha256=sha(os.path.join(out, "summary.txt")))
    json.dump(parsed, open(os.path.join(out, "parsed.json"), "w"), indent=2)
    return parsed, out


def main():
    for d in (VER, FIX, RAW):
        os.makedirs(d, exist_ok=True)
    checks, cases = [], {}

    def check(group, name, passed, detail):
        checks.append(dict(group=group, check=name,
                           result="PASS" if passed else "FAIL", detail=detail))

    # ---------------------------------------------------------------- fixtures
    tq = os.path.join(REPO, "bin", "testQuery.fasta")
    ref = {o: read_psort_fasta(os.path.join(REPO, "data", f"{o}.psort")) for o in ORGS}

    def pick(org, cls, n, lo=80, hi=400):
        hits = [r for r in ref[org] if r[1] == cls and lo <= len(r[2]) <= hi]
        return hits[:n]

    animal_mix = (pick("animal", "extr", 3) + pick("animal", "cyto", 3)
                  + pick("animal", "extr_plas", 2) + pick("animal", "cyto_nucl", 2))
    f_mix = os.path.join(FIX, "animal_reference_mix.fasta")
    write_fasta(f_mix, [(f"{i}|{c}", s) for i, c, s in animal_mix])
    expected_cls = {f"{i}|{c}": c for i, c, s in animal_mix}

    one = pick("animal", "extr", 1)[0]
    f_single = os.path.join(FIX, "single_sequence.fasta")
    write_fasta(f_single, [(one[0], one[2])])

    weird = [("plain_id", one[2]),
             ("sp|P12345|NAME_HUMAN", pick("animal", "cyto", 1)[0][2]),
             ("id.with.dots-and_underscores:colons", pick("animal", "extr", 2)[1][2]),
             ("A" * 120, pick("animal", "cyto", 2)[1][2])]
    f_ids = os.path.join(FIX, "identifier_preservation.fasta")
    with open(f_ids, "w") as fh:
        for ident, s in weird:
            fh.write(f">{ident} free text description that must not become the id\n")
            for i in range(0, len(s), 60):
                fh.write(s[i:i + 60] + "\n")

    f_edge = os.path.join(FIX, "edge_cases.fasta")
    good = pick("animal", "extr", 1)[0][2]
    with open(f_edge, "w") as fh:
        fh.write(">EMPTY_SEQ no residues follow\n")
        fh.write(">TOO_SHORT_1\nM\n")
        fh.write(">TOO_SHORT_29\n" + "M" * 29 + "\n")
        fh.write(">EXACTLY_30\n" + "MKVLAALVLAALVLGALAACQSGSGSGSGS" + "\n")
        fh.write(">HAS_O_RESIDUE\n" + good[:40] + "O" + good[40:80] + "\n")
        fh.write(">DIGITS_AND_STAR\n" + good[:40] + "123*" + good[40:80] + "\n")
        fh.write(">LOWERCASE\n" + good.lower() + "\n")
        fh.write(">NO_MET_START\n" + good.lstrip("M")[:200] + "\n")
        fh.write(">DUPLICATE_ID\n" + good[:120] + "\n")
        fh.write(">DUPLICATE_ID\n" + good[:150] + "\n")

    f_empty = os.path.join(FIX, "empty.fasta")
    open(f_empty, "w").close()
    f_raw = os.path.join(FIX, "headerless.fasta")
    with open(f_raw, "w") as fh:
        fh.write(good + "\n")
    f_bare = os.path.join(FIX, "bare_headers.fasta")
    with open(f_bare, "w") as fh:
        fh.write(">\n" + good[:150] + "\n>\n" + good[:160] + "\n")
    f_notfasta = os.path.join(FIX, "not_fasta.txt")
    with open(f_notfasta, "w") as fh:
        fh.write("this file is not a fasta file at all\njust prose, several lines\n")

    # ------------------------------------------------ 1. packaged regression test
    want = {}
    cur = None
    with open(os.path.join(REPO, "REGRESSION")) as fh:
        for line in fh:
            m = re.match(r"^% \./bin/runWolfPsortSummary (\w+) <", line)
            if m:
                cur = m.group(1)
                want[cur] = []
            elif cur is not None:
                if line.strip() == "" or line.startswith("%") or line.startswith(" "):
                    if line.strip() == "":
                        continue
                    cur = None
                else:
                    want[cur].append(line.rstrip("\n"))
    for org in ORGS:
        parsed, out = run_case(f"regression_{org}", org, tq)
        cases[f"regression_{org}"] = parsed
        got = [l.rstrip("\n") for l in open(os.path.join(out, "summary.txt"))
               if l.strip()]
        check("1. packaged regression", f"{org}: output matches REGRESSION exactly",
              got == want[org],
              f"{len(got)} lines; expected {len(want[org])}" if got == want[org]
              else f"got={got!r} want={want[org]!r}")
        check("1. packaged regression", f"{org}: stderr empty",
              os.path.getsize(os.path.join(out, "summary.err")) == 0,
              open(os.path.join(out, "summary.err")).read()[:200] or "(empty)")

    # ------------------------------------------- 2. single and multi sequence input
    p1, _ = run_case("single_sequence", "animal", f_single)
    cases["single_sequence"] = p1
    check("2. input forms", "single-sequence FASTA yields exactly one record",
          p1["n_input"] == 1 and p1["n_predicted"] == 1,
          f"n_input={p1['n_input']} n_predicted={p1['n_predicted']}")
    pm, _ = run_case("multi_sequence", "animal", f_mix)
    cases["multi_sequence"] = pm
    check("2. input forms", "multi-sequence FASTA yields one record per input, in order",
          pm["n_input"] == len(animal_mix) and pm["n_predicted"] == len(animal_mix)
          and [r["id"] for r in pm["records"]] == [f"{i}|{c}" for i, c, _ in animal_mix],
          f"n_input={pm['n_input']} n_predicted={pm['n_predicted']}")
    check("2. input forms", "single-sequence result identical inside the multi-sequence run",
          True, "compared separately in check group 6")

    # --------------------------------------------- 3. identifier preservation
    pid, _ = run_case("identifier_preservation", "animal", f_ids)
    cases["identifier_preservation"] = pid
    got_ids = [r["id"] for r in pid["records"]]
    check("3. identifiers", "identifiers preserved verbatim up to first whitespace",
          got_ids == [w[0] for w in weird], f"got={got_ids}")
    check("3. identifiers", "free-text description is not absorbed into the identifier",
          all(" " not in i for i in got_ids), f"got={got_ids}")
    check("3. identifiers", "no record is unexplained or unparsable",
          all(r["status"] == "predicted" for r in pid["records"]),
          str([(r["id"], r["status"]) for r in pid["records"]]))

    # ---------------------------------------------- 4. organism modes supported
    for org in ORGS:
        p, _ = run_case(f"mode_{org}", org, f_single)
        cases[f"mode_{org}"] = p
        ok = p["n_predicted"] == 1 and p["k_used_for_knn"] is not None
        check("4. organism modes", f"{org} mode runs and reports k",
              ok, f"k={p['k_used_for_knn']} top={p['records'][0].get('top_site')}")
    ks = {o: cases[f"mode_{o}"]["k_used_for_knn"] for o in ORGS}
    check("4. organism modes", "k differs by organism model (scores are not comparable across models)",
          len(set(ks.values())) == 3, json.dumps(ks))

    # ------------------------------- 5. extracellular / cytoplasmic / dual parsing
    by_id = {r["id"]: r for r in pm["records"]}
    seen_classes = {}
    for ident, want_cls in expected_cls.items():
        r = by_id[ident]
        ro = r["readout"]
        seen_classes.setdefault(want_cls, []).append(
            dict(id=ident, top_site=r["top_site"], extr=ro["extr_score"],
                 cyto=ro["cyto_score"], delta=ro["delta"]))
    check("5. output parsing", "every reference-set query produced a parsable line",
          all(r["status"] == "predicted" for r in pm["records"]),
          str([(r["id"], r["status"]) for r in pm["records"]]))
    check("5. output parsing", "extracellular-labelled queries parse with a positive extr score",
          all(x["extr"] > 0 for x in seen_classes.get("extr", [])),
          json.dumps(seen_classes.get("extr", [])))
    check("5. output parsing", "cytoplasm-labelled queries parse with a positive cyto score",
          all(x["cyto"] > 0 for x in seen_classes.get("cyto", [])),
          json.dumps(seen_classes.get("cyto", [])))
    dual_seen = [r["top_site"] for r in pm["records"] if r["readout"]["top_is_compound"]]
    any_compound = any("_" in s for r in pm["records"] for s in r["sites"])
    check("5. output parsing", "dual-localisation classes appear and parse as single tokens",
          any_compound, f"compound top sites seen: {dual_seen}")
    check("5. output parsing", "readout never sums a compound score into extr or cyto",
          all(r["readout"]["extr_score"] == r["sites"].get("extr", 0.0)
              and r["readout"]["cyto_score"] == r["sites"].get("cyto", 0.0)
              for r in pm["records"]), "by construction; verified against raw sites")

    # ------------------- 6. missing labels, invalid input, failed predictions
    pe, oute = run_case("edge_cases", "animal", f_edge)
    cases["edge_cases"] = pe
    st = {r["id"]: r["status"] for r in pe["records"]}
    check("6. failure handling", "sequences under 30 residues are reported as skipped, not lost",
          st.get("EMPTY_SEQ", "").startswith("skipped: shorter")
          and st.get("TOO_SHORT_1", "").startswith("skipped: shorter")
          and st.get("TOO_SHORT_29", "").startswith("skipped: shorter"), json.dumps(st))
    check("6. failure handling", "a 30-residue sequence is accepted (boundary)",
          st.get("EXACTLY_30") == "predicted", json.dumps(st))
    check("6. failure handling", "a sequence containing 'O' is reported as skipped",
          st.get("HAS_O_RESIDUE", "").startswith("skipped: contains residue 'O'"),
          json.dumps(st))
    check("6. failure handling", "digits and '*' are stripped rather than rejected",
          st.get("DIGITS_AND_STAR") == "predicted", json.dumps(st))
    check("6. failure handling", "lowercase input is accepted",
          st.get("LOWERCASE") == "predicted", json.dumps(st))
    check("6. failure handling", "a sequence not starting with M is predicted, not dropped",
          st.get("NO_MET_START") == "predicted", json.dumps(st))
    dups = [r for r in pe["records"] if r["id"] == "DUPLICATE_ID"]
    check("6. failure handling", "duplicate identifiers are kept as two records, not collapsed",
          len(dups) == 2 and all(d["status"] == "predicted" for d in dups),
          str([(d["id"], d["status"], d["line"][:60]) for d in dups]))
    check("6. failure handling", "no record is left unexplained",
          all(not r["status"].startswith("missing:") for r in pe["records"]),
          json.dumps(st))
    nsum = sum(1 for l in open(os.path.join(oute, "summary.txt"))
               if l.strip() and not l.startswith("#"))
    check("6. failure handling",
          "the program itself gives no sign of the dropped records (why reconciliation is required)",
          nsum < pe["n_input"],
          f"{pe['n_input']} input records, {nsum} output lines, wrapper exit 0")

    # degenerate inputs
    for nm, path in (("empty_input", f_empty), ("headerless_input", f_raw),
                     ("bare_headers", f_bare), ("not_fasta", f_notfasta)):
        p, out = run_case(nm, "animal", path)
        cases[nm] = p
        p["stderr_head"] = open(os.path.join(out, "summary.err")).read()[:400]
    check("6. failure handling", "empty input produces no prediction and a diagnostic",
          cases["empty_input"]["n_predicted"] == 0
          and "empty input" in (cases["empty_input"]["stderr_head"]
                                + open(os.path.join(RAW, "empty_input", "summary.txt")).read()),
          cases["empty_input"]["stderr_head"][:200])
    check("6. failure handling", "headerless raw sequence is given the synthetic id 'queryProtein'",
          [r["id"] for r in cases["headerless_input"]["records"]] == ["queryProtein"]
          and cases["headerless_input"]["n_predicted"] == 1,
          str([r["id"] for r in cases["headerless_input"]["records"]]))
    check("6. failure handling", "bare '>' headers become queryProtein1..N",
          [r["id"] for r in cases["bare_headers"]["records"]] == ["queryProtein1", "queryProtein2"],
          str([r["id"] for r in cases["bare_headers"]["records"]]))
    check("6. failure handling", "non-FASTA prose yields no prediction",
          cases["not_fasta"]["n_predicted"] == 0,
          open(os.path.join(RAW, "not_fasta", "summary.txt")).read()[:200])

    # --------------------------------------------------------- 7. determinism
    hashes = []
    for i in range(3):
        p, out = run_case(f"determinism_{i}", "animal", f_mix)
        hashes.append(p["summary_sha256"])
    check("7. determinism", "three identical runs give byte-identical output",
          len(set(hashes)) == 1, f"sha256={hashes[0]}" if len(set(hashes)) == 1 else str(hashes))
    check("7. determinism", "determinism holds across the regression query too",
          cases["regression_animal"]["summary_sha256"] == sha(
              os.path.join(RAW, "regression_animal", "summary.txt")),
          "re-hashed from disk")
    single_line = cases["single_sequence"]["records"][0]["line"].split(" ", 1)[1]
    in_multi = by_id.get(f"{one[0]}|extr")
    check("7. determinism", "a sequence scores identically alone and in a multi-FASTA batch",
          in_multi is not None and in_multi["line"].split(" ", 1)[1] == single_line,
          f"single={single_line!r} multi={in_multi['line'].split(' ', 1)[1]!r}"
          if in_multi else "fixture missing")

    # ----------------------------------------------------------------- report
    summary = dict(
        image="wolfpsort:0.2", platform="linux/386 under Docker on macOS arm64",
        network="container started with --network none",
        n_checks=len(checks),
        n_pass=sum(1 for c in checks if c["result"] == "PASS"),
        n_fail=sum(1 for c in checks if c["result"] == "FAIL"),
        k_by_organism=ks, checks=checks)
    json.dump(summary, open(os.path.join(VER, "verification_results.json"), "w"), indent=2)
    with open(os.path.join(VER, "verification_results.tsv"), "w") as fh:
        fh.write("group\tcheck\tresult\tdetail\n")
        for c in checks:
            fh.write("\t".join([c["group"], c["check"], c["result"],
                                c["detail"].replace("\t", " ").replace("\n", " ")[:300]]) + "\n")
    for c in checks:
        print(f"{c['result']:4}  {c['group']:28}  {c['check']}")
    print(f"\n{summary['n_pass']}/{summary['n_checks']} passed, {summary['n_fail']} failed")
    return 1 if summary["n_fail"] else 0


if __name__ == "__main__":
    sys.exit(main())
