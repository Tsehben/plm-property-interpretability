"""Submit frozen FASTA batches sequentially and verify AGGRESCAN results.

Uses the authors' server, not the local reimplementation. No solubility labels
are loaded. Completed jobs are reused; uncertain POSTs are never auto-retried.
"""

import argparse
import csv
from datetime import datetime, timezone
import fcntl
import gzip
import hashlib
from html.parser import HTMLParser
import io
import json
from pathlib import Path
import re
import subprocess
import time
from urllib.parse import urlencode, urljoin

ROOT = Path(__file__).resolve().parent
UPLOADS = ROOT / "server_uploads/aggrescan"
WORK = ROOT / "aggregation_server"
RESULTS = ROOT / "results/aggrescan"
BASE = "http://andromeda.uab.cat"
ENDPOINT = BASE + "/bioinf/cgi-bin/aap/aap_ov.pl"
METRICS = ["a3vSA", "nHS", "NnHS", "AAT", "THSA", "TA", "AATr", "THSAr", "Na4vSS"]
PROFILE_FIELDS = [
    "record_id",
    "position",
    "amino_acid",
    "a4v",
    "HSA",
    "NHSA",
    "a4vAHS",
    "in_hotspot",
]


def sha(data):
    return hashlib.sha256(data).hexdigest()


def now():
    return datetime.now(timezone.utc).isoformat()


def save_json(path, data):
    temp = path.with_suffix(path.suffix + ".tmp")
    temp.write_text(json.dumps(data, indent=2) + "\n")
    temp.replace(path)


def write_table(path, rows, fields):
    buffer = io.StringIO(newline="")
    writer = csv.DictWriter(
        buffer, fieldnames=fields, delimiter="\t", extrasaction="ignore"
    )
    writer.writeheader()
    writer.writerows(rows)
    data = buffer.getvalue().encode()
    path.write_bytes(gzip.compress(data, mtime=0) if path.suffix == ".gz" else data)


def read_fasta(path):
    records = {}
    for block in path.read_text().split(">")[1:]:
        lines = block.splitlines()
        uid, sequence = lines[0], "".join(lines[1:])
        assert uid not in records
        assert re.fullmatch(r"NESG_raw:\d+", uid), uid
        assert set(sequence) <= set("ACDEFGHIKLMNPQRSTVWY")
        records[uid] = sequence
    assert records
    return records


class ResultTables(HTMLParser):
    """Keep cells from each nested table separate, including red hotspot flags."""

    def __init__(self):
        super().__init__(convert_charrefs=True)
        self.stack = []
        self.tables = []
        self.red = []

    def handle_starttag(self, tag, attrs):
        attrs = dict(attrs)
        if tag == "table":
            self.stack.append({"cells": [], "cell": None})
        elif tag == "td" and self.stack:
            self.stack[-1]["cell"] = {"lines": [""], "red": [False]}
        elif tag == "font":
            self.red.append(attrs.get("color", "").lower() in ("red", "#ff0000"))
        elif tag == "br" and self.stack and self.stack[-1]["cell"] is not None:
            self.stack[-1]["cell"]["lines"].append("")
            self.stack[-1]["cell"]["red"].append(False)

    def handle_endtag(self, tag):
        if tag == "font" and self.red:
            self.red.pop()
        elif tag == "td" and self.stack:
            cell = self.stack[-1]["cell"]
            if cell is not None:
                pairs = [
                    (line.strip(), red)
                    for line, red in zip(cell["lines"], cell["red"])
                    if line.strip()
                ]
                self.stack[-1]["cells"].append(pairs)
                self.stack[-1]["cell"] = None
        elif tag == "table" and self.stack:
            self.tables.append(self.stack.pop()["cells"])

    def handle_data(self, data):
        if self.stack and self.stack[-1]["cell"] is not None:
            cell = self.stack[-1]["cell"]
            cell["lines"][-1] += data
            if data.strip() and any(self.red):
                cell["red"][-1] = True


def parse_result(html, expected):
    parser = ResultTables()
    parser.feed(html)
    # This transformation is observed on this server; sequence verification is mandatory.
    returned_to_original = {uid.replace(":", "3A"): uid for uid in expected}
    assert len(returned_to_original) == len(expected)
    summaries, residues, regions = [], [], []
    seen = set()
    for cells in parser.tables:
        texts = [" ".join(value for value, _ in cell) for cell in cells]
        if not texts or texts[0] not in returned_to_original:
            continue
        returned_id = texts[0]
        uid = returned_to_original[returned_id]
        assert uid not in seen, f"Duplicate returned protein: {uid}"
        seen.add(uid)
        sequence = expected[uid]
        metrics = dict(zip(METRICS, [float(v) for v in texts[2:11]]))
        assert len(metrics) == len(METRICS)
        start = texts.index("#")
        assert texts[start : start + 6] == ["#", "AA", "a4v", "HSA", "NHSA", "a4vAHS"]
        columns = cells[start + 6 : start + 12]
        assert len(columns) == 6 and all(len(c) == len(sequence) for c in columns), uid
        positions = [int(value) for value, _ in columns[0]]
        assert positions == list(range(1, len(sequence) + 1)), uid
        returned_sequence = "".join(value for value, _ in columns[1])
        assert returned_sequence == sequence, f"Sequence mismatch: {uid}"
        mask = [red for _, red in columns[1]]
        spans = []
        left = None
        for pos, active in enumerate(mask + [False], 1):
            if active and left is None:
                left = pos
            elif not active and left is not None:
                spans.append((left, pos - 1))
                left = None
        assert metrics["nHS"].is_integer() and int(metrics["nHS"]) == len(spans), uid
        assert abs(metrics["NnHS"] - 100 * len(spans) / len(sequence)) <= 0.00051, uid
        for left, right in spans:
            assert 1 <= left <= right <= len(sequence) and right - left + 1 >= 5
            assert "P" not in sequence[left - 1 : right]
            regions.append(
                dict(record_id=uid, start=left, end=right, length=right - left + 1)
            )
        profile = [float(v) for v, _ in columns[2]]
        # Both exported values are rounded: a4v to 0.001, Na4vSS to 0.1.
        assert abs(100 * sum(profile) / len(sequence) - metrics["Na4vSS"]) <= 0.101, uid
        for pos in range(len(sequence)):
            residues.append(
                dict(
                    record_id=uid,
                    position=pos + 1,
                    amino_acid=sequence[pos],
                    **{
                        key: float(columns[col][pos][0])
                        for col, key in enumerate(PROFILE_FIELDS[3:7], 2)
                    },
                    in_hotspot=mask[pos],
                )
            )
        summaries.append(
            dict(
                record_id=uid,
                server_id=returned_id,
                sha256=sha(sequence.encode()),
                length=len(sequence),
                **metrics,
                apr_fraction=sum(mask) / len(sequence),
                longest_apr=max((b - a + 1 for a, b in spans), default=0),
                status="ok" if spans else "zero_features",
            )
        )
    assert seen == set(expected), (
        f"Missing result IDs: {sorted(set(expected) - seen)[:10]}"
    )
    return summaries, residues, regions


def verify_spreadsheet(data, summaries):
    rows = list(csv.reader(io.StringIO(data.decode()), delimiter="\t"))
    assert len(rows) == 10 and rows[0][0] == "Sequence Name"
    by_id = {r["server_id"]: r for r in summaries}
    found = set()
    for col, name in enumerate(rows[0][1:], 1):
        if name not in by_id:
            continue
        assert name not in found
        found.add(name)
        for index, key in enumerate(METRICS, 1):
            assert abs(float(rows[index][col]) - by_id[name][key]) < 1e-8
    assert found == set(by_id)


def catalog():
    manifest_rows, jobs = [], []
    for folder in [UPLOADS, UPLOADS / "test"]:
        provenance = json.loads((folder / "export_provenance.json").read_text())
        manifest = folder / "upload_manifest.tsv"
        assert (
            sha(manifest.read_bytes())
            == provenance["outputs_sha256"]["upload_manifest.tsv"]
        )
        with manifest.open() as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        manifest_rows.extend(rows)
        for name in sorted({r["batch_file"] for r in rows if r["batch_file"]}):
            path = folder / name
            assert sha(path.read_bytes()) == provenance["outputs_sha256"][name]
            expected = read_fasta(path)
            selected = {r["record_id"]: r for r in rows if r["batch_file"] == name}
            assert set(expected) == set(selected)
            assert all(
                sha(s.encode()) == selected[uid]["sha256"]
                for uid, s in expected.items()
            )
            split = next(iter(selected.values()))["split"]
            jobs.append(dict(name=path.stem, split=split, path=path, expected=expected))
    assert len(manifest_rows) == len({r["record_id"] for r in manifest_rows}) == 3762
    return manifest_rows, jobs


def curl(url, output, body=None):
    assert url.startswith(BASE + "/bioinf/"), url
    command = [
        "curl",
        "--fail",
        "--silent",
        "--show-error",
        "--max-time",
        "180",
        "--output",
        str(output),
        "--dump-header",
        str(output.with_suffix(".headers.txt")),
    ]
    if body:
        command += [
            "--header",
            "Content-Type: application/x-www-form-urlencoded",
            "--data-binary",
            "@" + str(body),
        ]
    result = subprocess.run(command + [url], capture_output=True, text=True)
    if result.returncode:
        raise RuntimeError(
            f"curl failed ({result.returncode}): {result.stderr.strip()}"
        )


def run_job(job):
    folder = WORK / "jobs" / job["name"]
    folder.mkdir(parents=True, exist_ok=True)
    state_path = folder / "job.json"
    input_hash = sha(job["path"].read_bytes())
    if state_path.exists():
        state = json.loads(state_path.read_text())
        assert state["input_sha256"] == input_hash
        if state["status"] == "complete":
            for name, checksum in state["outputs_sha256"].items():
                assert sha((folder / name).read_bytes()) == checksum
            print(f"REUSE {job['name']} ({len(job['expected'])} proteins)", flush=True)
            return
        raise RuntimeError(
            f"{job['name']} has an uncertain prior request; inspect {state_path} before resubmitting."
        )
    payload = urlencode(
        {"sequence": job["path"].read_text().replace("\n", "\r\n")}
    ).encode()
    body = folder / "request.form"
    body.write_bytes(payload)
    state = dict(
        status="submitting",
        started_UTC=now(),
        endpoint=ENDPOINT,
        input_fasta=str(job["path"].relative_to(ROOT)),
        input_sha256=input_hash,
        request_sha256=sha(payload),
        split=job["split"],
        N=len(job["expected"]),
        settings="Server defaults; browser-style CRLF in sequence field",
        script_sha256=sha(Path(__file__).read_bytes()),
    )
    save_json(state_path, state)
    print(f"SUBMIT {job['name']} ({state['N']} proteins)", flush=True)
    start = time.monotonic()
    try:
        html_path = folder / "response.html"
        curl(ENDPOINT, html_path, body)
        raw = html_path.read_bytes()
        html = raw.decode("utf-8")
        summaries, residues, regions = parse_result(html, job["expected"])
        links = re.findall(r'href="([^"\s]+\.txt\.xls)"', html, flags=re.I)
        assert len(set(links)) == 1, "Missing or ambiguous result spreadsheet"
        spreadsheet_url = urljoin(BASE, links[0])
        spreadsheet = folder / "results.txt.xls"
        curl(spreadsheet_url, spreadsheet)
        verify_spreadsheet(spreadsheet.read_bytes(), summaries)
        compressed = folder / "response.html.gz"
        compressed.write_bytes(gzip.compress(raw, mtime=0))
        assert gzip.decompress(compressed.read_bytes()) == raw
        html_path.unlink()  # Exact received bytes retained in response.html.gz.
        save_json(folder / "protein_results.json", summaries)
        write_table(folder / "residue_scores.tsv.gz", residues, PROFILE_FIELDS)
        write_table(
            folder / "hotspots.tsv", regions, ["record_id", "start", "end", "length"]
        )
        state.update(
            status="complete",
            completed_UTC=now(),
            seconds=round(time.monotonic() - start, 2),
            spreadsheet_url=spreadsheet_url,
            sequence_hashes_verified=True,
            residue_coordinates_verified=True,
            summary_spreadsheet_agrees=True,
            outputs_sha256={
                p.name: sha(p.read_bytes())
                for p in folder.iterdir()
                if p.is_file() and p != state_path
            },
        )
        save_json(state_path, state)
        print(
            f"DONE {job['name']} in {state['seconds']}s; all sequences verified",
            flush=True,
        )
    except Exception as error:
        state.update(status="needs_review", error=str(error), stopped_UTC=now())
        save_json(state_path, state)
        raise


def compile_outputs(manifest, jobs):
    RESULTS.mkdir(exist_ok=True)
    lookup, completed = {}, []
    for job in jobs:
        folder = WORK / "jobs" / job["name"]
        if not (folder / "job.json").exists():
            continue
        state = json.loads((folder / "job.json").read_text())
        if state["status"] != "complete":
            continue
        for name, checksum in state["outputs_sha256"].items():
            assert sha((folder / name).read_bytes()) == checksum
        assert state["input_sha256"] == sha(job["path"].read_bytes())
        for row in json.loads((folder / "protein_results.json").read_text()):
            assert row["record_id"] not in lookup
            lookup[row["record_id"]] = row | {"job": job["name"]}
        completed.append(folder)
    rows = []
    for original in manifest:
        uid = original["record_id"]
        row = dict(
            record_id=uid,
            split=original["split"],
            sha256=original["sha256"],
            length=original["length"],
            status="not_run" if original["status"] == "ready" else "unsupported_input",
            unsupported_residues=original["nonstandard_residues"],
            method="AGGRESCAN authors' web server",
            evidence_type="predicted_aggregation_propensity",
        )
        if uid in lookup:
            assert lookup[uid]["sha256"] == original["sha256"]
            assert lookup[uid]["length"] == int(original["length"])
            row.update(lookup[uid])
        rows.append(row)
    fields = (
        [
            "record_id",
            "split",
            "sha256",
            "length",
            "status",
            "unsupported_residues",
            "method",
            "evidence_type",
            "job",
            "server_id",
        ]
        + METRICS
        + ["apr_fraction", "longest_apr"]
    )
    write_table(RESULTS / "protein_annotations.tsv", rows, fields)
    write_table(RESULTS / "annotation_status.tsv", rows, fields[:8])
    # Concatenate validated job outputs without loading all residue rows at once.
    for filename in ["residue_scores.tsv.gz", "hotspots.tsv"]:
        target = RESULTS / (filename if filename.endswith(".gz") else filename + ".gz")
        with (
            target.open("wb") as raw,
            gzip.GzipFile(fileobj=raw, mode="wb", filename="", mtime=0) as compressed,
        ):
            for index, folder in enumerate(completed):
                opener = gzip.open if filename.endswith(".gz") else open
                with opener(folder / filename, "rb") as source:
                    header = source.readline()
                    if index == 0:
                        compressed.write(header)
                    for line in source:
                        compressed.write(line)
    counts = {}
    for split in ["train", "validation", "reserved_test"]:
        subset = [r for r in rows if r["split"] == split]
        counts[split] = {
            status: sum(r["status"] == status for r in subset)
            for status in ["ok", "zero_features", "unsupported_input", "not_run"]
        }
    save_json(
        RESULTS / "run_summary.json",
        dict(
            updated_UTC=now(),
            counts=counts,
            completed_batches=len(completed),
            expected_batches=len(jobs),
            solubility_labels_used=False,
            original_datasets_modified=False,
            server=BASE,
            script_sha256=sha(Path(__file__).read_bytes()),
            inputs_sha256={
                str(job["path"].relative_to(ROOT)): sha(job["path"].read_bytes())
                for job in jobs
            },
            outputs_sha256={
                p.name: sha(p.read_bytes())
                for p in RESULTS.iterdir()
                if p.suffix in (".tsv", ".gz")
            },
        ),
    )
    print(json.dumps(counts, indent=2), flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--run",
        action="store_true",
        help="Submit pending jobs; otherwise compile saved results only",
    )
    parser.add_argument("--limit", type=int, help="Maximum number of new jobs")
    parser.add_argument("--split", choices=["train", "validation", "reserved_test"])
    parser.add_argument("--validate-pilot", action="store_true")
    args = parser.parse_args()
    WORK.mkdir(exist_ok=True)
    if args.validate_pilot:
        pilot = WORK / "pilot"
        summaries, residues, regions = parse_result(
            (pilot / "crlf_submission.html").read_text(),
            read_fasta(pilot / "input.fasta"),
        )
        verify_spreadsheet((pilot / "results.txt.xls").read_bytes(), summaries)
        print(json.dumps(summaries, indent=2))
        print(
            f"Verified {len(residues)} residues and {len(regions)} hotspot intervals."
        )
        return
    with (WORK / "runner.lock").open("w") as handle:
        fcntl.flock(handle, fcntl.LOCK_EX | fcntl.LOCK_NB)
        manifest, jobs = catalog()
        if args.run:
            submitted = 0
            for job in jobs:
                if args.split and job["split"] != args.split:
                    continue
                state = WORK / "jobs" / job["name"] / "job.json"
                complete = (
                    state.exists()
                    and json.loads(state.read_text())["status"] == "complete"
                )
                if args.limit is not None and submitted >= args.limit and not complete:
                    continue
                run_job(job)
                if not complete:
                    submitted += 1
                    time.sleep(
                        2
                    )  # One active job; avoid timestamp collisions and server bursts.
        compile_outputs(manifest, jobs)


if __name__ == "__main__":
    main()
