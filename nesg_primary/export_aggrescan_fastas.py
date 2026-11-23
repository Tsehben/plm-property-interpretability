"""Export NESG sequences for AGGRESCAN; --test exports the authorized test set.

Only upload eligibility is filtered. Original records, sequences and labels stay
unchanged. Uses frozen sequence artifacts. Test export reads sequences, not phenotype records.
"""

import argparse
import csv
import gzip
import hashlib
import json
from pathlib import Path
import zipfile

ROOT = Path(__file__).resolve().parent
RESULTS = ROOT / "results"
OUT = ROOT / "server_uploads/aggrescan"
AA = set("ACDEFGHIKLMNPQRSTVWY")
BATCH_SIZE = 100  # Operational choice; not a verified current server limit.


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def fasta(records):
    return "".join(f">{r['record_id']}\n{r['sequence']}\n" for r in records)


def write_tsv(path, rows, fields):
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--test",
        action="store_true",
        help="Export reserved-test sequences separately (authorized annotation use).",
    )
    args = parser.parse_args()
    out = OUT / "test" if args.test else OUT
    source = RESULTS / (
        "reserved_test/test_sequences.tsv.gz"
        if args.test
        else "development_records.tsv.gz"
    )
    before = digest(source)
    lock = json.loads((RESULTS / "test_lock.json").read_text())
    pin = (
        lock["restricted_checksums"]["reserved_test/test_sequences.tsv.gz"]
        if args.test
        else lock["rerun_safe_checksums"][source.name]
    )
    assert before == pin
    with gzip.open(source, "rt") as handle:
        # Solubility labels and scores do not enter export selection or headers.
        records = [
            {
                key: row.get(key, "")
                for key in ("record_id", "sid", "sha256", "split", "sequence")
            }
            for row in csv.DictReader(handle, delimiter="\t")
        ]
    if args.test:
        for row in records:
            row["split"] = "reserved_test"
    expected_splits = (
        [("reserved_test", 710)]
        if args.test
        else [("train", 2361), ("validation", 691)]
    )
    assert len(records) == sum(n for _, n in expected_splits)
    assert len({r["record_id"] for r in records}) == len(records)
    assert len({r["sha256"] for r in records}) == len(records)
    assert {r["split"] for r in records} == {s for s, _ in expected_splits}
    for row in records:
        sequence = row["sequence"]
        assert hashlib.sha256(sequence.encode()).hexdigest() == row["sha256"]
        assert not any(c.isspace() for c in row["record_id"])
        row["length"] = len(sequence)
        row["nonstandard_residues"] = "".join(sorted(set(sequence) - AA))
        row["status"] = (
            "nonstandard_residues"
            if row["nonstandard_residues"]
            else "unsupported_length"
            if not 5 <= len(sequence) <= 2000
            else "ready"
        )
        row["fasta_id"] = row["record_id"]
        row["batch_file"] = ""

    out.mkdir(parents=True, exist_ok=True)
    (out / "batches").mkdir(exist_ok=True)
    files = []
    batches = []
    counts = {}
    for split, expected in expected_splits:
        subset = sorted(
            [r for r in records if r["split"] == split], key=lambda r: r["record_id"]
        )
        assert len(subset) == expected
        ready = [r for r in subset if r["status"] == "ready"]
        combined = out / f"nesg_{split}_aggrescan.fasta"
        combined.write_text(fasta(ready))
        files.append(combined)
        for start in range(0, len(ready), BATCH_SIZE):
            batch = ready[start : start + BATCH_SIZE]
            path = out / "batches" / f"{split}_{start // BATCH_SIZE + 1:02d}.fasta"
            for row in batch:
                row["batch_file"] = str(path.relative_to(out))
            path.write_text(fasta(batch))
            files.append(path)
            batches.append(path)
        counts[split] = dict(
            total=len(subset),
            ready=len(ready),
            unsupported=len(subset) - len(ready),
            batches=(len(ready) + BATCH_SIZE - 1) // BATCH_SIZE,
        )

    fields = [
        "record_id",
        "fasta_id",
        "sid",
        "sha256",
        "split",
        "length",
        "status",
        "nonstandard_residues",
        "batch_file",
    ]
    manifest = [{k: r[k] for k in fields} for r in records]
    write_tsv(out / "upload_manifest.tsv", manifest, fields)
    write_tsv(
        out / "unsupported_records.tsv",
        [r for r in manifest if r["status"] != "ready"],
        fields,
    )

    # Verify every exported sequence exactly and each batch membership once.
    lookup = {r["record_id"]: r for r in records}
    batch_ids = []
    for path in files:
        lines = path.read_text().splitlines()
        assert len(lines) % 2 == 0
        ids = []
        for header, sequence in zip(lines[::2], lines[1::2]):
            assert header.startswith(">")
            uid = header[1:]
            assert uid in lookup and lookup[uid]["status"] == "ready"
            assert sequence == lookup[uid]["sequence"]
            ids.append(uid)
        assert len(ids) == len(set(ids))
        if path in batches:
            assert len(ids) <= BATCH_SIZE
            batch_ids.extend(ids)
    assert (
        len(batch_ids)
        == len(set(batch_ids))
        == sum(c["ready"] for c in counts.values())
    )
    assert set(batch_ids) == {r["record_id"] for r in records if r["status"] == "ready"}
    assert digest(source) == before

    instructions = """NESG development sequences for AGGRESCAN

Server: http://bioinf.uab.es/aggrescan/
Reference: https://doi.org/10.1186/1471-2105-8-65

Training: 2,359 ready; 2 unsupported (X).
Validation: 690 ready; 1 unsupported (X).
No reserved-test sequences are included. Unsupported proteins remain in NESG.

Use default server settings. Keep the FASTA identifiers unchanged.
If the page has a text box rather than file upload, paste the FASTA file contents.
The combined files and batches contain the SAME sequences: choose one approach.
Start with batches/train_01.fasta to check the output format. Batch size 100 is
an operational choice, not a verified current server limit. The server could
not be reached from the assistant's browser during preparation.

The ZIP contains batches, manifest and these instructions (no combined FASTAs).
Save each job under its batch name and retain all downloadable result files,
especially per-protein summaries and residue profiles/hotspot coordinates.
Record the run date, settings, server version if shown, and any errors.
Do not replace X, truncate sequences or treat failed jobs as zero aggregation.
Predictions are aggregation annotations, not experimental property labels.

Reproduce from repository root:
  python nesg_primary/export_aggrescan_fastas.py
"""
    if args.test:
        c = counts["reserved_test"]
        instructions = f"""NESG reserved-test sequences for AGGRESCAN

{c["ready"]} of {c["total"]} test sequences ready; {c["unsupported"]} unsupported.
Server: http://bioinf.uab.es/aggrescan/
Use the same default settings as training and validation. Keep identifiers unchanged.
Paste FASTA contents if the server provides a text box instead of file upload.
Choose either the combined FASTA or the batches; both contain the same sequences.
The ZIP contains {c["batches"]} batches of at most {BATCH_SIZE} proteins plus manifests.
Batch size is an operational choice, not a verified server limit.
Save complete outputs, including residue scores and hotspot coordinates, by batch name.
Record run date, settings and any errors. Failed jobs are not zero aggregation.
The manifest uses original record IDs and exact sequence hashes. The sid column is
blank because the sequence-only source does not include those secondary identifiers.
No phenotype records are read or included. Original sequences and splits are unchanged.

Reproduce: python nesg_primary/export_aggrescan_fastas.py --test
"""
    (out / "UPLOAD_INSTRUCTIONS.txt").write_text(instructions)
    archive = out / (
        "nesg_aggrescan_test_batches.zip" if args.test else "nesg_aggrescan_batches.zip"
    )
    archive_files = batches + [
        out / "upload_manifest.tsv",
        out / "unsupported_records.tsv",
        out / "UPLOAD_INSTRUCTIONS.txt",
    ]
    with zipfile.ZipFile(archive, "w", compression=zipfile.ZIP_DEFLATED) as handle:
        for path in archive_files:
            entry = zipfile.ZipInfo(
                str(path.relative_to(out)), date_time=(1980, 1, 1, 0, 0, 0)
            )
            entry.compress_type = zipfile.ZIP_DEFLATED
            handle.writestr(entry, path.read_bytes())
    with zipfile.ZipFile(archive) as handle:
        assert handle.testzip() is None
        for path in archive_files:
            assert handle.read(str(path.relative_to(out))) == path.read_bytes()
    report = dict(
        input_path=str(source.relative_to(ROOT)),
        input_sha256=before,
        script_sha256=digest(Path(__file__)),
        counts=counts,
        batch_size=BATCH_SIZE,
        batch_limit_is_server_verified=False,
        eligibility="20 standard amino acids, 5–2000 residues; no sequence edits",
        selection_uses_solubility_labels=False,
        reserved_test_opened=args.test,
        test_phenotype_records_opened=False,
        test_export_authorization="User requested test FASTAs for fixed aggregation annotation"
        if args.test
        else None,
        submitted_to_server=False,
        original_sequences_unchanged=True,
        fasta_roundtrip_verified=True,
        archive_verified=True,
        outputs_sha256={
            str(p.relative_to(out)): digest(p)
            for p in files + archive_files[len(batches) :] + [archive]
        },
    )
    (out / "export_provenance.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(counts, indent=2))
    print(f"Verified {len(batch_ids):,} sequences in {len(batches)} batches; {archive}")


if __name__ == "__main__":
    main()
