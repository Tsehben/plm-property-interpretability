"""Validate complete server annotations against exports; no phenotype access."""

import csv
import gzip
import hashlib
import json
import math
from pathlib import Path

from run_aggrescan_server import METRICS, RESULTS, ROOT, catalog, save_json, sha


def main():
    manifest, jobs = catalog()
    source = {r["record_id"]: r for r in manifest}
    sequences = {
        uid: sequence for job in jobs for uid, sequence in job["expected"].items()
    }
    with (RESULTS / "protein_annotations.tsv").open() as handle:
        proteins = list(csv.DictReader(handle, delimiter="\t"))
    by_id = {r["record_id"]: r for r in proteins}
    assert len(proteins) == len(by_id) == len(source) == 3762
    assert set(by_id) == set(source)
    for uid, row in by_id.items():
        assert all(
            row[key] == source[uid][key] for key in ("sha256", "split", "length")
        )
        if uid in sequences:
            assert row["status"] in ("ok", "zero_features")
            assert all(math.isfinite(float(row[key])) for key in METRICS)
            assert 0 <= float(row["apr_fraction"]) <= 1
        else:
            assert row["status"] == "unsupported_input"
            assert all(row[key] == "" for key in METRICS)

    accumulators = {
        uid: dict(
            count=0, hasher=hashlib.sha256(), score_sum=0.0, hotspot_positions=set()
        )
        for uid in sequences
    }
    residue_count = 0
    with gzip.open(RESULTS / "residue_scores.tsv.gz", "rt") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            uid = row["record_id"]
            accum = accumulators[uid]
            position = int(row["position"])
            assert position == accum["count"] + 1
            assert row["amino_acid"] == sequences[uid][position - 1]
            assert all(
                math.isfinite(float(row[key]))
                for key in ("a4v", "HSA", "NHSA", "a4vAHS")
            )
            assert row["in_hotspot"] in ("True", "False")
            accum["hasher"].update(row["amino_acid"].encode())
            accum["count"] += 1
            accum["score_sum"] += float(row["a4v"])
            if row["in_hotspot"] == "True":
                accum["hotspot_positions"].add(position)
            residue_count += 1
    intervals = {uid: [] for uid in sequences}
    with gzip.open(RESULTS / "hotspots.tsv.gz", "rt") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            uid = row["record_id"]
            start, end = int(row["start"]), int(row["end"])
            assert 1 <= start <= end <= len(sequences[uid])
            assert end - start + 1 == int(row["length"]) >= 5
            assert "P" not in sequences[uid][start - 1 : end]
            intervals[uid].append((start, end))
    for uid, accum in accumulators.items():
        row, sequence = by_id[uid], sequences[uid]
        assert accum["count"] == len(sequence)
        assert accum["hasher"].hexdigest() == row["sha256"]
        assert (
            abs(100 * accum["score_sum"] / len(sequence) - float(row["Na4vSS"]))
            <= 0.101
        )
        spans = intervals[uid]
        assert len(spans) == float(row["nHS"])
        positions = {i for a, b in spans for i in range(a, b + 1)}
        assert sum(b - a + 1 for a, b in spans) == len(positions)
        assert positions == accum["hotspot_positions"]
        assert abs(len(positions) / len(sequence) - float(row["apr_fraction"])) < 1e-12
        assert int(row["longest_apr"]) == max((b - a + 1 for a, b in spans), default=0)

    lock = json.loads((ROOT / "results/test_lock.json").read_text())
    for name, checksum in lock["rerun_safe_checksums"].items():
        assert sha((ROOT / "results" / name).read_bytes()) == checksum, name
    counts = {}
    for split in ("train", "validation", "reserved_test"):
        subset = [r for r in proteins if r["split"] == split]
        counts[split] = dict(
            total=len(subset),
            annotated=sum(r["status"] in ("ok", "zero_features") for r in subset),
            zero_features=sum(r["status"] == "zero_features" for r in subset),
            unsupported=sum(r["status"] == "unsupported_input" for r in subset),
        )
    assert sum(r["annotated"] for r in counts.values()) == 3759
    assert sum(r["unsupported"] for r in counts.values()) == 3
    checks = dict(
        status="PASS",
        counts=counts,
        residue_rows=residue_count,
        hotspot_intervals=sum(map(len, intervals.values())),
        sequence_hashes_verified=True,
        coordinates_and_summary_consistency_verified=True,
        frozen_development_artifacts_unchanged=True,
        test_phenotype_records_opened=False,
        validator_sha256=sha(Path(__file__).read_bytes()),
    )
    save_json(RESULTS / "quality_checks.json", checks)
    print(json.dumps(checks, indent=2))


if __name__ == "__main__":
    main()
