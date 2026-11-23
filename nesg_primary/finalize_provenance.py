"""Hash safe outputs; carry forward creation-time restricted hashes without reopening rows."""

from pathlib import Path
import json
import hashlib

ROOT = Path(__file__).resolve().parent
if __name__ == "__main__":
    lock = json.loads((ROOT / "results/test_lock.json").read_text())
    safe = {}
    for p in sorted(ROOT.rglob("*")):
        if (
            not p.is_file()
            or "__pycache__" in str(p)
            or p.is_relative_to(ROOT / "results/reserved_test")
            or p.is_relative_to(ROOT / "results/restricted_prelock")
            or p.name == "output_checksums.json"
        ):
            continue
        safe[str(p.relative_to(ROOT))] = hashlib.sha256(p.read_bytes()).hexdigest()
    (ROOT / "results/output_checksums.json").write_text(
        json.dumps(
            {
                "safe_files": safe,
                "restricted_files_creation_time_only": lock["restricted_checksums"],
            },
            indent=2,
        )
        + "\n"
    )
    print(
        "Hashed",
        len(safe),
        "safe files; retained creation-time hashes for",
        len(lock["restricted_checksums"]),
        "restricted files.",
    )
