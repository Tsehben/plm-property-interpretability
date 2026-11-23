#!/bin/bash
# Run WoLF PSORT 0.2 locally, in the isolated linux/386 container image.
#
#   ./wolfpsort.sh <animal|plant|fungi> <input.fasta> <output_dir>
#
# Writes, all raw and unmodified:
#   <output_dir>/summary.txt    runWolfPsortSummary stdout
#   <output_dir>/summary.err    runWolfPsortSummary stderr
#   <output_dir>/preflight.txt  checkFastaInput.pl diagnostics, which the normal
#                               pipeline swallows (they go to stdout and are then
#                               consumed by the downstream feature computer)
#   <output_dir>/run.json       provenance
#
# Nothing is sent to any web server. The container has no network.
set -euo pipefail

IMAGE="${WOLFPSORT_IMAGE:-wolfpsort:0.2}"
PLATFORM="linux/386"

if [ "$#" -ne 3 ]; then
    echo "usage: $0 <animal|plant|fungi> <input.fasta> <output_dir>" >&2
    exit 2
fi
ORG="$1"; IN="$2"; OUT="$3"
case "$ORG" in animal|plant|fungi) ;; *) echo "bad organism: $ORG" >&2; exit 2;; esac
[ -r "$IN" ] || { echo "cannot read $IN" >&2; exit 2; }
mkdir -p "$OUT"

# Optional extra runWolfPsortSummary flags, e.g. WOLFPSORT_FLAGS=--print-neighbors
read -r -a EXTRA <<< "${WOLFPSORT_FLAGS:-}"

# --network none: the predictor is purely local; this makes that enforceable.
docker run --rm -i --platform "$PLATFORM" --network none "$IMAGE" "${EXTRA[@]}" "$ORG" \
    < "$IN" > "$OUT/summary.txt" 2> "$OUT/summary.err" || true

docker run --rm -i --platform "$PLATFORM" --network none \
    --entrypoint /opt/wolfpsort/bin/checkFastaInput.pl "$IMAGE" \
    < "$IN" 2>&1 | grep '^#' > "$OUT/preflight.txt" || true

python3 - "$ORG" "$IN" "$OUT" "$IMAGE" <<'PY'
import hashlib, json, os, subprocess, sys
org, inp, out, image = sys.argv[1:5]
def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()
img = subprocess.run(["docker", "image", "inspect", image, "--format", "{{.Id}}"],
                     capture_output=True, text=True).stdout.strip()
json.dump(dict(organism=org, image=image, image_id=img,
               input=os.path.abspath(inp), input_sha256=sha(inp),
               summary_sha256=sha(os.path.join(out, "summary.txt")),
               preflight_sha256=sha(os.path.join(out, "preflight.txt"))),
          open(os.path.join(out, "run.json"), "w"), indent=2)
PY
