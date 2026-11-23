#!/usr/bin/env bash
# Homology audit: MMseqs2 clustering at three identity thresholds.
# Configuration is frozen in homology_audit/v1/frozen_specification.json and was
# written BEFORE any cluster output was produced or examined.
#
# Usage: run_clustering.sh <repo_root>
set -euo pipefail
export PATH="/opt/homebrew/bin:$PATH"
REPO="$1"
V="$REPO/homology_audit/v1"
FA="$V/development_10139.fasta"
TMPROOT="${TMPDIR:-/tmp}/mmseqs_homology_audit"
mkdir -p "$V/clusters" "$TMPROOT"

mmseqs version

for ID in 30 50 80; do
  MSI=$(python3 -c "print($ID/100)")
  echo "=== min-seq-id $MSI  -c 0.8  --cov-mode 0 ==="
  mmseqs easy-cluster "$FA" "$V/clusters/id${ID}" "$TMPROOT/id${ID}" \
      --min-seq-id "$MSI" \
      -c 0.8 \
      --cov-mode 0 \
      --seq-id-mode 0 \
      --cluster-mode 0 \
      --cluster-reassign 1 \
      --alignment-mode 3 \
      -s 7.5 \
      -e 0.001 \
      --threads 12
  rm -rf "$TMPROOT/id${ID}"
  # easy-cluster leaves the representative/all-seqs FASTAs; keep only the TSV
  rm -f "$V/clusters/id${ID}_all_seqs.fasta" "$V/clusters/id${ID}_rep_seq.fasta"
done

echo
echo "=== cluster TSVs ==="
for ID in 30 50 80; do
  printf 'id%s  lines=%s  clusters=%s\n' "$ID" \
    "$(wc -l < "$V/clusters/id${ID}_cluster.tsv")" \
    "$(cut -f1 "$V/clusters/id${ID}_cluster.tsv" | sort -u | wc -l)"
done
