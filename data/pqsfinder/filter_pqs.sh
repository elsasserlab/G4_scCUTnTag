#!/usr/bin/env bash
set -euo pipefail

# Filter PQS_scores.mm10.bed by pqsfinder score (column 5, strictly greater
# than the threshold) and write the result to PQS_scores.min<THRESH>.mm10.bed
# in the same directory as this script.

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BED="$SCRIPT_DIR/PQS_scores.mm10.bed"

if [[ ! -f "$BED" ]]; then
  echo "Error: $BED not found" >&2
  exit 1
fi

for thresh in 40 50; do
  out="$SCRIPT_DIR/PQS_scores.min${thresh}.mm10.bed"
  awk -v t="$thresh" '$5 > t' "$BED" > "$out"
  echo "Wrote $out ($(wc -l < "$out") records)"
done
