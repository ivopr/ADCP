#!/bin/sh
# Score every finished target with the reference scorer. Idempotent: a target
# that already has score_ref_nc.txt is skipped.
# Usage: score_all.sh <results_dir> [nc|rmsd] [cutoff]
set -eu
RESULTS="${1:?Usage: $0 <results_dir> [nc|rmsd] [cutoff]}"
MODE="${2:-nc}"
CUTOFF="${3:-0.8}"
SCRIPT_DIR=$(cd "$(dirname "$0")" && pwd)
ok=0; skip=0; fail=0
for d in "${RESULTS}"/*/; do
	p=$(basename "$d")
	[ -f "${d}/done" ] || { echo "skip ${p} (no done marker)"; skip=$((skip+1)); continue; }
	[ -f "${d}/score_ref_${MODE}.txt" ] && { echo "have ${p}"; skip=$((skip+1)); continue; }
	printf '%-6s ... ' "$p"
	if "${SCRIPT_DIR}/score_ref.sh" "$d" "$MODE" "$CUTOFF" > "${d}/score_ref.log" 2>&1; then
		echo "$(grep -c '^ ' "${d}/score_ref_${MODE}.txt" 2>/dev/null || echo 0) clusters"
		ok=$((ok+1))
	else
		echo "FAILED -- see ${d}/score_ref.log"
		fail=$((fail+1))
	fi
done
echo
echo "scored ${ok}, skipped ${skip}, failed ${fail}"
