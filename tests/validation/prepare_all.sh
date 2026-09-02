#!/bin/sh
# Prepare every target in the manifest. Idempotent: a target that already has a
# staged map set is left alone, so this can be re-run after a network blip.
#
# Usage: prepare_all.sh <targets_dir> [set_filter]
set -eu
TARGETS="${1:?Usage: $0 <targets_dir> [set_filter]}"
ONLY="${2:-ABCDEFG}"
SCRIPT_DIR=$(cd "$(dirname "$0")" && pwd)

ok=0; skip=0; fail=0
grep -v '^#' "${SCRIPT_DIR}/manifest.tsv" | awk -F'\t' 'NF>1 && $1!="set"' |
while IFS='	' read SET PDB APO LEN CYC METRIC REF SOURCE; do
	case "$ONLY" in *"$SET"*) ;; *) continue ;; esac
	if [ -f "${TARGETS}/${PDB}/stage/rigidReceptor.C.map" ]; then
		echo "have ${PDB}"; continue
	fi
	printf '%-6s %s ... ' "$SET" "$PDB"
	if "${SCRIPT_DIR}/prepare_target.sh" "$PDB" "$LEN" "${TARGETS}/${PDB}" \
	   > "${TARGETS}/${PDB}.prep.log" 2>&1; then
		echo "$(grep -m1 '^OK:' "${TARGETS}/${PDB}.prep.log" || echo ok)"
	else
		echo "FAILED -- see ${TARGETS}/${PDB}.prep.log"
		grep -m1 -E '^(FAIL|SKIP)' "${TARGETS}/${PDB}.prep.log" 2>/dev/null | sed 's/^/       /'
	fi
done
