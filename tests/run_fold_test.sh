#!/bin/sh
# =============================================================================
# ADCP folding regression test.
#
# Asserts three things the old exit-code-only tests did not:
#   1. same seed -> bit-identical trajectory (catches uninitialised reads,
#      unseeded RNG use, and anything else that makes runs irreproducible)
#   2. a final energy is actually reported
#   3. that energy is finite and physically plausible
#
# Determinism is asserted rather than an exact stored energy because rand()
# differs between glibc and BSD libc; an exact baseline would fail on macOS
# for reasons unrelated to ADCP.
#
# Usage: run_fold_test.sh <adcp_binary> <workdir> <ramaprob.data>
# =============================================================================
set -eu

ADCP="${1:?Usage: $0 <adcp_binary> <workdir> <ramaprob.data>}"
WD="${2:?}"
RAMA="${3:?}"

rm -rf "$WD"; mkdir -p "$WD"
cp "$RAMA" "$WD/ramaprob.data"
cd "$WD"

SEQ=APGVGVAPGVGV
SEED=12345

"$ADCP" -r 5000x1000 -t 2 -s "$SEED" "$SEQ" -p Bias=NULL -o run1.pdb > run1.log 2>&1
"$ADCP" -r 5000x1000 -t 2 -s "$SEED" "$SEQ" -p Bias=NULL -o run2.pdb > run2.log 2>&1

# 1. determinism
grep '^ATOM' run1.pdb > atoms1.txt
grep '^ATOM' run2.pdb > atoms2.txt
if ! cmp -s atoms1.txt atoms2.txt; then
	echo "FAIL: seed $SEED produced different trajectories across two runs"
	exit 1
fi
NATOM=$(wc -l < atoms1.txt | tr -d ' ')
if [ "$NATOM" -lt 100 ]; then
	echo "FAIL: only $NATOM ATOM records written, expected a full trajectory"
	exit 1
fi

# 2. an energy was reported
E=$(grep 'REMARK ENERGY' run1.pdb | tail -1 | awk '{print $3}')
if [ -z "$E" ]; then
	echo "FAIL: no 'REMARK ENERGY' record in output"
	exit 1
fi

# 3. finite and physically plausible for this peptide
case "$E" in
	*nan*|*NaN*|*inf*|*Inf*) echo "FAIL: non-finite final energy: $E"; exit 1 ;;
esac
if ! awk -v e="$E" 'BEGIN { exit !(e > -500 && e < 500) }'; then
	echo "FAIL: final energy out of plausible range: $E"
	exit 1
fi

echo "PASS: deterministic over 2 runs, $NATOM atoms, final energy $E"
