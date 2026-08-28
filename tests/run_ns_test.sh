#!/bin/sh
# =============================================================================
# ADCP nested-sampling regression test.
#
# Until this existed, nothing in the suite passed -n, so nested.cpp (1325 lines)
# and checkpoint_io.cpp (991 lines) had zero coverage — including store_chain's
# realloc-grow loop over the live NS population, which Phase 2 step 3 converts
# to std::vector. Feeding a multi-model PDB is what drives that loop: read_in_
# from_pdb calls pdbin once per MODEL and store_chain grows cpoints each time.
#
# Asserts:
#   1. every snapshot in the fixture is read and stored (exercises the growth)
#   2. same seed -> bit-identical NS output across two runs
#   3. NS reports finite, plausible energies
#
# Usage: run_ns_test.sh <adcp_binary> <workdir> <ramaprob.data> <multimodel.pdb>
# =============================================================================
set -eu

ADCP="${1:?Usage: $0 <adcp_binary> <workdir> <ramaprob.data> <multimodel.pdb>}"
WD="${2:?}"
RAMA="${3:?}"
SNAPSHOTS="${4:?}"

rm -rf "$WD"; mkdir -p "$WD"
cp "$RAMA" "$WD/ramaprob.data"

# 20 models is enough to grow cpoints 20 times and still run in seconds.
NMODELS=20
awk -v n="$NMODELS" '/^MODEL/{c++} c<=n{print} c>n{exit}' "$SNAPSHOTS" > "$WD/snapshots.pdb"
cd "$WD"

GOT=$(grep -c '^MODEL' snapshots.pdb || true)
if [ "$GOT" -ne "$NMODELS" ]; then
	echo "FAIL: fixture has $GOT MODEL records, expected $NMODELS"
	exit 1
fi

SEED=12345
run() {
	"$ADCP" -n -f snapshots.pdb -r 2x10 -s "$SEED" -p Bias=NULL -o "$1.pdb" > "$1.log" 2>&1
}
run run1
run run2

# 1. every snapshot was actually read and handed to store_chain.
#    read_in_from_pdb logs one "next chain:" per accepted snapshot.
STORED=$(grep -c 'next chain: NAA' run1.log || true)
if [ "$STORED" -ne "$NMODELS" ]; then
	echo "FAIL: $STORED snapshots stored, expected $NMODELS (cpoints growth loop)"
	exit 1
fi

# NS writes its samples to <outfile>_<iteration>, not to <outfile>.
cat run1.pdb_* 2>/dev/null | grep '^ATOM' > atoms1.txt || true
cat run2.pdb_* 2>/dev/null | grep '^ATOM' > atoms2.txt || true

NATOM=$(wc -l < atoms1.txt | tr -d ' ')
if [ "$NATOM" -lt 100 ]; then
	echo "FAIL: only $NATOM ATOM records written, expected NS samples"
	exit 1
fi

# 2. determinism
if ! cmp -s atoms1.txt atoms2.txt; then
	echo "FAIL: seed $SEED produced different NS output across two runs"
	exit 1
fi

# 3. energies are finite, plausible, and monotonically decreasing.
#    Monotonicity is the defining invariant of nested sampling: each iteration
#    discards the worst-likelihood point, so the reported energy can only fall.
#    It is asserted instead of an exact baseline because these energies are
#    post-MC and therefore rand()-dependent, which differs between glibc and
#    BSD libc -- the same reason run_fold_test.sh stores no baseline energy.
grep -h 'REMARK ENERGY' run1.pdb_* 2>/dev/null | awk '{print $3}' > energies.txt || true
NE=$(wc -l < energies.txt | tr -d ' ')
if [ "$NE" -lt 2 ]; then
	echo "FAIL: only $NE 'REMARK ENERGY' records in NS output, expected the sample sequence"
	exit 1
fi
if ! awk '
	/nan|NaN|inf|Inf/ { print "FAIL: non-finite NS energy: " $1; exit 1 }
	$1 <= -5000 || $1 >= 5000 { print "FAIL: NS energy out of plausible range: " $1; exit 1 }
	NR > 1 && $1 > prev { print "FAIL: NS energy rose from " prev " to " $1; exit 1 }
	{ prev = $1 }
' energies.txt; then
	exit 1
fi

E=$(tail -1 energies.txt)
# The sequence is printed so a conversion can be diffed against the parent
# commit on the same machine. These assertions catch crashes, miscounts and
# nondeterminism, but NOT a small deterministic corruption of the stored
# population -- for that, compare this line across builds.
echo "NS energy sequence: $(tr '\n' ' ' < energies.txt)"
echo "PASS: deterministic NS over 2 runs, $NMODELS snapshots stored, $NATOM atoms, last energy $E"
