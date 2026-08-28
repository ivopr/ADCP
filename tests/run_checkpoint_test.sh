#!/bin/sh
# =============================================================================
# ADCP checkpoint write/restart regression test.
#
# The only coverage of read_in_from_checkpoint, read_checkpoint_entry,
# print_checkpoint_entry and the single-object Chaint/Biasmap reallocs in
# checkpoint_io.cpp. Phase 2 step 3b replaces those reallocs and turns cpoints
# into a std::vector<Chain>; without this test that happens blind.
#
# Note on what is NOT asserted: a restarted run does not reproduce the
# uninterrupted run's later trajectory, because the checkpoint stores the
# sample population but not the RNG state, so the two diverge after the restart
# point. Verified empirically — do not "fix" that by asserting equality here.
# What is deterministic, and what this asserts, is restart-from-the-same-
# checkpoint.
#
# Asserts:
#   1. checkpoint files are written
#   2. a restart reads them and resumes at the checkpointed iteration
#   3. two restarts from the same checkpoint are bit-identical
#
# Usage: run_checkpoint_test.sh <adcp_binary> <workdir> <ramaprob.data> <multimodel.pdb>
# =============================================================================
set -eu

ADCP="${1:?Usage: $0 <adcp_binary> <workdir> <ramaprob.data> <multimodel.pdb>}"
WD="${2:?}"
RAMA="${3:?}"
SNAPSHOTS="${4:?}"

rm -rf "$WD"; mkdir -p "$WD"
cp "$RAMA" "$WD/ramaprob.data"

NMODELS=20
awk -v n="$NMODELS" '/^MODEL/{c++} c<=n{print} c>n{exit}' "$SNAPSHOTS" > "$WD/snapshots.pdb"
cd "$WD"

SEED=12345
COMMON="-n -f snapshots.pdb -r 2x10 -s $SEED -p Bias=NULL -C 2,ckpt"

# 1. write checkpoints
# shellcheck disable=SC2086
"$ADCP" $COMMON -o write.pdb > write.log 2>&1

NCKPT=$(ls ckpt_* 2>/dev/null | wc -l | tr -d ' ')
if [ "$NCKPT" -lt 1 ]; then
	echo "FAIL: no checkpoint files written"
	exit 1
fi
if [ ! -s ckpt_0 ]; then
	echo "FAIL: checkpoint ckpt_0 is empty"
	exit 1
fi

# 2 + 3. restart twice from the same checkpoint; both must resume and agree
RESTART_AT=2
restart() {
	# shellcheck disable=SC2086
	"$ADCP" $COMMON -R "$RESTART_AT" -o "$1.pdb" > "$1.log" 2>&1
}
restart rst1
restart rst2

# Resuming means output starts at the checkpointed iteration, not at 0: a run
# that silently ignored -R would rewrite rst1.pdb_0.
if [ -f rst1.pdb_0 ]; then
	echo "FAIL: restart wrote iteration 0, so it did not resume from the checkpoint"
	exit 1
fi

cat rst1.pdb_* 2>/dev/null | grep '^ATOM' > r1.txt || true
cat rst2.pdb_* 2>/dev/null | grep '^ATOM' > r2.txt || true

NATOM=$(wc -l < r1.txt | tr -d ' ')
if [ "$NATOM" -lt 50 ]; then
	echo "FAIL: restart produced only $NATOM ATOM records"
	exit 1
fi

if ! cmp -s r1.txt r2.txt; then
	echo "FAIL: two restarts from the same checkpoint diverged"
	exit 1
fi

echo "PASS: $NCKPT checkpoints written, restart resumed at iteration $RESTART_AT, $NATOM atoms reproducible"
