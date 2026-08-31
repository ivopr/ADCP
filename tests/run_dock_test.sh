#!/bin/sh
# =============================================================================
# ADCP docking smoke test -- exercises the AutoDock grid path.
#
# Every other test in this suite runs the folding path (-p Bias=NULL), where the
# receptor term is literally 'extE 0.000000'. Nothing touched the docking code
# until this test existed, which is how a stack buffer overflow that crashed
# every production docking run (main.c swapChains[], fixed) went unnoticed.
#
# Runs the native 3Q47 peptide against the published target, using the exact
# parameter form scripts/runADCP.py:143-171 uses in production.
#
# Note on what is NOT asserted: that extE is negative. Measured across seeds at
# 20k steps it ranges +55.8 to -16.1, and some seeds stay positive even at 200k.
# Individual short runs genuinely fail to find the pocket -- that is why the
# documented protocol runs 50 of them. Asserting a favourable pose here would be
# a flaky test wearing a physics costume. See run_redock_validation.sh for the
# real scientific check.
#
# Usage: run_dock_test.sh <adcp_binary> <staged_target_dir> <workdir>
# =============================================================================
set -eu

ADCP="${1:?Usage: $0 <adcp_binary> <staged_target_dir> <workdir>}"
STAGE="${2:?}"
WD="${3:?}"

if [ ! -f "${STAGE}/rigidReceptor.C.map" ]; then
	echo "SKIP: docking target not staged in ${STAGE}"
	exit 77
fi

rm -rf "${WD}"
mkdir -p "${WD}"
cp "${STAGE}"/rigidReceptor.*.map "${STAGE}/transpoints" "${STAGE}/con" \
   "${STAGE}/ramaprob.data" "${WD}/"
cd "${WD}"

# The native 3Q47 ligand: ASN-PRO-ILE-SER-ASP-VAL-ASP.
SEQ=npisdvd
SEED=4242
OPTS='Bias=NULL,external=5,con,1.0,1.0,Opt=1,0.25,0.75,0.0'

# 250k steps, not a token 20k: the swapChains overflow only manifested past
# swapMutateSteps (200000), so a shorter run passes happily with the bug still
# in place. Verified -- at 20k this test passes against the broken binary, at
# 250k it fails with exit 139. Costs ~4s.
run() {
	# $1 output pdb, $2 log, $3 seed
	"$ADCP" "$SEQ" -r 1x250000 -p "$OPTS" -s "$3" -o "$1" > "$2" 2>&1
}

# --- 1. exit status ---------------------------------------------------------
# On its own this catches the swapChains overflow: crashed runs still wrote a
# plausible 50-model PDB with sensible energies, so only the exit code and the
# missing 'best target energy' line revealed the failure.
set +e
run dock1.pdb dock1.log "$SEED"
rc=$?
set -e
if [ "$rc" -ne 0 ]; then
	# 139 = SIGSEGV, the swapChains overflow signature.
	echo "FAIL: adcp exited $rc on the docking path"
	tail -5 dock1.log
	exit 1
fi

if ! grep -q 'successfully finished' dock1.log; then
	echo "FAIL: run did not reach normal termination"
	tail -5 dock1.log
	exit 1
fi

# --- 2. the receptor was actually loaded ------------------------------------
grep -q 'grid box initialise' dock1.log || {
	echo "FAIL: grid maps were not loaded"; exit 1; }
# AGFR's AutoSite pocket scan is sensitive to the ADFRsuite build/platform that
# generated the target -- ADFRsuite 1.0 on this machine found 104 fill points
# for 3Q47, not the 106 the reference target and this comment used to pin.
# Assert the transpoints file loaded and wasn't the degenerate single-point
# centerXYZ fallback (main.c: "if 0 transpts, use centerXYZ"), not an exact count.
NPTS=$(grep -o 'transpoints initialise success with [0-9]*' dock1.log | grep -o '[0-9]*$')
if [ -z "$NPTS" ] || [ "$NPTS" -lt 50 ]; then
	echo "FAIL: expected at least 50 translation points, got '${NPTS:-none}'"
	grep -i transpoint dock1.log || true
	exit 1
fi

# --- 3. the receptor was actually felt --------------------------------------
# The discriminator against the folding path, which reports exactly 0.000000.
EXT=$(grep 'extE' dock1.log | tail -1 | sed 's/.*extE *//' | awk '{print $1}')
[ -n "$EXT" ] || { echo "FAIL: no external energy reported"; exit 1; }
case "$EXT" in
	*nan*|*NaN*|*inf*|*Inf*) echo "FAIL: non-finite external energy: $EXT"; exit 1 ;;
esac
if awk -v e="$EXT" 'BEGIN { exit !(e == 0) }'; then
	echo "FAIL: extE is exactly 0 -- grid maps are not contributing"
	exit 1
fi

# --- 4. determinism ---------------------------------------------------------
run dock2.pdb dock2.log "$SEED"
grep '^ATOM' dock1.pdb > atoms1.txt
grep '^ATOM' dock2.pdb > atoms2.txt
cmp -s atoms1.txt atoms2.txt || {
	echo "FAIL: seed $SEED produced different docking trajectories"; exit 1; }

# --- 5. a finite total energy was reported ----------------------------------
E=$(grep 'REMARK ENERGY' dock1.pdb | tail -1 | awk '{print $3}')
[ -n "$E" ] || { echo "FAIL: no 'REMARK ENERGY' record in output"; exit 1; }
case "$E" in
	*nan*|*NaN*|*inf*|*Inf*) echo "FAIL: non-finite energy: $E"; exit 1 ;;
esac

echo "PASS: docking ran deterministically, extE $EXT, final energy $E"
