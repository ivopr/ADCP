#!/bin/sh
# =============================================================================
# Score one finished target with the ADCP authors' own clusterADCP.
#
# This is the reference scorer, not a reimplementation: clusterADCP.py ships in
# ADFRsuite and is what produced the numbers in the papers. Our binary's output
# is already in the format it parses -- it reads the
#   "Energy = totalE X ( diagnolE Y extE Z firstlastE W) Rotamers: ..."
# lines that src/probe.cpp:1866 writes, field for field.
#
# What it computes, and the constants it uses (read from the source, not guessed):
#   ranking energy   0.25*totalE + 0.75*extE, reported as kcal/mol (x 0.59219)
#   contact cluster  peptide CB (CA for GLY) vs receptor CB within 8 A,
#                    Jaccard >= cutoff, leader clustering seeded by best energy
#   fnc              non-hydrogen atom pairs within 5 A of the reference,
#                    |intersection| / |reference pairs|
#   side chains      rebuilt from the Rotamers: indices before contacts are
#                    counted -- the docked output is backbone + pseudo-gamma
#                    only, so skipping this undercounts every contact
#
# Usage: score_ref.sh <result_dir> [mode] [cutoff]
#   mode    nc (contact clustering, for fnc)  |  rmsd (backbone RMSD clustering)
#   cutoff  0.8 for nc, 2.5 for rmsd
# =============================================================================
set -eu

DIR="${1:?Usage: $0 <result_dir> [nc|rmsd] [cutoff]}"
MODE="${2:-nc}"
CUTOFF="${3:-}"

ADFR="${ADFR_HOME:-$HOME/ADFRsuite-1.0}"
CLUSTER="${ADFR}/CCSBpckgs/ADCP/clusterADCP.py"

[ -x "${ADFR}/bin/pythonsh" ] || { echo "SKIP: no ADFRsuite at ${ADFR}"; exit 77; }
[ -f "$CLUSTER" ] || { echo "SKIP: no clusterADCP.py in ${ADFR}"; exit 77; }
[ -d "$DIR" ] || { echo "FAIL: no such result dir $DIR"; exit 1; }

case "$MODE" in
	nc)   FLAG="-nc";   [ -n "$CUTOFF" ] || CUTOFF=0.8 ;;
	rmsd) FLAG="-rmsd"; [ -n "$CUTOFF" ] || CUTOFF=2.5 ;;
	*) echo "FAIL: mode must be nc or rmsd"; exit 1 ;;
esac

cd "$DIR"
[ -f native_pep.pdb ] || { echo "FAIL: no native_pep.pdb in $DIR"; exit 1; }
[ -f receptor.pdbqt ] || cp ../../targets/"$(basename "$DIR")"/receptor.pdbqt . 2>/dev/null || true
[ -f receptor.pdbqt ] || { echo "FAIL: no receptor.pdbqt in $DIR"; exit 1; }

# clusterADCP takes ONE multi-model pdb, so concatenate the replicas the way
# upstream runADCP.py does -- including its filter: a replica whose best energy
# is more than 20 (internal units) worse than the best replica is dropped before
# clustering. Reproducing that filter matters; without it the cluster ranks
# differ from the published ones.
BEST=$(grep -h '^best target energy' run_*.log 2>/dev/null \
	| awk '{print $4}' | sort -g | head -1)
[ -n "$BEST" ] || { echo "FAIL: no replica reported a best target energy"; exit 1; }

rm -f all_poses.pdb
KEPT=0
for log in run_*.log; do
	e=$(awk '/^best target energy/{v=$4} END{print v}' "$log")
	[ -n "$e" ] || continue
	pdb="${log%.log}.pdb"
	[ -f "$pdb" ] || continue
	if awk -v a="$e" -v b="$BEST" 'BEGIN{exit !(a < b + 20)}'; then
		cat "$pdb" >> all_poses.pdb
		KEPT=$((KEPT + 1))
	fi
done
echo "clustering ${KEPT} replicas (best target energy ${BEST}, keeping < best+20)"

"${ADFR}/bin/pythonsh" "$CLUSTER" \
	-i all_poses.pdb -rec receptor.pdbqt -ref native_pep.pdb \
	"$FLAG" "$CUTOFF" 2>&1 | tee "score_ref_${MODE}.txt"
