#!/bin/sh
# =============================================================================
# Run the validation manifest at one of three protocol tiers.
#
# Tiers, and what each is honestly good for:
#
#   smoke    8 replicas x 200k steps. Proves every target still runs end to end
#            on every code path. Says nothing about docking quality.
#   reduced  half the published replicas, a quarter of the published steps
#            (= 1/8 of the CPU). Enough to catch a regression between commits.
#            NOT comparable with the papers.
#   full     the published protocol, verbatim:
#              linear  80 replicas x 3M steps per amino acid
#              long    80 replicas x 7M steps per amino acid
#              cyclic 300 replicas x 1M steps per amino acid
#            The only tier whose numbers may be put next to a published one.
#
# Parallelism: replicas are independent processes, so the machine is filled at
# whichever level has enough work. A smoke target only has 8 replicas, which
# would leave half a 16-core box idle while targets queued up behind it -- so
# targets run concurrently too, sized to keep CORES busy without oversubscribing.
#
# Resumable: a target with results/<pdb>/done is skipped, so an interrupted
# multi-day run picks up where it stopped.
#
# Usage: run_set.sh <adcp_bin> <targets_dir> <results_dir> <tier> [set_filter]
#   set_filter  optional, e.g. "A" or "BC" -- only rows whose set is in it run
# Env:
#   ADCP_CORES            total cores to use (default: nproc)
#   ADCP_REPLICA_BUDGET   per-replica wall-clock cap in seconds
#   ADCP_MAX_ROTAMERS     rotamers tried per side chain (default 20, 0 = all)
# =============================================================================
set -eu

ADCP="${1:?Usage: $0 <adcp_bin> <targets_dir> <results_dir> <tier> [set_filter]}"
TARGETS="${2:?}"
RESULTS="${3:?}"
TIER="${4:?smoke|reduced|full}"
ONLY="${5:-ABCDEFG}"

SCRIPT_DIR=$(cd "$(dirname "$0")" && pwd)
MANIFEST="${SCRIPT_DIR}/manifest.tsv"
CORES="${ADCP_CORES:-$(nproc 2>/dev/null || echo 4)}"

[ -x "$ADCP" ] || { echo "FAIL: no adcp binary at $ADCP"; exit 1; }
[ -f "$MANIFEST" ] || { echo "FAIL: no manifest at $MANIFEST"; exit 1; }
case "$TIER" in smoke|reduced|full) ;; *) echo "FAIL: unknown tier $TIER"; exit 1 ;; esac

mkdir -p "$RESULTS"
STARTED=$(date -u +%Y-%m-%dT%H:%M:%SZ)

# replicas and steps for a row: <set> <length>
plan() {
	set_id="$1"; len="$2"
	case "$set_id" in
		A)     rep=80;  per=7000000 ;;   # long peptides
		B|C|G) rep=300; per=1000000 ;;   # cyclic
		*)     rep=80;  per=3000000 ;;   # linear
	esac
	steps=$((per * len))
	case "$TIER" in
		smoke)   rep=8; steps=200000 ;;
		reduced) rep=$(( (rep + 1) / 2 )); steps=$((steps / 4)) ;;
	esac
	echo "$rep $steps"
}

# How many targets to run at once, and how many cores each gets. Only the smoke
# tier has fewer replicas per target than cores; reduced and full already have
# 40-300 replicas each, so one target at a time already saturates the machine.
if [ "$TIER" = smoke ]; then
	TPAR=$(( CORES / 8 ))
	[ "$TPAR" -lt 1 ] && TPAR=1
else
	TPAR=1
fi
SLOT_CORES=$(( CORES / TPAR ))
[ "$SLOT_CORES" -lt 1 ] && SLOT_CORES=1

# ---------------------------------------------------------------------------
# Everything for one target. Runs in its own background subshell.
# ---------------------------------------------------------------------------
run_one_target() {
	SET="$1"; PDB="$2"; LEN="$3"; CYC="$4"; METRIC="$5"; REF="$6"; SOURCE="$7"

	tdir="${TARGETS}/${PDB}"
	odir="${RESULTS}/${PDB}"

	SEQ=$(cat "${tdir}/peptide.seq")
	OBS=$(cut -d' ' -f2 "${tdir}/peptide.chain" 2>/dev/null || echo "$LEN")
	set -- $(plan "$SET" "$OBS")
	REPLICAS=$1; STEPS=$2

	# Option string mirrors scripts/runADCP.py:196-202 exactly.
	OPTS='Bias=NULL,external=5,con,1.0,1.0'
	case "$CYC" in
		backbone)    OPTS="${OPTS},external2=4,con,1.0,1.0" ;;
		ss)          OPTS="${OPTS},SSbond=80,2.2,20,0.5" ;;
		backbone+ss) OPTS="${OPTS},external2=4,con,1.0,1.0,SSbond=80,2.2,20,0.5" ;;
	esac
	OPTS="${OPTS},Opt=1,0.25,0.75,0.0"

	# Cap on rotamers tried per side chain. 20 is upstream's `cyclic` value and
	# buys ~4-5x, which is what makes iterating on this benchmark practical at
	# all -- the full sweep drops from ~80 min to ~20. It is a speed-for-accuracy
	# trade, not a free win: the rotamer subset is drawn at RANDOM per call, so
	# the energy becomes an even noisier function of the coordinates. Recorded in
	# run_info.txt so no report can silently mix capped and uncapped runs.
	# ADCP_MAX_ROTAMERS=0 restores the exhaustive scan.
	MAXROT="${ADCP_MAX_ROTAMERS:-20}"
	if [ "$MAXROT" -gt 0 ]; then
		OPTS="${OPTS},MaxRotamers=${MAXROT}"
	fi

	rm -rf "$odir"; mkdir -p "$odir"
	cp "${tdir}/stage/"* "$odir/"
	cp "${tdir}/peptide.pdb" "${odir}/native_pep.pdb"
	cp "${tdir}/receptor.pdbqt" "$odir/" 2>/dev/null || true
	cp "${tdir}/provenance.txt" "$odir/" 2>/dev/null || true

	# Per-replica wall-clock cap. A replica can hang at 100% CPU with no output:
	# reproduced on 1SFI, 3P8F and 4KEL (all backbone+ss) where individual seeds
	# spun for nearly two hours on a 200k-step job. Without a cap one such
	# replica stalls a multi-day run indefinitely.
	#
	# Derived from measurement, not a guess. Measured here: a healthy 200k-step
	# replica takes 278 s alone and 290 s with eight in parallel -- ~700 steps/s
	# wall once map loading (60-120 s, and independent of step count) is paid. An
	# earlier budget of 60 + steps/4000 came from a wrong baseline and killed
	# every replica of five healthy targets; those failures looked like engine
	# bugs and were not.
	#
	# At smoke tier every replica runs the same 200k steps, but the cost of those
	# steps is not the same: measured across all 49 targets, linear peptides up
	# to 20 aa finish a whole 8-replica target in 27-232 s, while 4K1E -- a
	# healthy 14-mer with backbone+disulfide cycles, 0 hung -- needed 950 s. The
	# cyclic terms cost roughly 5-10x per step.
	#
	# A flat 400 s was tried and was wrong: it killed nearly every healthy cyclic
	# replica and inflated the hung count on 1SFI from 2/8 to 7/8, turning "slow"
	# into what looked like "hung". The budget therefore follows the cyclization,
	# which is what actually drives the cost.
	case "$TIER" in
		smoke)
			case "$CYC" in
				none) BUDGET=400 ;;
				*)    BUDGET=1800 ;;
			esac ;;
		*)  BUDGET=$(( 4 * (120 + STEPS / 700) ))
			case "$CYC" in none) ;; *) BUDGET=$((BUDGET * 4)) ;; esac ;;
	esac
	BUDGET="${ADCP_REPLICA_BUDGET:-$BUDGET}"
	[ "$BUDGET" -gt 172800 ] && BUDGET=172800

	t0=$(date +%s)
	(
		cd "$odir"
		i=1; running=0
		while [ "$i" -le "$REPLICAS" ]; do
			# set +e inside the subshell is load-bearing: under the script's
			# errexit a non-zero timeout/crash exit kills the subshell at the
			# timeout command itself, so run_N.rc is never written and a real
			# hang gets miscounted as a generic failure with hung=0.
			( set +e
			  timeout -k 10 "$BUDGET" "$ADCP" "$SEQ" -r "1x${STEPS}" -p "$OPTS" \
				-s "$i" -o "run_${i}.pdb" > "run_${i}.log" 2>&1
			  echo $? > "run_${i}.rc" ) &
			running=$((running + 1))
			if [ "$running" -ge "$SLOT_CORES" ]; then wait; running=0; fi
			i=$((i + 1))
		done
		wait
	)
	t1=$(date +%s)

	bad=0; hung=0; i=1
	while [ "$i" -le "$REPLICAS" ]; do
		rc=$(cat "${odir}/run_${i}.rc" 2>/dev/null || echo missing)
		case "$rc" in
			0)        ;;
			124|137)  hung=$((hung + 1)); bad=$((bad + 1)) ;;
			missing)  hung=$((hung + 1)); bad=$((bad + 1)) ;;
			*)        bad=$((bad + 1)) ;;
		esac
		i=$((i + 1))
	done

	{
		echo "set          ${SET}"
		echo "pdb          ${PDB}"
		echo "sequence     ${SEQ}"
		echo "length       observed=${OBS} published=${LEN}"
		echo "cyclization  ${CYC}"
		echo "tier         ${TIER}"
		echo "replicas     ${REPLICAS}"
		echo "steps        ${STEPS}"
		echo "options      ${OPTS}"
		echo "max_rotamers ${MAXROT}"
		echo "failed_runs  ${bad}"
		echo "hung_runs    ${hung}"
		echo "budget_s     ${BUDGET}"
		echo "seconds      $((t1 - t0))"
		echo "reference    ${REF}  (${SOURCE})"
		echo "metric       ${METRIC}"
	} > "${odir}/run_info.txt"

	if [ "$bad" -eq 0 ]; then
		: > "${odir}/done"
		printf '%-6s %-6s %-20s %3s aa %-12s ok in %ss\n' \
			"$SET" "$PDB" "$SEQ" "$OBS" "$CYC" "$((t1 - t0))"
	elif [ "$hung" -gt 0 ] && [ "$bad" -eq "$hung" ] && [ "$hung" -lt "$REPLICAS" ]; then
		# Some replicas hung but others produced poses. Keep the target and score
		# what survived, but record it -- a silently smaller ensemble would look
		# like a worse engine rather than a hang.
		: > "${odir}/done"
		printf '%-6s %-6s %-20s %3s aa %-12s ok in %ss (%s/%s HUNG, killed)\n' \
			"$SET" "$PDB" "$SEQ" "$OBS" "$CYC" "$((t1 - t0))" "$hung" "$REPLICAS"
	else
		printf '%-6s %-6s %-20s %3s aa %-12s FAILED %s/%s (%s hung)\n' \
			"$SET" "$PDB" "$SEQ" "$OBS" "$CYC" "$bad" "$REPLICAS" "$hung"
	fi
}

echo "tier=${TIER} cores=${CORES} sets=${ONLY} targets_in_parallel=${TPAR} cores_per_target=${SLOT_CORES}"
echo

# ---------------------------------------------------------------------------
# Dispatch. The manifest is read into a list first: a `while read` on the far
# side of a pipe runs in a subshell, which would make job accounting invisible.
# ---------------------------------------------------------------------------
QUEUE=$(mktemp)
grep -v '^#' "$MANIFEST" | awk -F'\t' 'NF>1 && $1!="set"' > "$QUEUE"

inflight=0
while IFS='	' read -r SET PDB APO LEN CYC METRIC REF SOURCE; do
	case "$ONLY" in *"$SET"*) ;; *) continue ;; esac

	if [ -f "${RESULTS}/${PDB}/done" ]; then
		echo "SKIP ${PDB} (already done)"; continue
	fi
	if [ ! -f "${TARGETS}/${PDB}/stage/rigidReceptor.C.map" ] \
	   || [ ! -f "${TARGETS}/${PDB}/peptide.seq" ]; then
		echo "SKIP ${PDB} (target not prepared -- run prepare_target.sh first)"
		continue
	fi

	run_one_target "$SET" "$PDB" "$LEN" "$CYC" "$METRIC" "$REF" "$SOURCE" &
	inflight=$((inflight + 1))
	if [ "$inflight" -ge "$TPAR" ]; then wait; inflight=0; fi
done < "$QUEUE"
wait
rm -f "$QUEUE"

echo
echo "started ${STARTED}, finished $(date -u +%Y-%m-%dT%H:%M:%SZ)"
echo "results in ${RESULTS}"
