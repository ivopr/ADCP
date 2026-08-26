#!/bin/sh
# =============================================================================
# ADCP redocking validation -- the scientific test.
#
# Replicates the documented protocol (ccsb.scripps.edu/adcpv11/documentation):
# many independent Monte Carlo searches, ranked by the target energy, with the
# top-ranked pose compared against the crystallographic ligand.
#
# The system is 3Q47 and its native peptide npisdvd (ASN459-PRO460-ILE461-
# SER462-ASP463-VAL464-ASP465). Because ADCP docks into a fixed receptor frame
# and the grid box was built around that ligand, the docked pose and the native
# are already in the same coordinate frame -- so RMSD is a direct comparison and
# needs no superposition. That is also the standard definition of docking RMSD.
#
# Ranking mirrors scripts/runADCP.py:293-297: sort by the 'best target energy'
# line each run prints, and report kcal/mol as targetE * 0.59219.
#
# Usage: run_redock_validation.sh <adcp_bin> <staged_target> <workdir> \
#                                 [nruns] [steps] [rmsd_threshold]
# =============================================================================
set -eu

ADCP="${1:?Usage: $0 <adcp_bin> <staged_target> <workdir> [nruns] [steps] [rmsd]}"
STAGE="${2:?}"
WD="${3:?}"
NRUNS="${4:-16}"
STEPS="${5:-2500000}"
RMSD_MAX="${6:-2.5}"

if [ ! -f "${STAGE}/rigidReceptor.C.map" ] || [ ! -f "${STAGE}/native_pep.pdbqt" ]; then
	echo "SKIP: docking target not staged in ${STAGE}"
	exit 77
fi

rm -rf "${WD}"
mkdir -p "${WD}"
cp "${STAGE}"/rigidReceptor.*.map "${STAGE}/transpoints" "${STAGE}/con" \
   "${STAGE}/ramaprob.data" "${STAGE}/native_pep.pdbqt" "${WD}/"
cd "${WD}"

SEQ=npisdvd
OPTS='Bias=NULL,external=5,con,1.0,1.0,Opt=1,0.25,0.75,0.0'

echo "Redocking ${SEQ} into 3Q47: ${NRUNS} runs x ${STEPS} steps"

# Independent searches, in parallel. Seeds are fixed so the whole validation is
# reproducible rather than a different experiment on every invocation.
i=1
while [ "$i" -le "$NRUNS" ]; do
	( "$ADCP" "$SEQ" -r "1x${STEPS}" -p "$OPTS" -s "$i" -o "run_${i}.pdb" \
	    > "run_${i}.log" 2>&1; echo $? > "run_${i}.rc" ) &
	i=$((i + 1))
done
wait

# Any crash invalidates the whole result. A segfaulting ADCP still writes a
# plausible truncated PDB, so this must be checked explicitly.
FAILED=0
i=1
while [ "$i" -le "$NRUNS" ]; do
	rc=$(cat "run_${i}.rc" 2>/dev/null || echo 1)
	if [ "$rc" -ne 0 ]; then
		echo "FAIL: run $i exited $rc"
		FAILED=$((FAILED + 1))
	fi
	i=$((i + 1))
done
[ "$FAILED" -eq 0 ] || { echo "FAIL: ${FAILED}/${NRUNS} runs did not complete"; exit 1; }

python3 - "$NRUNS" "$RMSD_MAX" <<'PY'
import sys, glob, re

nruns = int(sys.argv[1])
rmsd_max = float(sys.argv[2])

def backbone(path, from_model=False):
    """Backbone atoms as {(resnum, name): (x,y,z)} plus residues in sequence order.

    Residues are sorted numerically, NOT taken in order of appearance: PDBQT
    writes atoms in torsion-tree order, so the native ligand's records start at
    residue 462 and jump around. Pairing by first appearance silently compares
    the wrong residues and inflates the RMSD.
    """
    lines = open(path).read().splitlines()
    if from_model:
        cur, keep = [], []
        for l in lines:
            if l.startswith('MODEL'):
                cur = []
            elif l.startswith('ATOM'):
                cur.append(l)
            elif l.startswith('ENDMDL'):
                keep = cur
        lines = keep or cur
    else:
        lines = [l for l in lines if l[:4] in ('ATOM', 'HETA')]
    out = {}
    for l in lines:
        name = l[12:16].strip()
        if name not in ('N', 'CA', 'C', 'O'):
            continue
        res = int(l[22:26])
        out[(res, name)] = (float(l[30:38]), float(l[38:46]), float(l[46:54]))
    return out, sorted({r for r, _ in out})

# Rank exactly as runADCP.py does: by the 'best target energy' each run reports.
energies = []
for i in range(1, nruns + 1):
    e = None
    for line in open('run_%d.log' % i):
        if line.startswith('best target energy'):
            e = float(line.split()[3])
    if e is None:
        sys.exit('FAIL: run %d never printed "best target energy"' % i)
    energies.append((e, i))
energies.sort()

print('  rank  run   targetE     kcal/mol')
for rank, (e, i) in enumerate(energies[:5], 1):
    print('  %4d  %3d  %9.4f  %9.2f' % (rank, i, e, e * 0.59219))

nat, nat_order = backbone('native_pep.pdbqt')
best_e, best_run = energies[0]
pose, pose_order = backbone('run_%d.pdb' % best_run, from_model=True)

if len(pose_order) != len(nat_order):
    sys.exit('FAIL: %d docked residues vs %d native' % (len(pose_order), len(nat_order)))

def rmsd(pose, pose_order):
    tot = n = 0
    for rn, rd in zip(nat_order, pose_order):
        for nm in ('N', 'CA', 'C', 'O'):
            a, b = nat.get((rn, nm)), pose.get((rd, nm))
            if a and b:
                tot += sum((x - y) ** 2 for x, y in zip(a, b))
                n += 1
    return (tot / n) ** 0.5, n

best_rmsd, natoms = rmsd(pose, pose_order)
print('  top-ranked pose: run %d, targetE %.4f (%.2f kcal/mol)'
      % (best_run, best_e, best_e * 0.59219))
print('  backbone RMSD to native: %.2f A over %d atoms' % (best_rmsd, natoms))

# Informational: the best RMSD anywhere in the ensemble. If this is good while
# the top-ranked pose is not, the search works and the scoring/ranking is what
# is failing -- a materially different diagnosis.
all_r = []
for _, i in energies:
    p, o = backbone('run_%d.pdb' % i, from_model=True)
    if len(o) == len(nat_order):
        all_r.append((rmsd(p, o)[0], i))
all_r.sort()
print('  best RMSD in ensemble: %.2f A (run %d)' % all_r[0])

if best_rmsd > rmsd_max:
    print('FAIL: top-ranked pose RMSD %.2f A exceeds %.2f A' % (best_rmsd, rmsd_max))
    sys.exit(1)
print('PASS: redocked to %.2f A (threshold %.2f A)' % (best_rmsd, rmsd_max))
PY
