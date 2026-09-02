#!/bin/sh
# =============================================================================
# Prepare one ADCP docking target from an RCSB entry.
#
# Produces the same thing tests/fetch_target.sh gets pre-baked for 3Q47, but for
# an arbitrary PDB id: a .trg (AutoGrid maps + AutoSite translation points) plus
# the crystallographic peptide to score against.
#
# The box follows both ADCP papers: smallest box around the crystallographic
# peptide, 4 A padding on every side (Bioinformatics 2019 p.372; JCTC 2019
# "The docking box was defined using the crystallographic peptide with a padding
# of 4 A on every side").
#
# Targets built here are NOT byte-identical to the ones the authors used:
# AutoSite's fill-point count depends on the ADFRsuite build and platform (see
# the note in tests/CMakeLists.txt, where this machine finds 104 points for 3Q47
# and the published target has 106). provenance.txt records what this machine
# produced so a later discrepancy is traceable rather than mysterious.
#
# Usage: prepare_target.sh <pdb_id> <expected_peptide_length> <out_dir>
# Exit:  0  prepared
#        77 unavailable (no ADFRsuite / no network) -- CTest treats this as SKIP
#        1  prepared nothing, and says why
# =============================================================================
set -eu

PDB_ID="${1:?Usage: $0 <pdb_id> <expected_peptide_length> <out_dir>}"
WANT_LEN="${2:?}"
OUT="${3:?}"

ADFR="${ADFR_HOME:-$HOME/ADFRsuite-1.0}"
PADDING="${ADCP_BOX_PADDING:-4}"
SCRIPT_DIR=$(cd "$(dirname "$0")" && pwd)

[ -x "${ADFR}/bin/agfr" ] || { echo "SKIP: no ADFRsuite at ${ADFR}"; exit 77; }
command -v python3 > /dev/null 2>&1 || { echo "SKIP: python3 not available"; exit 77; }

ID=$(printf '%s' "$PDB_ID" | tr 'A-Z' 'a-z')
mkdir -p "$OUT"
cd "$OUT"

# --- 1. fetch --------------------------------------------------------------
if [ ! -s "${ID}.pdb" ]; then
	URL="https://files.rcsb.org/download/$(printf '%s' "$ID" | tr 'a-z' 'A-Z').pdb"
	if command -v curl > /dev/null 2>&1; then
		curl -sSfL --connect-timeout 20 --max-time 300 -o "${ID}.pdb.part" "$URL" || {
			rm -f "${ID}.pdb.part"; echo "SKIP: could not download ${ID} (offline?)"; exit 77; }
	elif command -v wget > /dev/null 2>&1; then
		wget -q --timeout=20 -O "${ID}.pdb.part" "$URL" || {
			rm -f "${ID}.pdb.part"; echo "SKIP: could not download ${ID} (offline?)"; exit 77; }
	else
		echo "SKIP: neither curl nor wget"; exit 77
	fi
	mv "${ID}.pdb.part" "${ID}.pdb"
fi

# --- 2. split receptor / peptide -------------------------------------------
# The peptide chain is identified by residue count, then checked against the
# length the paper published. Getting this wrong silently would compare against
# the wrong molecule for the rest of the benchmark's life, so a mismatch is
# fatal and prints every chain it saw.
python3 - "${ID}.pdb" "$WANT_LEN" "$PADDING" <<'PY'
import sys, collections

path, want = sys.argv[1], int(sys.argv[2])
AA3 = set("ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR "
          "TRP TYR VAL MSE HID HIE HIP CYX SEC PYL".split())

model, chains = 0, collections.OrderedDict()
for line in open(path):
    if line.startswith('ENDMDL'):
        break                      # biological unit: first model only
    if not line.startswith('ATOM') and not line.startswith('HETATM'):
        continue
    alt = line[16]
    if alt not in (' ', 'A'):
        continue                   # keep one altloc, the first
    resn = line[17:20].strip()
    ch = line[21]
    key = (line[22:27], resn)      # resseq + icode
    rec = chains.setdefault(ch, {'res': [], 'std': 0, 'lines': []})
    if key not in rec['res']:
        rec['res'].append(key)
        if resn in AA3:
            rec['std'] += 1
    rec['lines'].append(line)

report = ', '.join('%s:%d(%d std)' % (c, len(v['res']), v['std'])
                   for c, v in chains.items())

# Prefer the chain whose length matches the published one, then fall back to the
# shortest peptide-sized chain.
#
# The published length is NOT always the chain's residue count. JCTC 2019 Table 1
# lists the group-I HIV integrase binders as length 6, but 3AV9 chain X is the
# 8-mer SAKIDNLD -- the paper appears to be counting the macrocycle, not the
# chain. So a mismatch is reported and recorded, never fatal: dying here would
# reject a target over a bookkeeping difference in someone else's table.
peptide_like = [c for c, v in chains.items()
                if 0 < v['std'] <= 30 and len(chains) > 1]
cand = [c for c in peptide_like if chains[c]['std'] == want]
if not cand:
    cand = sorted(peptide_like, key=lambda c: (chains[c]['std'], c))
if not cand:
    sys.exit('FAIL: no peptide-sized chain (<=30 standard residues). chains: %s'
             % report)
if len(cand) > 1 and chains[cand[0]]['std'] == chains[cand[1]]['std']:
    sys.stderr.write('NOTE: %d chains of equal length (%s); taking %s\n'
                     % (len(cand), ','.join(cand), cand[0]))
pep = cand[0]
obs = chains[pep]['std']
if obs != want:
    sys.stderr.write('NOTE: chain %s has %d standard residues, manifest says %d '
                     '-- using the structure, recording both\n' % (pep, obs, want))

# Blank the altLoc column on the conformer we keep. Leaving the 'A' in place
# makes prepare_receptor derive 5-character hydrogen names from it (HE2A -> the
# H names HE2A1/HE2A2), which overflow the atom-name field and push the residue
# name one column left. agfr then dies with "invalid or missing coordinate(s)"
# hundreds of lines from the cause. This is what prepare_receptor's own
# deleteAltB cleanup means by "rename XX@A atoms -> XX".
def one_conformer(l):
    return l[:16] + ' ' + l[17:]

# Only the standard amino acids of that chain are the peptide. Waters and ions
# sharing its chain id are not, and letting them through inflates the docking
# box -- 3Q47 carries 9 such residues on the peptide chain.
pep_lines = [one_conformer(l) for l in chains[pep]['lines']
             if l[17:20].strip() in AA3]
rec_lines = [one_conformer(l) for c, v in chains.items() if c != pep
             for l in v['lines']]
if not rec_lines:
    sys.exit('FAIL: chain %s is the whole structure, nothing left as receptor' % pep)

with open('peptide.pdb', 'w') as fh:
    fh.writelines(pep_lines); fh.write('END\n')
with open('receptor_raw.pdb', 'w') as fh:
    fh.writelines(rec_lines); fh.write('END\n')

ONE = {'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLN':'Q','GLU':'E',
       'GLY':'G','HIS':'H','ILE':'I','LEU':'L','LYS':'K','MET':'M','MSE':'M',
       'PHE':'F','PRO':'P','SER':'S','THR':'T','TRP':'W','TYR':'Y','VAL':'V'}
seq = ''.join(ONE.get(r[1], 'X') for r in chains[pep]['res'] if r[1] in AA3)
with open('peptide.seq', 'w') as fh:
    fh.write(seq + '\n')
with open('peptide.chain', 'w') as fh:
    fh.write('%s %d %d\n' % (pep, obs, want))   # chain, observed len, published len

# The box is computed here rather than handed to agfr as "-b ligand": that mode
# wants a ligand PDBQT with a torsion tree, and building one for a 20-mer means
# running the peptide through prepare_ligand's rotatable-bond detection for no
# benefit. The box is a bounding box plus padding either way, so compute it.
xs = [float(l[30:38]) for l in pep_lines]
ys = [float(l[38:46]) for l in pep_lines]
zs = [float(l[46:54]) for l in pep_lines]
pad = float(sys.argv[3])
box = []
for v in (xs, ys, zs):
    box.append(((min(v) + max(v)) / 2.0, (max(v) - min(v)) + 2 * pad))
with open('box.txt', 'w') as fh:
    fh.write('%.3f %.3f %.3f %.3f %.3f %.3f\n'
             % (box[0][0], box[1][0], box[2][0], box[0][1], box[1][1], box[2][1]))

print('peptide chain %s, %d residues, sequence %s' % (pep, want, seq))
print('chains: %s' % report)
print('box center %.3f %.3f %.3f  size %.3f %.3f %.3f'
      % (box[0][0], box[1][0], box[2][0], box[0][1], box[1][1], box[2][1]))
PY

# --- 3. receptor -> pdbqt ---------------------------------------------------
# prepare_receptor's own defaults do the cleanup (-U nphs_lps_waters_nonstdres):
# strip waters, merge non-polar hydrogens, drop all-nonstandard chains.
#
# Two passes. prepare_receptor does not fail on a heteroatom it cannot type -- it
# prints a warning and writes a PDBQT line with no charge, which agfr then
# rejects with a coordinate parse error hundreds of lines away from the cause
# (3WNF's cadmium ions do exactly this). So: run it, ask its own log which
# residues it could not type, drop those, run it again. Using the tool's
# diagnosis beats maintaining a guessed whitelist of acceptable metals -- and the
# dropped atoms have no AutoGrid map anyway, so they could not contribute to the
# docking even if they survived.
prep() {
	"${ADFR}/bin/prepare_receptor" -r "$1" -A checkhydrogens \
		-o receptor.pdbqt > prepare_receptor.log 2>&1
}
prep receptor_raw.pdb || {
	echo "FAIL: prepare_receptor failed for ${ID}"; tail -5 prepare_receptor.log; exit 1; }

UNTYPED=$(sed -n 's/.*no Gasteiger parameters available for atom [^:]*:[^:]*: *\([A-Za-z]\{1,3\}\)[0-9]*:.*/\1/p' \
	prepare_receptor.log | sort -u | tr '\n' ' ')
if [ -n "$(printf '%s' "$UNTYPED" | tr -d ' ')" ]; then
	echo "NOTE: dropping untypeable hetero residues: ${UNTYPED}"
	python3 - receptor_raw.pdb receptor_clean.pdb $UNTYPED <<'DROPPY'
import sys
src, dst, drop = sys.argv[1], sys.argv[2], {r.upper() for r in sys.argv[3:]}
kept = 0
with open(src) as fi, open(dst, 'w') as fo:
    for l in fi:
        if l.startswith('HETATM') and l[17:20].strip().upper() in drop:
            continue
        fo.write(l)
        kept += l.startswith(('ATOM', 'HETATM'))
print('  receptor kept %d atoms' % kept)
DROPPY
	mv receptor_clean.pdb receptor_raw.pdb
	prep receptor_raw.pdb || {
		echo "FAIL: prepare_receptor failed for ${ID} after dropping ${UNTYPED}"
		tail -5 prepare_receptor.log; exit 1; }
fi
echo "$UNTYPED" > receptor_dropped_het.txt

# --- 4. maps ----------------------------------------------------------------
AGFR_CMD="${ADFR}/bin/agfr -r receptor.pdbqt -b user $(cat box.txt) -m all -o ${ID}"
$AGFR_CMD > agfr.log 2>&1 || { echo "FAIL: agfr failed for ${ID}"; tail -20 agfr.log; exit 1; }
[ -s "${ID}.trg" ] || { echo "FAIL: agfr produced no ${ID}.trg"; tail -20 agfr.log; exit 1; }

# --- 5. stage ---------------------------------------------------------------
# Unpack the .trg into the flat layout the engine expects in its working
# directory -- the same nine map types scripts/runADCP.py stages. Doing it once
# here means run_set.sh only has to copy a directory.
python3 - "${ID}.trg" stage <<'STAGEPY'
import sys, os, struct, zipfile

trg, stage = sys.argv[1], sys.argv[2]
ELEMENTS = ['C', 'A', 'SA', 'N', 'NA', 'OA', 'HD', 'd', 'e']
os.makedirs(stage, exist_ok=True)

with zipfile.ZipFile(trg) as z:
    names = z.namelist()
    root = os.path.commonprefix([n for n in names if '/' in n]).split('/')[0]
    for e in ELEMENTS:
        member = '%s/rigidReceptor.%s.map' % (root, e)
        if member not in names:
            sys.exit('missing %s in %s' % (member, trg))
        open(os.path.join(stage, 'rigidReceptor.%s.map' % e), 'wb').write(z.read(member))
    npy = z.read('%s/translationPoints.npy' % root)

# Dtype is NOT fixed across ADFRsuite builds: the published 3Q47 target carries
# float64 points (which is all tests/fetch_target.sh knows how to read) while
# the agfr on this machine emits float32. Read whichever is there.
if npy[:6] != b'\x93NUMPY':
    sys.exit('translationPoints.npy is not a .npy file')
hlen = struct.unpack('<H', npy[8:10])[0]
header = npy[10:10 + hlen].decode()
if 'True' in header.split('fortran_order')[1][:12]:
    sys.exit('fortran-ordered translationPoints.npy is not supported: %s' % header.strip())
if "'<f8'" in header:
    code, width = 'd', 8
elif "'<f4'" in header:
    code, width = 'f', 4
else:
    sys.exit('unexpected translationPoints.npy dtype: %s' % header.strip())
body = npy[10 + hlen:]
vals = struct.unpack('<%d%s' % (len(body) // width, code), body)
pts = [vals[i:i + 3] for i in range(0, len(vals), 3)]

with open(os.path.join(stage, 'transpoints'), 'w') as fh:
    fh.write('%s\n' % len(pts))
    for p in pts:
        fh.write('%7.3f %7.3f %7.3f\n' % p)

# peptide.cpp stop()s without this file; type-5 docking scores every residue anyway
open(os.path.join(stage, 'con'), 'w').write('1\n')
open('translation_points.count', 'w').write('%d\n' % len(pts))
print('staged %d maps and %d translation points' % (len(ELEMENTS), len(pts)))
STAGEPY

cp peptide.pdb stage/native_pep.pdb
RAMA="${SCRIPT_DIR}/../../data/ramaprob.data"
[ -f "$RAMA" ] && cp "$RAMA" stage/ramaprob.data

# Every staged map is checked against its own header before this target is
# declared usable. agfr does not always fail cleanly: it can die inside
# addGradientToMaps *after* writing a .trg, leaving maps that are truncated or
# padded. Testing "the .trg exists and is non-empty" accepted 6 of 49 targets
# whose maps disagreed with their own NELEMENTS by thousands of values. A
# corrupt map does not crash the docking -- it silently changes the potential,
# which is the worst possible failure mode for a benchmark.
python3 - stage <<'CHECKPY'
import sys, os
stage = sys.argv[1]
ELEM = ['C', 'A', 'SA', 'N', 'NA', 'OA', 'HD', 'd', 'e']
bad = []
for e in ELEM:
    p = os.path.join(stage, 'rigidReceptor.%s.map' % e)
    if not os.path.exists(p):
        bad.append('%s.map missing' % e); continue
    with open(p) as fh:
        head = [next(fh, '') for _ in range(6)]
        try:
            nx, ny, nz = (int(x) for x in head[4].split()[1:4])
        except Exception:
            bad.append('%s.map has an unreadable NELEMENTS header' % e); continue
        want = (nx + 1) * (ny + 1) * (nz + 1)
        got = sum(1 for _ in fh)
    if got != want:
        bad.append('%s.map holds %d values, its header declares %d' % (e, got, want))
if bad:
    sys.stderr.write('FAIL: agfr produced corrupt maps:\n')
    for b in bad:
        sys.stderr.write('  %s\n' % b)
    sys.exit(1)
CHECKPY
if [ $? -ne 0 ]; then
	rm -rf stage "${ID}.trg"
	echo "FAIL: ${ID} maps did not validate; stage and .trg removed so a partial"
	echo "      target cannot be mistaken for a good one. See agfr.log."
	exit 1
fi

# --- 6. provenance ----------------------------------------------------------
# Everything a later reader needs to explain a number that does not reproduce.
NPTS=$(cat translation_points.count 2>/dev/null || echo unknown)
{
	echo "pdb_id            ${ID}"
	echo "peptide_length    published=${WANT_LEN} observed=$(cut -d' ' -f2 peptide.chain)"
	echo "peptide_chain     $(cut -d' ' -f1 peptide.chain)"
	echo "peptide_sequence  $(cat peptide.seq)"
	echo "dropped_het       $(cat receptor_dropped_het.txt 2>/dev/null || echo none)"
	echo "box               ligand + ${PADDING} A padding"
	echo "translation_pts   ${NPTS}"
	echo "adfrsuite         ${ADFR}"
	echo "agfr_version      $(${ADFR}/bin/agfr --version 2>&1 | grep -iEo '[0-9]+\.[0-9.]+' | head -1 || echo unknown)"
	echo "agfr_cmd          ${AGFR_CMD}"
	echo "prepared_utc      $(date -u +%Y-%m-%dT%H:%M:%SZ)"
} > provenance.txt

echo "OK: ${ID}.trg ready (${NPTS} translation points, peptide $(cat peptide.seq))"
