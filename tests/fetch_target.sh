#!/bin/sh
# =============================================================================
# Fetch and stage the published 3Q47 docking target.
#
# ADCP's docking path needs AutoGrid maps that normally come from AGFR (part of
# the ADFRsuite). The ADCP documentation publishes an already-prepared target
# for 3Q47, so the suite is not required:
#
#   https://ccsb.scripps.edu/adcpv11/documentation/
#
# The .trg is a zip. Staging mirrors scripts/runADCP.py:143-171 exactly -- same
# map list, same transpoints format, same 'con' file -- so the tests exercise
# the real production layout rather than an invented one.
#
# translationPoints.npy is read with struct rather than numpy: it is a plain
# float64 array behind a short header, and requiring numpy just to run the test
# suite would be a needless dependency.
#
# Usage: fetch_target.sh <cache_dir>
# Exit:  0  staged successfully
#        77 unavailable (no network / no python3) -- CTest treats this as SKIP
# =============================================================================
set -eu

CACHE="${1:?Usage: $0 <cache_dir>}"
TRG_URL="https://ccsb.scripps.edu/mamba/examples/3Q47.trg"
PEP_URL="https://ccsb.scripps.edu/adcpv11/download/471/"
STAGE="${CACHE}/stage"
TRG="${CACHE}/3Q47.trg"
PEP="${CACHE}/3Q47_pep.pdbqt"

mkdir -p "${CACHE}"

if ! command -v python3 > /dev/null 2>&1; then
	echo "SKIP: python3 not available, cannot unpack the target"
	exit 77
fi

fetch() {
	# $1 url, $2 destination
	if command -v curl > /dev/null 2>&1; then
		curl -sSfL --connect-timeout 20 --max-time 600 -o "$2" "$1"
	elif command -v wget > /dev/null 2>&1; then
		wget -q --timeout=20 -O "$2" "$1"
	else
		return 1
	fi
}

if [ ! -s "${TRG}" ]; then
	echo "Fetching ${TRG_URL} (~10 MB) ..."
	if ! fetch "${TRG_URL}" "${TRG}.part"; then
		rm -f "${TRG}.part"
		echo "SKIP: could not download the 3Q47 target (offline?)"
		exit 77
	fi
	mv "${TRG}.part" "${TRG}"
else
	echo "Using cached ${TRG}"
fi

if [ ! -s "${PEP}" ]; then
	if ! fetch "${PEP_URL}" "${PEP}.part"; then
		rm -f "${PEP}.part"
		echo "SKIP: could not download the native 3Q47 peptide (offline?)"
		exit 77
	fi
	mv "${PEP}.part" "${PEP}"
fi

rm -rf "${STAGE}"
mkdir -p "${STAGE}"

python3 - "${TRG}" "${STAGE}" <<'PY'
import sys, os, struct, zipfile

trg, stage = sys.argv[1], sys.argv[2]

# Same nine atom types runADCP.py stages. main.c picks a subset based on the
# sequence (SA for CYS, A for aromatics, NA for HIS) and falls back to C.map
# for the rest, so staging all nine keeps any test sequence working.
ELEMENTS = ['C', 'A', 'SA', 'N', 'NA', 'OA', 'HD', 'd', 'e']

with zipfile.ZipFile(trg) as z:
    names = z.namelist()
    root = os.path.commonprefix([n for n in names if '/' in n]).split('/')[0]

    for e in ELEMENTS:
        member = '%s/rigidReceptor.%s.map' % (root, e)
        if member not in names:
            sys.exit('missing %s in target archive' % member)
        with open(os.path.join(stage, 'rigidReceptor.%s.map' % e), 'wb') as fh:
            fh.write(z.read(member))

    npy = z.read('%s/translationPoints.npy' % root)
    with open(os.path.join(stage, '3Q47_recH.pdbqt'), 'wb') as fh:
        fh.write(z.read('%s/3Q47_recH.pdbqt' % root))

# .npy: b'\x93NUMPY', major, minor, uint16 header length, then raw little-endian
# float64. Guard the assumptions rather than trusting the file blindly.
if npy[:6] != b'\x93NUMPY':
    sys.exit('translationPoints.npy is not a .npy file')
hlen = struct.unpack('<H', npy[8:10])[0]
header = npy[10:10 + hlen].decode()
if "'<f8'" not in header or 'True' in header.split('fortran_order')[1][:12]:
    sys.exit('unexpected translationPoints.npy layout: %s' % header.strip())
body = npy[10 + hlen:]
vals = struct.unpack('<%dd' % (len(body) // 8), body)
pts = [vals[i:i + 3] for i in range(0, len(vals), 3)]

# Byte-for-byte the format runADCP.py writes: count line, then '%7.3f' columns.
with open(os.path.join(stage, 'transpoints'), 'w') as fh:
    fh.write('%s\n' % len(pts))
    for p in pts:
        fh.write('%7.3f %7.3f %7.3f\n' % p)

# runADCP.py writes a constraint file containing just '1'. The file must exist
# (peptide.c stop()s otherwise) but type-5 docking scores every residue anyway.
with open(os.path.join(stage, 'con'), 'w') as fh:
    fh.write('1\n')

print('staged %d maps and %d translation points' % (len(ELEMENTS), len(pts)))
PY

cp "${PEP}" "${STAGE}/native_pep.pdbqt"

# Required for both folding and docking; energy.c looks for it in the CWD.
RAMA_SRC="$(dirname "$0")/../data/ramaprob.data"
[ -f "${RAMA_SRC}" ] || RAMA_SRC="${2:-}/ramaprob.data"
cp "${RAMA_SRC}" "${STAGE}/ramaprob.data"

echo "Target staged in ${STAGE}"
