#!/bin/sh
# =============================================================================
# ADCP -o overflow regression test.
#
# nestedsampling() builds its output filename with sprintf into a 1010-byte
# malloc'd buffer ("%s_%d", outfile_name, checkpoint_counter). outfile_name
# comes from -o via copy_string, which is unbounded, so a long enough -o value
# overflowed that buffer. Fixed with snprintf + a stop() bounds check
# (src/nested.cpp, the unconditional output-file-init block near the top of
# nestedsampling()) -- this asserts the fix, not the overflow: the process
# must stop() with a clear message instead of crashing.
#
# -o itself must survive main.cpp's own fopen(argv[i],"w") check, so a single
# 1000+ byte path component won't do (ENAMETOOLONG, NAME_MAX=255 on ext4) --
# build a real, deeply-nested but per-component-short directory tree instead,
# the way an unusually deep working-directory tree could occur legitimately.
#
# Usage: run_checkpoint_longname_test.sh <adcp_binary> <workdir> <ramaprob.data> <multimodel.pdb>
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

# 21 levels x 51 bytes ("/" + 50 chars) = 1071 bytes of directory path, each
# component well under NAME_MAX -- comfortably past the 1010-byte buffer.
COMPONENT=$(printf 'x%.0s' $(seq 1 50))
DEEP="."
i=1
while [ "$i" -le 21 ]; do
	DEEP="$DEEP/$COMPONENT"
	i=$((i + 1))
done
mkdir -p "$DEEP"
LONGNAME="$DEEP/out"

set +e
OUT=$("$ADCP" -n -f snapshots.pdb -r 1x1 -s 12345 -p Bias=NULL -o "$LONGNAME" 2>&1)
STATUS=$?
set -e

if [ "$STATUS" -eq 0 ]; then
	echo "FAIL: adcp exited 0 with an overlong -o value; expected stop()"
	echo "$OUT"
	exit 1
fi

if ! echo "$OUT" | grep -q "too long to build a checkpoint output name"; then
	echo "FAIL: expected the -o-too-long stop() message, got:"
	echo "$OUT"
	exit 1
fi

echo "PASS: overlong -o rejected by stop() instead of overflowing"
