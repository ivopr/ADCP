# Protocol for Comparison 3

Every command that produced a number in [3.md](3.md), in the order they were
run, so the whole comparison can be reproduced or extended. This protocol is
structurally identical to [2-protocol.md](2-protocol.md) — same tooling, same
steps, same order — because the point of re-running it is to see whether
anything *changed*, not to invent a new methodology. Deviations from
2-protocol.md are called out explicitly where they happen. Paths are as used
during this session; substitute your own tmp locations as needed.

## 0. Environment

```sh
date -u +"%Y-%m-%d %H:%M:%S UTC"   # 2026-08-30 14:52:01 UTC (session start)
gcc --version        # gcc (Ubuntu 13.3.0-6ubuntu2~24.04.1) 13.3.0
clang --version       # Ubuntu clang version 18.1.3 (1ubuntu1)
lscpu | grep -E "Model name|^CPU\(s\)"   # 12th Gen Intel Core i7-12650H, 16 threads
nproc                 # 16
git rev-parse HEAD    # 1fb8818b833e2f248e96e71e731f92e8ba0f9cb5
```

Same machine (same CPU model) as comparison 2. Pristine is the same commit,
`1c1a330` — it never changes. Current HEAD is `1fb8818`, 11 commits ahead of
comparison 2's `4ed5d6a` (`git log --oneline 4ed5d6a..1fb8818`); see 3.md for
what those 11 commits are and why none of them were expected to move any
number here.

**The docking target cache was gone this time** (`build-dock/tests/target`
had been deleted earlier in the session it turns out doesn't matter for
reproducing this protocol — just re-fetch it if it's missing):

```sh
sh tests/fetch_target.sh build-dock/tests/target
# Fetching https://ccsb.scripps.edu/mamba/examples/3Q47.trg (~10 MB) ...
# staged 9 maps and 106 translation points
```

**A `noclobber`-related gotcha for this specific shell setup**: this
session's zsh has `noclobber` set, so a plain `cmd > existing_file.log`
fails with a cryptic `file exists` error if the target already exists from a
prior attempt — it is not a real command failure, just zsh refusing to
clobber. Use `>|` (force-clobber) for any log redirect you might re-run, or
always redirect to a fresh filename. This cost some confusion mid-session
(see the MSan build below) before being identified.

**A `diff` gotcha**: this environment defines a shell function named `diff`
(from a Claude Code shell snapshot) that shadows `/usr/bin/diff` and breaks
on ordinary two-file diffs. Use `command diff` to bypass it.

## 1. Build pristine (`1c1a330`)

Use a git worktree so pristine and current HEAD coexist without disturbing
the main checkout:

```sh
git worktree add /tmp/<your-tmp>/adcp-pristine 1c1a330
cd /tmp/<your-tmp>/adcp-pristine
mkdir -p build_pristine build_pristine_tools
```

Pristine's tree is flat (no `src/`/`tools/` split) and predates CMake. Build
the main binary with the same optimization level CMake's `Release` preset
uses, plus `-fcommon` (required on gcc ≥10 to link pristine's tentative
global definitions):

```sh
gcc -std=c99 -O3 -DNDEBUG -fcommon \
  nested.c aadict.c energy.c main.c metropolis.c flex.c peptide.c probe.c \
  rotation.c vector.c params.c error.c checkpoint_io.c vdw.c canonicalAA.c \
  -lm -o build_pristine/adcp_pristine
```

Compiles with warnings only (implicit-declaration, unused fscanf/fgets
return values — all pre-existing, none new). Tool binaries link a different
subset of the same core files (pristine has no shared static library, so
each tool restates its own dependency list):

```sh
gcc -std=c99 -O3 -DNDEBUG -fcommon bfactor.c peptide.c aadict.c params.c \
  error.c vector.c rotation.c canonicalAA.c vdw.c -lm -o build_pristine_tools/bfactor

gcc -std=c99 -O3 -DNDEBUG -fcommon -include stdio.h cm.c aadict.c error.c \
  -lm -o build_pristine_tools/cm
  # -include stdio.h works around cm.c including params.h (which uses FILE*)
  # before stdio.h itself.

gcc -std=c99 -O3 -DNDEBUG -fcommon ramachandran.c peptide.c aadict.c params.c \
  error.c vector.c rotation.c canonicalAA.c vdw.c -lm -o build_pristine_tools/ramachandran

gcc -std=c99 -O3 -DNDEBUG -fcommon dssp2cm.c -lm -o build_pristine_tools/dssp2cm

gcc -std=c99 -O3 -DNDEBUG -fcommon cdlearn.c peptide.c aadict.c params.c \
  error.c vector.c rotation.c canonicalAA.c vdw.c energy.c probe.c \
  metropolis.c flex.c checkpoint_io.c nested.c -lm -o build_pristine_tools/cdlearn

gcc -std=c99 -O3 -DNDEBUG -fcommon pauling.c peptide.c aadict.c params.c \
  error.c vector.c rotation.c canonicalAA.c vdw.c -lm -o build_pristine_tools/pauling

gcc -std=c99 -O3 -DNDEBUG -fcommon mergie.c -lm -o build_pristine_tools/mergie
gcc -std=c99 -O3 -DNDEBUG -fcommon statistics.c -lm -o build_pristine_tools/statistics
gcc -std=c99 -O3 -DNDEBUG -fcommon oops.c -lm -o build_pristine_tools/oops
```

All 9 build clean.

## 2. Build current HEAD

```sh
cmake -B /tmp/<your-tmp>/adcp-current-build -DADCP_TOOLS=ON \
  -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTING=OFF <repo-root>
cmake --build /tmp/<your-tmp>/adcp-current-build -j"$(nproc)"
```

Fresh build directory — never reuse a stale tree across source renames (see
MIGRATION.md's own "Verification" section for why). Note the current tool
binary for Ramachandran is named `rama`, not `ramachandran` — the CMake
target name differs from pristine's source file name; keep that straight
when scripting comparisons across both trees.

## 3. CLI parity

```sh
<pristine>/build_pristine/adcp_pristine -h > h_pristine.txt 2>&1
<current>/src/adcp -h > h_current.txt 2>&1
command diff -u h_pristine.txt h_current.txt   # only the Usage: argv[0] path differs
```

## 4. Folding determinism + timing

```sh
/usr/bin/time -v sh tests/run_fold_test.sh <bin> <workdir> data/ramaprob.data
```

Run once per binary; the script runs the same seed twice internally and
asserts the two trajectories are byte-identical before reporting PASS.

**Beyond 2-protocol.md**: also diff the two binaries' fold output against
*each other*, not just each internally against itself — the script's
internal PASS only proves determinism, not cross-binary equality:

```sh
md5sum <pristine-workdir>/atoms1.txt <current-workdir>/atoms1.txt
command diff -q <pristine-workdir>/atoms1.txt <current-workdir>/atoms1.txt
```

Both md5s should equal `6a438d0a673006235fccd2b1b7007ba3` — MIGRATION.md's
own documented reference value for this fixture. Getting the exact same
hash as a value recorded in a totally different document, from an
independent fresh build, is a good sanity check that nothing about the
environment (compiler version, libm, etc.) has silently drifted.

## 5. Docking — 16-seed sweep at a fixed budget

`dock16.sh <bin> <label> <steps>` (reproduced below in full — this session
wrote it fresh rather than reusing a script from a prior session's `/tmp`,
since those don't survive between sessions):

```sh
#!/bin/sh
set -eu
BIN="${1:?Usage: $0 <bin> <label> <steps>}"; LABEL="${2:?}"; STEPS="${3:?}"
STAGE=<repo-root>/build-dock/tests/target/stage
WD=/tmp/<your-tmp>/dock16_${LABEL}_${STEPS}
rm -rf "$WD"; mkdir -p "$WD"
SEQ=npisdvd
OPTS='Bias=NULL,external=5,con,1.0,1.0,Opt=1,0.25,0.75,0.0'
echo "seed targetE md5 wall_s exit"
for seed in $(seq 1 16); do
	SD="$WD/seed$seed"; mkdir -p "$SD"
	cp "$STAGE"/rigidReceptor.*.map "$STAGE/transpoints" "$STAGE/con" "$STAGE/ramaprob.data" "$SD/"
	t0=$(date +%s.%N)
	set +e
	( cd "$SD" && "$BIN" "$SEQ" -r 1x"$STEPS" -p "$OPTS" -s "$seed" -o "seed${seed}.pdb" > run.log 2>&1 )
	rc=$?
	set -e
	t1=$(date +%s.%N)
	wall=$(awk -v a="$t0" -v b="$t1" 'BEGIN{printf "%.2f", b-a}')
	targetE=$(grep "best target energy" "$SD/run.log" 2>/dev/null | awk '{print $NF}')
	md5=$(grep '^ATOM' "$SD/seed${seed}.pdb" 2>/dev/null | md5sum | awk '{print $1}')
	echo "$seed ${targetE:-N/A} ${md5:-N/A} $wall $rc"
done
```

Run at `steps=10000` for both binaries. Result this session: 16/16 seeds
identical `targetE` *and* ATOM-record md5, matching 2.md's table exactly,
value for value — confirming the `sqrt`-precision fix still holds and
nothing since has moved this budget.

## 6. Crash-rate mapping

Same shape as `dock16.sh` but without `set -e` bracketing the docking call,
so a non-zero exit (SIGSEGV) is recorded, not treated as script-fatal:

```sh
#!/bin/sh
BIN="${1:?}"; LABEL="${2:?}"; STEPS="${3:?}"
STAGE=<repo-root>/build-dock/tests/target/stage
WD=/tmp/<your-tmp>/crashmap_${LABEL}_${STEPS}
rm -rf "$WD"; mkdir -p "$WD"
SEQ=npisdvd
OPTS='Bias=NULL,external=5,con,1.0,1.0,Opt=1,0.25,0.75,0.0'
crashes=0
for seed in $(seq 1 16); do
	SD="$WD/seed$seed"; mkdir -p "$SD"
	cp "$STAGE"/rigidReceptor.*.map "$STAGE/transpoints" "$STAGE/con" "$STAGE/ramaprob.data" "$SD/"
	set +e
	( cd "$SD" && "$BIN" "$SEQ" -r 1x"$STEPS" -p "$OPTS" -s "$seed" -o "seed${seed}.pdb" > run.log 2>&1 )
	rc=$?
	set -e
	[ "$rc" -ne 0 ] && crashes=$((crashes+1))
	echo "seed=$seed rc=$rc"
done
echo "CRASHES: $crashes/16 at steps=$STEPS ($LABEL)"
```

Run for pristine at `steps` = 50000, 150000, 250000; for current HEAD at
250000 only (0 crashes expected and confirmed — no need to sweep budgets
that can't fail). Every pristine number this session **exactly matched**
2.md's crash-rate table, seed-for-seed survival count, from a completely
independent fresh build — see 3.md's table for the numbers.

## 7. Seed-10 divergence (the pre-existing pristine UB, re-confirmed not re-derived)

Binary search over step count, both binaries, fixed seed 10, redocking
option set, comparing `best target energy`:

```sh
run_probe() {
  bin="$1"; steps="$2"
  WD=/tmp/<your-tmp>/seed10_probe_$$
  rm -rf "$WD"; mkdir -p "$WD"
  cp "$STAGE"/rigidReceptor.*.map "$STAGE/transpoints" "$STAGE/con" "$STAGE/ramaprob.data" "$WD/"
  ( cd "$WD" && "$bin" npisdvd -r 1x"$steps" \
      -p 'Bias=NULL,external=5,con,1.0,1.0,Opt=1,0.25,0.75,0.0' -s 10 -o out.pdb > run.log 2>&1 )
  grep "best target energy" "$WD/run.log" | awk '{print $NF}'
}
for steps in 15000 20000 25000 30000 35000 40000; do
  echo "$steps $(run_probe "$PRISTINE" $steps) $(run_probe "$CURRENT" $steps)"
done
```

Result: identical to 2.md — flip between 20,000 (`-11.9092`/`-11.9092`,
match) and 25,000 (`-15.3654`/`-11.9092`, diverge). **This session did not
re-derive the isolation** (patching pristine's `swapChains` array from
`[swapLength]` to `[swapLength+1]` alone, nothing else, to reproduce current
HEAD's value) — that mechanism was already nailed down in
[2-protocol.md](2-protocol.md) §7, and since pristine is a fixed, unchanging
commit, there is nothing for a re-derivation to find that the first
derivation didn't already establish. Re-run it only if you have reason to
doubt the original isolation, not as a matter of routine.

## 8. Full-length redocking protocol

```sh
{ time sh tests/run_redock_validation.sh <bin> build-dock/tests/target/stage \
    <workdir> 16 2500000 2.5 ; } > log 2>&1
```

Run once per binary. The script fans all 16 seeds out in parallel (`&` +
`wait`), so wall time is roughly one seed's worth of compute on a
16-thread machine, not 16×. Current HEAD: `PASS`, wall 1:39.02, top-ranked
pose run 15, targetE -28.5508 (-16.91 kcal/mol), RMSD 0.66 Å over 28 atoms —
matches MIGRATION.md's own documented redock result exactly. Pristine:
`FAIL`, wall 1:10.85 (short because most runs die early, not because it's
fast), 15/16 runs `exited 1`; per-seed `.log`/`.rc` files show 15
`Segmentation fault (core dumped)` and one clean completion (seed 10 — the
one seed the crash-rate map already showed surviving at this budget).

## 9. Nested sampling and checkpoint

```sh
sh tests/run_ns_test.sh <bin> <workdir> data/ramaprob.data tests/data/output.pdb
sh tests/run_checkpoint_test.sh <bin> <workdir> data/ramaprob.data tests/data/output.pdb
```

NS energy sequence reproduced **exactly** against MIGRATION.md's own
recorded before/after values (not just against 2.md — against the original
source of truth), on both sides:

- pristine: `40.058829 31.225769 10.764044 7.463047 6.393763`
- current: `40.058829 688.128472 81.995422 10.222738 7.086147`

This is `35fb3fb`'s deliberate NS point-seeding fix, unrelated to anything
in this session's own commits.

**Cross-binary checkpoint format check**: copy every `ckpt_*` file
pristine's run produced (not just `ckpt_0`) into a fresh directory with the
same `snapshots.pdb`/`ramaprob.data` fixtures, then run current HEAD's
`adcp` with `-R 2` pointed at them, using the exact `COMMON` options
`run_checkpoint_test.sh` uses internally (`grep '^COMMON=' tests/run_checkpoint_test.sh`
to get them verbatim rather than retyping from memory — they must match
exactly or the resume will legitimately fail on option mismatch, not format
mismatch):

```sh
cp <pristine-workdir>/ckpt_* <pristine-workdir>/snapshots.pdb <workdir>/
cp data/ramaprob.data <workdir>/
( cd <workdir> && <current>/src/adcp -n -f snapshots.pdb -r 2x10 -s 12345 \
    -p Bias=NULL -C 2,ckpt -R 2 -o cross.pdb )
```

`rc=0`, writes only `cross.pdb_3`/`cross.pdb_4` (not `_0`/`_1`/`_2`),
confirming it resumed past the checkpointed iteration rather than
restarting — format-compatible across the whole migration.

## 10. Tools/ correctness sweep

```sh
# byte-for-byte
<pristine>/cm    < tests/data/peptide12.pdb > out_p; <current>/cm    < tests/data/peptide12.pdb > out_c; command cmp out_p out_c

# ramachandran / rama -- note the binary name differs between trees
timeout 5 <pristine>/ramachandran < tests/data/peptide12.pdb   # rc=124, hangs
timeout 5 <current>/rama          < tests/data/peptide12.pdb   # rc=0, completes

# bfactor: multi-model fixture with mismatched residue counts, explicit timeout
timeout 10 <pristine>/bfactor < tests/data/bfactor_multimodel.pdb   # rc=124 THIS session
timeout 10 <current>/bfactor  < tests/data/bfactor_multimodel.pdb   # rc=1, clean stop()

# cdlearn: baseline, run from a directory that HAS ramaprob.data --
# current HEAD now checks for it before reaching the -L check (that's the
# fix from this session's own work, "cdlearn never calls
# ramaprob_initialise()"), so running from a directory without it produces
# a DIFFERENT (also correct, but not apples-to-apples) early exit. Put
# ramaprob.data next to the binary first.
cp data/ramaprob.data <somewhere>/
( cd <somewhere> && <pristine>/cdlearn < /dev/null )   # rc=1, "No list file given"
( cd <somewhere> && <current>/cdlearn  < /dev/null )   # rc=1, same message, plus "ramaprob initialise success"

# cdlearn: the long-token regression fixture
( cd <somewhere> && <pristine>/cdlearn -L tests/data/cdlearn_longtoken.txt -l V )
  # rc=134, glibc "*** buffer overflow detected ***", core dumped
( cd <somewhere> && <current>/cdlearn -L tests/data/cdlearn_longtoken.txt -l V )
  # rc=1, "ERROR! cdlearn: PDB base name from list file is too long."
```

**The `bfactor` result is more precise than 2-protocol.md's characterization
this time.** 2-protocol.md's session found the corruption "didn't visibly
crash on the plain build" and needed ASan/an explicit `strlen` probe to
surface it as a heap-buffer-overflow. This session's plain (non-ASan) run,
on the exact same fixture, hung indefinitely instead (`rc=124` under a
10-second `timeout`) — a different symptom of the same underlying undefined
behavior, both consistent with "reads past the end of an under-sized
buffer," neither more or less "correct" than the other since UB has no
required manifestation. Confirm the actual defect with ASan regardless of
which symptom you see first (below) — don't rely on the plain-build symptom
alone, since it can differ run to run, machine to machine, or (as here)
session to session on the identical binary and identical input.

**ASan confirmation for `bfactor`**:

```sh
gcc -std=c99 -O0 -g -fcommon -fsanitize=address,undefined \
  bfactor.c peptide.c aadict.c params.c error.c vector.c rotation.c \
  canonicalAA.c vdw.c -lm -o bfactor_asan
timeout 10 ./bfactor_asan < tests/data/bfactor_multimodel.pdb
```

Result: identical to 2-protocol.md — `heap-buffer-overflow`, `READ of size
8`, `#0 subtract vector.c:37`, `#1 update bfactor.c:71`.

## 11. MemorySanitizer, adcp itself

Same as 2-protocol.md's §10 MSan section, run fresh:

```sh
clang -std=c99 -O1 -g -fcommon -fsanitize=memory -fno-omit-frame-pointer \
  -Wno-implicit-function-declaration -Wno-return-type -Wno-error \
  nested.c aadict.c energy.c main.c metropolis.c flex.c peptide.c probe.c \
  rotation.c vector.c params.c error.c checkpoint_io.c vdw.c canonicalAA.c \
  -lm -o build_msan/adcp_msan
```

```sh
cmake -B build-msan-current -DADCP_TOOLS=OFF -DCMAKE_BUILD_TYPE=Debug \
  -DBUILD_TESTING=OFF -DCMAKE_CXX_COMPILER=clang++ \
  -DCMAKE_CXX_FLAGS="-fsanitize=memory -fno-omit-frame-pointer -O1" \
  -DCMAKE_EXE_LINKER_FLAGS="-fsanitize=memory" <repo-root>
cmake --build build-msan-current --target adcp
```

**Fold, current HEAD only, with the CI-proven workaround**:

```sh
MSAN_OPTIONS="halt_on_error=0:abort_on_error=0:exitcode=0" \
  timeout 120 <msan-current>/src/adcp APGVGV -r 1000x500 -t 2 -s 12345 \
  -p Bias=NULL -o fold.pdb
```

`rc=0`, exactly **1** warning: `biasmap_initialise`, `energy.cpp:222`, the
known glibc-`fopen`-interceptor false positive documented in
`.github/workflows/ci.yml`'s `msan` job. Nothing else.

**Docking, seed 10, redocking options, 25,000 steps, both binaries**:

```sh
MSAN_OPTIONS="halt_on_error=0:abort_on_error=0:exitcode=0" timeout 60 \
  <bin> npisdvd -r 1x25000 -p 'Bias=NULL,external=5,con,1.0,1.0,Opt=1,0.25,0.75,0.0' \
  -s 10 -o dock.pdb
```

Pristine: **1** warning — `crankshaft`, `metropolis.c:732` — the chi1
finding, reproduced at the exact documented source line. Current HEAD: **2**
warnings — `mark_aa_from_file` (`peptide.cpp:1996`) and `biasmap_initialise`
(`energy.cpp:222`), both the same known `fopen`-adjacent false positives as
the fold run — and critically **no `crankshaft` report**. This is the one
number in this whole comparison that *should* differ from 2.md (2.md's own
session found the `crankshaft` report on current HEAD, since it predated
the fix): confirms `5dd2a23` (`fix: transmutate/transmove/transopt never
copied chi1/chi2, corrupting docking`) actually closed it, independently,
via the same tool that originally found it.

Don't expect the exact same warning *count* every run — this session got 1
report on pristine where 2-protocol.md's session got 3 before terminating
(MSan's continue-past-error is best-effort, and process-layout sensitivity
affects which of several genuinely-present issues get hit before exit). The
`crankshaft` report's *presence or absence*, not the total count, is the
signal that matters.

## 12. Compile-time and LOC

```sh
wc -l <pristine>/*.c <pristine>/*.h
wc -l src/*.cpp include/*.h tools/*.cpp

rm -rf build_timing && /usr/bin/time -v gcc -std=c99 -O3 -DNDEBUG -fcommon \
  <15 pristine sources> -lm -o build_timing/adcp

cmake -B adcp-current-build-timing -DADCP_TOOLS=OFF \
  -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTING=OFF <repo-root>
/usr/bin/time -v cmake --build adcp-current-build-timing --target adcp -j1
```

Both single-threaded, so the numbers are compiler-bound rather than
parallelism-bound; the CMake figure includes CMake's own generator/build
overhead, which the raw `gcc` figure has no equivalent of.

Also grabbed this session, cheap and not in 2-protocol.md's numbered steps:
binary size via `ls -la` on both `adcp` binaries.

## 13. Final sanity check

Not part of the pristine-vs-current comparison (pristine has no CMake/ctest
infrastructure to compare against), but worth doing before calling the
session done: confirm current HEAD's own test suite is still green after
all this — if the comparison somehow perturbed the working build, better to
find out now.

```sh
ctest --test-dir build --output-on-failure
```

## Cleanup

```sh
git worktree remove --force /tmp/<your-tmp>/adcp-pristine
git worktree prune
rm -rf /tmp/<your-tmp>/adcp-current-build /tmp/<your-tmp>/adcp-current-build-timing \
       /tmp/<your-tmp>/build-msan-current /tmp/<your-tmp>/dock16_* /tmp/<your-tmp>/crashmap_* \
       /tmp/<your-tmp>/seed10_probe_* /tmp/<your-tmp>/redock_* /tmp/<your-tmp>/ns_* \
       /tmp/<your-tmp>/ckpt_* /tmp/<your-tmp>/msan_*
```

The docking target cache (`build-dock/tests/target`) is worth keeping —
re-fetching it costs a network round-trip and ~10 MB every time it's
missing, which happened once already this session.
