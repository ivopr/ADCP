# Protocol for Comparison 2

Every command that produced a number in [2.md](2.md), in the order they were
run, so the whole comparison can be reproduced or extended. Paths are as used
during the session; substitute your own tmp locations as needed.

## 0. Environment

```
gcc --version        # gcc (Ubuntu 13.3.0-6ubuntu2~24.04.1) 13.3.0
lscpu                # 12th Gen Intel Core i7-12650H, 16 threads
nproc                # 16
```

The docking target (3Q47) was already staged at
`build-dock/tests/target/stage` from prior work in this repo — no network
fetch was needed. If it isn't present, run
`tests/fetch_target.sh <cache_dir>` first.

## 1. Build pristine (`1c1a330`)

```sh
git worktree add /tmp/adcp-pristine 1c1a330
cd /tmp/adcp-pristine
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

Tool binaries link a different subset of the same core files (pristine has
no shared static library, so each tool restates its own dependency list;
determined by resolving each linker error one at a time):

```sh
gcc -std=c99 -O3 -DNDEBUG -fcommon bfactor.c peptide.c aadict.c params.c \
  error.c vector.c rotation.c canonicalAA.c vdw.c -lm -o build_pristine_tools/bfactor

gcc -std=c99 -O3 -DNDEBUG -fcommon -include stdio.h cm.c aadict.c error.c \
  -lm -o build_pristine_tools/cm
  # -include stdio.h works around cm.c including params.h (which uses FILE*)
  # before stdio.h itself -- a real pristine header-ordering bug, harmless
  # here since it only needed forcing the declaration to appear first.

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

## 2. Build current HEAD

Fresh build directory (never reuse a stale tree across source renames — see
MIGRATION.md's own "Verification" section for why):

```sh
cmake -B /tmp/adcp-current-build -DADCP_TOOLS=ON -DCMAKE_BUILD_TYPE=Release \
  -DBUILD_TESTING=OFF .
cmake --build /tmp/adcp-current-build -j"$(nproc)"
```

## 3. CLI parity

```sh
/tmp/adcp-pristine/build_pristine/adcp_pristine -h > h_pristine.txt 2>&1
/tmp/adcp-current-build/src/adcp -h > h_current.txt 2>&1
diff -u h_pristine.txt h_current.txt   # only argv[0] differs
```

## 4. Folding determinism + timing

```sh
/usr/bin/time -v sh tests/run_fold_test.sh <bin> <workdir> data/ramaprob.data
```
Run once per binary; the script itself runs the same seed twice internally
and asserts the two trajectories are byte-identical before reporting PASS.

## 5. Docking — 16-seed sweep at a fixed budget

`/tmp/adcp-compare/dock16.sh <bin> <label> <steps>` (written for this
session): stages the 3Q47 target fresh per seed, runs
`<bin> npisdvd -r 1x<steps> -p 'Bias=NULL,external=5,con,1.0,1.0,Opt=1,0.25,0.75,0.0' -s <seed> -o seed<seed>.pdb`
for seeds 1–16, and records `best target energy`, an md5 of the `ATOM`
records, and per-run wall time via `date +%s.%N` deltas. Used at
`steps=10000` for the full 16-seed table in 2.md.

## 6. Crash-rate mapping

`/tmp/adcp-compare/crashmap.sh <bin> <steps>`: same options and sequence,
16 seeds, records exit code per seed without treating a non-zero exit as a
script-fatal error (`set +e`/`set -e` bracketing the docking call — a `set
-e` script would otherwise abort at the first SIGSEGV instead of recording
it). Run at steps = 50000, 150000, 250000 for pristine to reproduce the
crash-rate table.

## 7. Seed-10 divergence bisection (the headline finding)

Binary search over step count, both binaries, fixed seed 10, same option
set, comparing `best target energy`:

```sh
for steps in 15000 20000 25000 30000 35000 40000; do
  # run both binaries at $steps, seed 10, diff targetE
done
```

Found the flip between 20,000 (match) and 25,000 (diverge). Then `git
bisect` over commit range `1c1a330..a1b6718` (the commit compares/1.md's
"0/16 differ at 2.5M steps" claim was measured against — first confirmed
byte-identical to current HEAD's `adcp` binary via `cmp`, ruling out any of
the 14 commits after it as the cause) using a probe script
(`/tmp/adcp-compare/bisect_probe.sh`) that:

1. `cd`s into a bisect worktree, detects CMake vs flat-Makefile-era layout
   by the presence of `CMakeLists.txt`, and builds `adcp` accordingly
   (`exit 125` to skip a revision that won't build under either scheme —
   used for `efa13ff`, the one commit still mid-restructure).
2. Runs the seed-10/25,000-step reproduction in a fresh per-run directory.
3. `exit 0` if `best target energy` is pristine's `-15.3654`, `exit 1`
   otherwise.

First bad commit: `7e2192d` (confirmed by testing all 9 intermediate
commits individually rather than trusting `git bisect run`'s commit
selection blindly, since the range spans a C→CMake layout transition that
could confuse the automatic search).

**Isolation step**, once `7e2192d` was identified: copied pristine's
`main.c`, changed only

```diff
- Chain* swapChains[swapLength];
+ Chain* swapChains[swapLength + 1];
```

recompiled with the exact same pristine command line from step 1, reran the
seed-10/25,000-step probe. Result: `-11.9092`, matching current HEAD exactly
— confirming the array-size change alone (not `7e2192d`'s two `energy.c`
`NULL`-check additions, which never fire when input files are present) is
sufficient to reproduce the divergence.

## 8. Full-length redocking protocol

```sh
sh tests/run_redock_validation.sh <bin> build-dock/tests/target/stage \
  <workdir> 16 2500000 2.5
```
wrapped in `{ time ... ; } > log 2>&1` for wall/user/sys time. Run once per
binary. Pristine's run fails the script's own completion check (expected —
this is what it's measuring); the per-seed `.log` files under `<workdir>`
were inspected directly for `Segmentation fault (core dumped)` and the
`best target energy` line where present.

## 9. Nested sampling and checkpoint

```sh
sh tests/run_ns_test.sh <bin> <workdir> data/ramaprob.data tests/data/output.pdb
sh tests/run_checkpoint_test.sh <bin> <workdir> data/ramaprob.data tests/data/output.pdb
```

**Cross-binary checkpoint format check** (not part of the stock test
scripts): copied every `ckpt_*` file pristine's run produced (not just
`ckpt_0` — the restart flag `-R 2` reads `ckpt_2`, named by iteration, not
by write order) into a fresh directory alongside the same
`snapshots.pdb`/`ramaprob.data` fixtures, then ran current HEAD's `adcp`
with `-R 2` pointed at them. `rc=0` and correct resumption at iteration 3/4
confirms format compatibility independent of the content difference the NS
fix causes.

## 10. Tools/ correctness sweep

```sh
# byte-for-byte
<pristine bin>/cm    < tests/data/peptide12.pdb > out_p; <current>/cm    < tests/data/peptide12.pdb > out_c; cmp out_p out_c
timeout 5 <pristine>/ramachandran < tests/data/peptide12.pdb   # hangs, rc=124
<current>/rama                    < tests/data/peptide12.pdb   # completes

# bfactor: exact CMake smoke-test invocation
<bin>/bfactor < tests/data/bfactor_multimodel.pdb

# cdlearn: baseline and the long-token regression fixture
<bin>/cdlearn
<bin>/cdlearn -L tests/data/cdlearn_longtoken.txt -l V
```

**ASan confirmation for `bfactor`**, since the plain build didn't visibly
crash on the multi-model fixture (silent corruption, not a crash — the
interesting case):

```sh
gcc -std=c99 -O0 -g -fcommon -fsanitize=address,undefined \
  bfactor.c peptide.c aadict.c params.c error.c vector.c rotation.c \
  canonicalAA.c vdw.c -lm -o bfactor_asan
./bfactor_asan < tests/data/bfactor_multimodel.pdb
# AddressSanitizer: heap-buffer-overflow, vector.c:37 subtract, reading 8
# bytes past a 3016-byte allocation
```

**MemorySanitizer on the `adcp` binary itself**, following up on the initial
finding above. Built both binaries with clang; pristine needs `-fcommon`
(same reason as the gcc build) plus relaxing two now-hard-error diagnostics
clang enables by default for old C:

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
  -DCMAKE_EXE_LINKER_FLAGS="-fsanitize=memory" .
cmake --build build-msan-current --target adcp
```

Neither MSan build survives a plain run — both die on the first report.
**The fix is `MSAN_OPTIONS=halt_on_error=0:abort_on_error=0:exitcode=0`,
not `halt_on_error=0` alone**: MSan's own `-help=1` dump
(`MSAN_OPTIONS=help=1 <bin> -h`) says continue-after-report dedup is
"(asan only)", and empirically `halt_on_error=0` by itself still
terminates the process — `abort_on_error=0` is the flag that actually lets
execution continue past a report.

With that combination, current HEAD's fold reproduction runs to completion
past one report (`biasmap_initialise`, `energy.cpp:222`); the docking
reproduction (seed 10, redocking option set, 25,000 steps) surfaces three
distinct reports before terminating partway through (MSan's
continue-past-error support is best-effort, not guaranteed to reach the
end of a long-running program): `biasmap_initialise`,
`mark_aa_from_file` (`peptide.cpp:2002`), and — the one that matters —
`crankshaft` at `metropolis.cpp:727`.

**Ruling out the first two as tooling noise, not application bugs**: MSan's
own diagnostic under the `fopen`-adjacent reports reads
`Uninitialized bytes in fopen at offset 4 inside [...]` — a message from
MSan's own `fopen` interceptor about its internal state, not about
anything ADCP touches. This matches this repo's own `.github/workflows`
MSan job notes, which independently documents a confirmed false positive
from plain glibc's `fopen()` on this toolchain (found via `biasmap_initialise`
in that CI job too). Attempts to suppress it directly via
`-fsanitize-ignorelist=<file>` were tried (`fun:name`, `fun:*name*`,
`fun:*name*=uninstrumented`, `[memory]`-scoped sections, `src:<absolute
path>`) and none suppressed the specific reports — verified the ignorelist
mechanism itself works via a minimal reproducer and a blanket `fun:*`
pattern (both of which correctly suppressed everything), so the failure is
in pattern-matching specifics not worth chasing further once the
`Uninitialized bytes in fopen` message already gives a more direct answer.

**The `crankshaft` report is not tooling noise** — no libc/libstdc++ call
anywhere near the flagged expression, in either the C or the C++ version:
`chain->aa[i].chi1 != DBL_MAX`, guarded by `chain->aa[i].id != 'G' &&
chain->aa[i].id != 'A'` (glycine/alanine have no `chi1`, so this is meant
to skip them). Pristine's report is at the identical logical statement
(`metropolis.c:732`), reached through the identical call chain
(`crankshaft → move → simulate → main`), confirming the bug is
byte-identical across the whole C99 → C++17 migration — present before it,
untouched by it.

## 11. Timing methodology

- `/usr/bin/time -v` for anything single-shot, wall/user/sys/RSS all read
  from its output.
- Shell `time` (`{ time cmd; } > log 2>&1`) for the parallel 16-way redock
  runs, since `/usr/bin/time` only accounts the direct child, not the
  backgrounded `&` fan-out `run_redock_validation.sh` uses internally.
- `date +%s.%N` deltas inside `dock16.sh` for per-seed timing within a
  sequential sweep.
- No repeated-trial averaging was done beyond what's noted inline (fold and
  the 16-seed 10k sweep were re-run to check determinism, not to average
  timing) — single-machine, single-run wall-clock numbers throughout,
  explicitly not claimed as statistically rigorous benchmarks. Where a
  delta is within ~1%, 2.md calls it noise rather than a real effect.

## 12. Compile-time and LOC

```sh
wc -l /tmp/adcp-pristine/*.c /tmp/adcp-pristine/*.h
wc -l src/*.cpp include/*.h tools/*.cpp

rm -rf build_timing && /usr/bin/time -v gcc -std=c99 -O3 -DNDEBUG -fcommon \
  <15 pristine sources> -lm -o build_timing/adcp

rm -rf /tmp/adcp-current-build-timing
cmake -B /tmp/adcp-current-build-timing -DADCP_TOOLS=OFF \
  -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTING=OFF .
/usr/bin/time -v cmake --build /tmp/adcp-current-build-timing --target adcp -j1
```

Both single-threaded, so the numbers are compiler-bound rather than
parallelism-bound; the CMake figure includes CMake's own generator/build
overhead on top of the `g++` invocations, which the raw `gcc` figure does
not have an equivalent of — noted as a caveat in 2.md rather than corrected
for, since it reflects what a from-scratch build actually costs a user
either way.

## Cleanup

```sh
git worktree remove --force /tmp/adcp-pristine
git worktree remove --force /tmp/adcp-a1b6718
git worktree remove --force /tmp/adcp-bisect-wt   # if still registered
git worktree prune
```
