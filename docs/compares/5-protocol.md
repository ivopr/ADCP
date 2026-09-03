# Protocol for Comparison 5

Every command that produced a number in [5.md](5.md), in the order run.
Same structure as [4-protocol.md](4-protocol.md), with two new arms —
upstream's `cyclic` branch and the ADFRsuite 1.0 binary — and one new step
(§7a, the per-commit attribution bisect). Paths are as used during this
session; substitute your own tmp locations.

## 0. Environment

```sh
date -u +"%Y-%m-%d %H:%M:%S UTC"   # 2026-09-02 14:06:53 UTC (session start)
gcc --version        # gcc (Ubuntu 13.3.0-6ubuntu2~24.04.1) 13.3.0
clang --version      # Ubuntu clang version 18.1.3 (1ubuntu1)
nproc                # 16
git rev-parse HEAD   # d409acb6b117d6c8f3ecfaaed52dced1ea06e1dc
```

Third arm — the vendor binary:

```sh
file ~/ADFRsuite-1.0/bin/adcp_Linux-x86_64
# ELF 64-bit LSB executable, x86-64, dynamically linked, with debug_info,
# not stripped; mtime 2019-07-15, 788,549 B
```

Fourth arm — upstream's `cyclic` branch. It is not fetched by default; the
tip is 13 commits `upstream/master` never merged:

```sh
git fetch upstream cyclic:refs/remotes/upstream/cyclic
git log --oneline -1 upstream/cyclic          # 1d34035 fix HillClimb
git log --oneline upstream/master..upstream/cyclic | wc -l   # 13
```

The docking target cache (`build-dock/tests/target`) was gone again — fourth
comparison running in which this has happened. `sh tests/fetch_target.sh
build-dock/tests/target` re-fetches in ~5s; run it in parallel with the
pristine `git worktree add`.

The `noclobber` and shadowed-`diff` gotchas from
[3-protocol.md](3-protocol.md) §0 both bit again this session. Use `>|` for
every log redirect and `command diff` for every diff.

**Do not `pkill -f` on a pattern that could match your own shell** — one
attempt to kill a stale wait-loop this session killed the calling shell too
(exit 144). Match on the script name, or just let the loop time out.

## 1-2. Build pristine, current HEAD and `cyclic`

Identical to [3-protocol.md](3-protocol.md) §1-2: the same 15-source
`gcc -std=c99 -O3 -DNDEBUG -fcommon` invocation for `adcp_pristine` (409,728 B,
byte-size unchanged since comparison 2), and CMake with `-DADCP_TOOLS=ON
-DCMAKE_BUILD_TYPE=Release -DBUILD_TESTING=OFF` for current HEAD. The ADFR
arm needs no build.

`cyclic` builds with **pristine's exact command line, unmodified** — same 15
sources, same flags, same warnings, no Makefile edits:

```sh
git worktree add -f <tmp>/cyclic upstream/cyclic
cd <tmp>/cyclic && mkdir -p build_cyclic
gcc -std=c99 -O3 -DNDEBUG -fcommon \
  nested.c aadict.c energy.c main.c metropolis.c flex.c peptide.c probe.c \
  rotation.c vector.c params.c error.c checkpoint_io.c vdw.c canonicalAA.c \
  -lm -o build_cyclic/adcp_cyclic
```

Zero errors; 409,728 B, the same size as pristine's binary. Ignore
`metropolis_old.c` — it is present in the branch but not in the build.

The tools/ correctness sweep ([4-protocol.md](4-protocol.md) §10) was **not**
re-run this session: it has been identical across comparisons 2, 3 and 4, no
commit since touched a tool source, and ADFRsuite ships no tool binaries, so
there is no third-arm question to answer. `ctest` (§13) covers the tools that
have tests.

## 3. CLI parity

```sh
for b in "$PRISTINE" "$CURRENT" "$CYCLIC" "$ADFR"; do "$b" -h > h_<label>.txt 2>&1; done
command diff -u h_pristine.txt h_current.txt   # only Usage: argv[0]
command diff -u h_pristine.txt h_adfr.txt      # banner (2 lines vs 3) + argv[0]
command diff -u h_cyclic.txt   h_adfr.txt      # only Usage: argv[0]
```

`cyclic` and ADFR both exit **1** on `-h` and print
`PepTide 2.0, Copyright (c) 2004 - 2010 Alexei Podtelezhnikov`; pristine and
current exit 0 and print the two-line `ADCP 0.1` banner. That split — banner
*and* exit code — is the first of the two lineage signals. The entire
options block is identical on all four arms, so every flag this protocol
uses is portable everywhere.

## 4. Folding determinism

```sh
sh tests/run_fold_test.sh <bin> <workdir> data/ramaprob.data
md5sum <workdir>/atoms1.txt
```

All four PASS their internal 2-run determinism check, in two pairs:

| arm | ATOM md5 | final energy |
|---|---|---|
| pristine, current HEAD | `6a438d0a673006235fccd2b1b7007ba3` | 5.936278 |
| `cyclic`, ADFRsuite | `abf915bc58b0eaedf4e4d3a1353f723b` | 9.274271 |

This is the second and stronger lineage signal: a byte-identical fold
trajectory between the public `cyclic` tip and the shipped binary. It is
also the cheapest check in this whole protocol (~2 min) — run it first if
you only want to re-confirm the lineage claim.

Budget note: the ADFR fold run is slower than either local `master` build
and does not finish inside a 2-minute tool timeout — give it its own call
with a longer budget rather than looping arms in one command.

## 5. Docking — 16-seed sweep at 10,000 steps

Same `dock16.sh` as [3-protocol.md](3-protocol.md) §5, verbatim, run for all
four arms (they are independent — run them concurrently, the whole sweep
finishes in about 15s wall). All 64 runs exit 0.

Pristine: 16/16 `targetE` **and** ATOM md5 equal to compares/2-4. Current
HEAD: all 16 moved. `cyclic`: 10/16 md5-identical to current HEAD. ADFR:
distinct from every other arm.

Count the agreement mechanically rather than by eye — `targetE` can match
while the pose differs, and that distinction carries information:

```sh
join <(tail -n+2 dock16_current.txt) <(tail -n+2 dock16_cyclic.txt) |
  awk '{if($3==$7) m++; else if($2==$6) e++; else d++}
       END{print "md5-identical="m+0, "energy-equal-pose-diff="e+0, "different="d+0}'
# md5-identical=10 energy-equal-pose-diff=2 different=4
```

Field offsets matter here: `join` emits key + 4 fields from each file, so
the md5 columns are `$3`/`$7`, not `$3`/`$6`. Getting this wrong reports
16/16 different and quietly inverts the finding.

## 6. Crash-rate mapping

Same `crashmap.sh` as before. Run pristine, `cyclic` and ADFR at 50,000 /
150,000 / 250,000 and current at 250,000 — ten configurations, all
independent, all launchable concurrently.

Pristine 12/16, 14/16, 14/16 crashes — exact match to compares/2-4. Current,
`cyclic` and ADFR: 0/16 at every budget. Pristine is the only arm that
crashes at all.

## 7. Seed-10 boundary probe

Same as [4-protocol.md](4-protocol.md) §7, extended to four arms:

```sh
for steps in 15000 20000 22090 22095 25000 30000 40000; do
  # run each of the four binaries at $steps, seed 10, redocking options
done
```

Pristine's 22,090-clean / 22,095-diverges boundary reproduces exactly.
Current HEAD and `cyclic` are both flat across it **at the same value**
(-13.0192); ADFR is flat at its own (-11.5366). The divergence is
`master`-only.

## 7a. Per-commit attribution of the current-HEAD move (new this session)

The one genuinely new step. compares/4's headline was that HEAD matched
pristine; this session's sweep showed it no longer does, so the movement has
to be attributed to a commit rather than assumed.

```sh
for c in e28a5e3 c6dcf24 7ea27e1 58337c5 797eb72 7e38b98; do
  git worktree add -f wt_$c $c
  cmake -B b_$c -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTING=OFF wt_$c && cmake --build b_$c -j4
done
# then: seed 1 and seed 10, 10,000 steps, redocking options, on each build
```

Six builds, twelve runs, under two minutes total. Result: the step is
entirely between `c6dcf24` and `7ea27e1`. See 5.md's table.

This is worth keeping as a standing step whenever the 16-seed sweep moves:
it converts "something changed" into "this commit changed it" for the cost
of a handful of Release builds.

## 8. Full-length redocking protocol

```sh
{ time sh tests/run_redock_validation.sh <bin> build-dock/tests/target/stage \
    <workdir> 16 2500000 2.5 ; } >| log 2>&1
```

Run once per arm, **serially** — each invocation already fans 16 seeds out
in parallel and saturates the machine.

- current HEAD: `PASS`, 0.70 Å, run 4, targetE -30.3416, wall 1:23.71
- ADFR: `PASS`, 0.77 Å, run 6, targetE -28.2783, wall 1:37.33
- `cyclic`: `PASS`, 0.88 Å, run 1, targetE -29.1304, wall 1:28.16
- pristine: `FAIL`, 15/16 `Segmentation fault`, seed 10 sole survivor, wall 1:07.50

## 9. Nested sampling and checkpoint

```sh
sh tests/run_ns_test.sh <bin> <workdir> data/ramaprob.data tests/data/output.pdb
sh tests/run_checkpoint_test.sh <bin> <workdir> data/ramaprob.data tests/data/output.pdb
```

Pristine, `cyclic` and ADFR all give the identical NS sequence
(`40.058829 31.225769 10.764044 7.463047 6.393763`); current gives
`35fb3fb`'s post-fix sequence. All four checkpoint roundtrips PASS. Note NS
does **not** split along lineage lines the way fold does — `cyclic`'s
changes never touched it.

Cross-binary restart, run for all three foreign arms (options must match
`grep '^COMMON=' tests/run_checkpoint_test.sh` verbatim):

```sh
cp <foreign-workdir>/ckpt_* <foreign-workdir>/snapshots.pdb data/ramaprob.data <workdir>/
( cd <workdir> && <current>/src/adcp -n -f snapshots.pdb -r 2x10 -s 12345 \
    -p Bias=NULL -C 2,ckpt -R 2 -o cross.pdb )
```

Current HEAD resumes from pristine's checkpoints, from `cyclic`'s, and from
the ADFR binary's — all `rc=0`, all writing only
`cross.pdb_3`/`cross.pdb_4`. The checkpoint format is stable across all four
arms and both lineages.

## 10. Sanitizers

Not re-run this session. MSan and the UBSan `swapChains` confirmation
([4-protocol.md](4-protocol.md) §11, §11a) both target pristine's frozen
behaviour and current HEAD's fix of it; neither has changed, the §7 boundary
probe still shows the fix holding, and the vendor binary has no source to
instrument. Re-run them when the `swapChains` claim is in doubt, not as
routine.

## 11. Size, LOC, final sanity check

```sh
stat -c "%n %s" <pristine-bin> <current-bin> <cyclic-bin> \
  ~/ADFRsuite-1.0/bin/adcp_Linux-x86_64
cat src/*.cpp include/*.h tools/*.cpp | wc -l   # current: 22,906
cat <pristine>/*.c <pristine>/*.h | wc -l       # pristine: 24,441 (frozen)
cat <cyclic>/*.c <cyclic>/*.h | wc -l           # cyclic: 24,447
ctest --test-dir build                          # 22/22 pass
```

Note `ctest` is at 22 tests now, up from compares/4's 19 — the 49-target
validation set (`82b4c12`) and the smoke-test fixes since. LOC counts must
exclude `CMakeLists.txt` to stay comparable with earlier comparisons'
numbers.

## Cleanup

```sh
git worktree remove --force <pristine-worktree> <cyclic-worktree> \
  <each wt_* from §7a>
git worktree prune
rm -rf <all scratch build/work directories under your tmp>
```

Keep `refs/remotes/upstream/cyclic` — it costs nothing and saves a fetch
next time. Keep `build-dock/tests/target` if you can — it has now needed re-fetching in
comparisons 2, 3, 4 and 5, so budget the ~5s either way.
