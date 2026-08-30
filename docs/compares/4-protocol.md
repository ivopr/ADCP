# Protocol for Comparison 4

Every command that produced a number in [4.md](4.md), in the order run.
Identical in structure to [3-protocol.md](3-protocol.md) — the point of this
comparison was again to see whether anything *changed* (this time: the
seed-10 UB investigation and `temp` fix), not to invent a new methodology.
Paths are as used during this session; substitute your own tmp locations.

## 0. Environment

```sh
date -u +"%Y-%m-%d %H:%M:%S UTC"   # 2026-08-30 18:10:00 UTC (session start)
gcc --version        # gcc (Ubuntu 13.3.0-6ubuntu2~24.04.1) 13.3.0
clang --version      # Ubuntu clang version 18.1.3 (1ubuntu1)
git rev-parse HEAD    # e28a5e3aea677ffc48f43551fe85de140ad62af3
```

Same machine as comparisons 2 and 3. Pristine is the same commit, `1c1a330`.
Current HEAD is `e28a5e3`, one commit ahead of comparison 3's `9c1e645` — the
seed-10 UB diagnosis and `temp` defense-in-depth fix.

The docking target cache (`build-dock/tests/target`) was gone again at
session start (third time this has happened across comparisons 2-4) —
`tests/fetch_target.sh` re-fetches it in ~5s, run it in parallel with the
pristine `git worktree add` to save a few seconds; both are independent.

## 1-2. Build pristine and current HEAD

Identical to [3-protocol.md](3-protocol.md) §1-2 — same 15-source pristine
`gcc -std=c99 -O3 -DNDEBUG -fcommon` invocation for `adcp_pristine`, same 9
tool-binary invocations, same CMake invocation for current HEAD
(`-DADCP_TOOLS=ON -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTING=OFF`). All 10
pristine binaries and current HEAD build clean, same warnings as before
(implicit-declaration, unused return values — pre-existing, not new).

## 3. CLI parity

Same as always: only the `Usage: <argv[0]>` path differs.

## 4. Folding determinism + timing

Same as [3-protocol.md](3-protocol.md) §4. Both binaries: `PASS: deterministic
over 2 runs, 76000 atoms, final energy 5.936278`. Cross-binary ATOM md5:
`6a438d0a673006235fccd2b1b7007ba3` on both — the same value recorded in
MIGRATION.md and every prior comparison.

## 5. Docking — 16-seed sweep at 10,000 steps

Same `dock16.sh` as [3-protocol.md](3-protocol.md) §5 (rewritten fresh this
session rather than reusing a prior session's `/tmp` copy, which doesn't
survive between sessions). 16/16 seeds identical `targetE` and ATOM-record
md5 between pristine and current HEAD, and every value equal to compares/2.md
and compares/3.md's own tables too.

## 6. Crash-rate mapping

Same `crashmap.sh` as before. Pristine: 12/16, 14/16, 14/16 crashes at
50k/150k/250k steps — **exact match**, survivor-for-survivor, to compares/2.md
and compares/3.md. Current HEAD: 0/16 at 250k.

## 7. Seed-10 divergence — now including the exact boundary from the UB fix

This is the one step with new content this session. In addition to the
standard binary-search points, probe the exact boundary MIGRATION.md's
"CONFIRMED" section pinned down (22,090 clean / 22,095 diverges on pristine):

```sh
for steps in 15000 20000 22090 22095 25000 30000 40000; do
  # run both binaries at $steps, seed 10, redocking options, diff targetE
done
```

Result:

| steps | pristine | current |
|---|---|---|
| 15,000-22,090 | -11.9092 | -11.9092 |
| 22,095 | -12.0108 | -11.9092 |
| 25,000+ | -15.3654 | -11.9092 |

Current HEAD stays flat at `-11.9092` through the *exact* step where
pristine first diverges — direct confirmation that the `temp` fix
(`e28a5e3`) changed nothing observable, anywhere, including at the single
most sensitive point this whole investigation ever found. This is the
predicted outcome from MIGRATION.md's own verification table, re-derived
independently rather than trusted.

## 8. Full-length redocking protocol

Same as [3-protocol.md](3-protocol.md) §8. Current HEAD: `PASS`, top-ranked
pose run 15, targetE -28.5508 (-16.91 kcal/mol), RMSD 0.66 Å over 28 atoms —
identical to every prior measurement. Pristine: `FAIL`, 15/16 runs
`Segmentation fault`, seed 10 the sole survivor — same pattern as before.

## 9. Nested sampling and checkpoint

Same as before. NS sequences match MIGRATION.md's recorded before/after
exactly on both sides. Checkpoint roundtrip `PASS` on both; cross-binary
restart (current HEAD resuming from pristine's own `ckpt_2`) confirmed
format-compatible again (`rc=0`, writes only `cross.pdb_3`/`cross.pdb_4`).

## 10. Tools/ correctness sweep

Same as [3-protocol.md](3-protocol.md) §10, same commands. Every result
matches: `cm` byte-identical, `ramachandran`/`rama` hangs-vs-completes,
`bfactor` multimodel hangs on pristine (`rc=124`) / rejects cleanly on
current (`rc=1`), `cdlearn` baseline and long-token findings unchanged.
ASan confirmation for `bfactor`: identical `heap-buffer-overflow`,
`vector.c:37` in `subtract`.

## 11. MemorySanitizer

Same recipes as [3-protocol.md](3-protocol.md) §11 (plain `-fsanitize=memory`,
no origin tracking needed for this comparison — that was last session's
diagnostic tool, not a standing part of this protocol). Fold, current HEAD:
1 known false positive (`biasmap_initialise`, `energy.cpp:222`), completes.
Docking, seed 10, 25,000 steps: pristine shows `crankshaft`/`metropolis.c:732`
(frozen, always will); current HEAD shows only the two known `fopen`-adjacent
false positives, no `crankshaft` report — same as compares/3.md.

## 11a. UBSan confirmation of the `swapChains` mechanism (new this session)

Not part of comparisons 2/3's protocol — added after MIGRATION.md's
"CONFIRMED" section made this a standing, re-checkable claim. Build both
binaries with `-fsanitize=address,undefined -fno-sanitize-recover=all`
(pristine: `gcc -std=c99 -O0 -g -fcommon`, same source list; current: CMake
with `-DCMAKE_CXX_FLAGS="-fsanitize=address,undefined -fno-sanitize-recover=all -O0 -g"`).
Run the seed-10/22,095-step reproduction with `ASAN_OPTIONS=halt_on_error=0`
(needed on pristine to get past the unrelated, already-documented
`params.c:1121` `%256[` overflow and actually reach the docking loop).

Pristine: `main.c:237:14: runtime error: index 10 out of bounds for type
'Chain *[*]'` — reproduced exactly, same file:line as MIGRATION.md's
"CONFIRMED" section. Current HEAD: 0 errors, `rc=0`, completes normally.

## 12. Compile-time and LOC

Same recipes as [3-protocol.md](3-protocol.md) §12. LOC: pristine 24,441
(unchanged, frozen), current 22,778 (up 11 lines from compares/3.md's 22,767
— the `temp` fix's explanatory comment). Compile time, single-threaded:
pristine 3.98s wall / 3.78s user; current 7.22s wall / 6.62s user — both
within noise of compares/3.md's 3.99s/7.07s.

Also grabbed: binary size (`ls -la` on both `adcp` binaries) — pristine
409,728 B (frozen), current 497,248 B, unchanged from compares/3.md (the
`temp` fix is a declaration-line change, no new code paths).

## 13. Final sanity check

```sh
ctest --test-dir build --output-on-failure
```
19/19, confirming the comparison work didn't perturb the working tree.

## Cleanup

```sh
git worktree remove --force <pristine-worktree-path>
git worktree prune
rm -rf <all scratch build/work directories under your tmp>
```

The docking target cache (`build-dock/tests/target`) is worth keeping this
time too — it has now needed re-fetching in every one of comparisons 2, 3,
and 4. If this keeps happening, consider whether something is deleting it
between sessions, or just budget the ~5s re-fetch as a standing cost of this
protocol's step 0.
