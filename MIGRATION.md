# ADCP: C99 → Modern C++17 Migration Plan

## Context

ADCP (AutoDock CrankPep) is a 23.7k-line C99 peptide-docking Monte Carlo engine, distributed as a CLI binary (invoked as a subprocess by `runADCP.py`/AGFR, not linked as a library by any known external consumer). The codebase does heavy manual memory management (100 malloc, 76 realloc, 159 free calls) and has a proven production bug in this exact area: commit `7e2192d` fixed a stack-buffer-overflow VLA (`Chain* swapChains[swapLength]` sized one element short) in `src/main.c` that crashed every production docking run for 200k+ step searches. Be precise about what would have caught it: `std::vector` accessed through `operator[]` would not have, and this plan mandates `operator[]` in hot paths. An ASan build would have, which is why the sanitizer CI leg below is the load-bearing mitigation for this bug class, not the container conversions. Commit `5e94b0b` then added the only tests that exercise the docking energy path at all (`dock_3q47_smoke`, `val_3q47_redock`); everything else runs the folding path where the receptor energy term is literally `0.000000`.

Goal: migrate the whole tree to idiomatic modern C++17, incrementally, file-by-file, keeping every step compiling and green on the gcc+clang CI matrix, without ever breaking the docking path this codebase already got burned on once. (**Read the CI hardening section before trusting that phrase.** Throughout Phase 1 the "gcc+clang matrix" compiled zero C++ with clang — the workflow set `CC` and never `CXX`. Every "validated on gcc+clang" claim in this document and in the Phase 1 commit messages was false until the fix described below.) No exceptions as the error-handling idiom (matches current `stop()`/`exit()` style, avoids MPI/OpenMP unwinding hazards) — RAII is fine, `throw`/`catch` is not (one unavoidable exception: `std::vector`/`std::string` throw `bad_alloc` on OOM, which is accepted as behaviorally equivalent to the current fatal-OOM behavior).

## Standard & build setup (do once, first)

- Target **C++17**. Gives `if constexpr` (for zero-overhead numeric dispatch replacing function-pointer tables), `std::string_view`, structured bindings. C++20 buys nothing here (no generic-library surface needing concepts/ranges) — don't reach for it.
- Root [CMakeLists.txt](CMakeLists.txt): change `project(ADCP VERSION 0.1 LANGUAGES C)` → `LANGUAGES C CXX`, add `CMAKE_CXX_STANDARD 17` / `CMAKE_CXX_STANDARD_REQUIRED ON` / `CMAKE_CXX_EXTENSIONS OFF` alongside the existing C99 settings. Keep both C and C++ language settings until the last `.c` file is gone, then drop the C ones in a final cleanup commit.
- `-fno-common`/`ADCP_LEGACY_COMMON` was only ever applied to `CMAKE_C_FLAGS`, so it never reached a single C++ TU; once `src/` went all-C++ the option was a silent no-op. Phase 1 step 7 mirrors both branches into `CMAKE_CXX_FLAGS`. Note this is hygiene, not a behavior fix: g++ emits no COMMON symbols for C++ at all (verified with `nm`: zero `C`-type symbols either way), so the flag is inert for C++ TUs — it is set so the escape hatch stays honest, and it matters for `tools/*.c` until Phase 1 step 8 lands.
- Per-file conversion mechanics (repeat for every file, every phase): `git mv src/foo.c src/foo.cpp` (preserves blame — files carry LGPL attribution headers), update the one line in [src/CMakeLists.txt](src/CMakeLists.txt) or `tools/CMakeLists.txt`, build with gcc+clang, run tests, commit. CMake mixes `.c`/`.cpp` sources in one target natively — no "big bang" moment is structurally required, and the existing CI matrix job keeps validating the whole tree at every commit without workflow changes.
- ~~Recommended addition during Phase 1: MPI/OpenMP legs and an ASan+UBSan config.~~ **Done, but only after Phase 1 had already finished** — see "CI hardening" below. Recommending it and not building it is what let the `energy.cpp` regression through; the sanitizer leg found two real bugs within minutes of first existing.

## Phase 1 — whole-tree rename to `.cpp`, compile clean, no behavior change

Rename every `.c` in `src/` and `tools/` to `.cpp`, bottom-up by include-dependency order, fixing only what's needed to compile as standard C++ (mainly explicit casts on `malloc`/`void*`, which C++ requires and C doesn't):

1. `error.c`, `vector.c` — no local deps, prove the toolchain
2. `rotation.c`, `random16.c`
3. `aadict.c`, `canonicalAA.c`
4. `params.c` — heaviest allocator (57 calls), surfaces most cast issues early while changes are still "just compile," not "redesign"
5. `peptide.c` — owns the core structs (`Chain`, `AA`, `Biasmap`, `FLEX_data`)
6. `vdw.c`, `energy.c` (largest file, 3123 lines), `probe.c`, `metropolis.c`, `flex.c`, `checkpoint_io.c`, `nested.c`
7. `main.c` — two changes beyond casts. `double swapEnergy[swapLength + 1]; Chain* swapChains[swapLength + 1];` at [src/main.cpp:214-216](src/main.cpp#L214-L216) are VLAs, not valid standard C++ (gcc/clang accept them only as a non-portable extension, and the build sets `CMAKE_CXX_EXTENSIONS OFF`). **Fix: mark `swapLength` `const`**, not `std::vector`. `swapLength` is `10` and never reassigned, so `const` makes `swapLength + 1` a constant expression and both arrays become plain fixed-size arrays — the VLA is gone with a one-word diff, no heap allocation, and none of the ~40 indexing sites change. `std::vector` was the original prescription here on the grounds that bounded containers would have caught `7e2192d`'s off-by-one; that reasoning does not hold, because this plan also (correctly) mandates `operator[]` over `.at()` in hot paths, and `operator[]` is unchecked. ASan is what catches that bug class — see the Phase 1 CI recommendation above. Run `dock_3q47_smoke` (the regression guard added in `5e94b0b`) immediately after this file.
8. `tools/*.c`, simplest first (`pauling.c`, `bfactor.c`, `ramachandran.c`, `cm.c`, `dssp2cm.c`, `mergie.c`, `oops.c`, `statistics.c`), `cdlearn.c` last (only OpenMP consumer, none of `tools/` is on the docking test path per [tests/CMakeLists.txt](tests/CMakeLists.txt)).

## Phase 2 — idiomatic conversion, same order, one module (or tight group) at a time

Apply per category:

- **malloc/free arrays → `std::vector`**: `Chain`'s array members in `peptide.h`/`peptide.c`, `params.c`'s option-list arrays, `flex.c`'s index lists, `checkpoint_io.c`'s buffers. Use `operator[]` not `.at()` in hot paths (matches existing unchecked-access behavior).
- **Fixed C string buffers → `std::string`**: note that `main.cpp`'s `char swapname[12]` and `FILE *swapFile` ([src/main.cpp:209-213](src/main.cpp#L209-L213)) are **dead code** — every consumer is commented out, `swapname` is only ever written by one `sprintf` and `swapFile` never leaves `NULL`. Delete both rather than converting them.
- **Function-pointer dispatch → templates, only in hot paths**: `vdw.c`'s `vdw_fn` sites (lines 124, 219, 448, 751) and `probe.c:155` fire every Monte Carlo step — replace with a C++17 non-type template parameter (`template<auto Fn> ...`) or a flattened `switch` calling each candidate directly, whichever is the smaller diff once you read the call sites; either gives zero-overhead dispatch. `energy.c`'s `sumf`/`sumf_diag`/`dfdx` (~line 2720+) are cold diagnostic/finite-derivative code, not the hot loop — leave as raw function pointers unless converting is free while touching surrounding code anyway.
- **`FILE*` → RAII wrapper, not iostream**: don't rewrite the `scanf`/`fprintf`-based PDB-format parsing (high risk, no real benefit) — just wrap the handle lifetime (open/close pairing) in a small RAII class so early-return paths can't leak. Apply to `checkpoint_io.c` and `flex.c`. (`main.c`'s `swapFile` needs no wrapper — it is dead, see the `std::string` item above.)
- **3 `goto` sites** — none are cleanup patterns (verified: no `free`/`fclose` at either label):
  - `energy.c:291`, `energy.c:2876` — conditional-skip patterns, restructure as `if`/loop-flag.
  - `tools/ramachandran.c:165` — forward jump to reuse an already-read record; extract the loop body into a local lambda called once before the loop and once inside it. `smoke_rama_pdb` guards a related historical bug in this exact file — diff its output before/after.
- **MPI/OpenMP**: no `extern "C"` needed, both headers are C++-safe as-is. Verify empirically via the new CI legs, no code change expected.

Ordering within Phase 2, by risk × blast radius:

1. **Done.** Leaves: `error.c`, `vector.c`, `rotation.c`, `random16.c`, `aadict.c`, `canonicalAA.c`
2. **Done.** `params.c` — do before anything downstream depends on its idiomatic shape
3. `peptide.c` — owns the shared structs; converting `Chain`'s arrays touches every file that reads `chain->aa[i]`, so do it once, early, before those callers are themselves converted
4. `main.c` idiomatic pass (string/RAII file wrapper — the VLA fix already landed in Phase 1)
5. `vdw.c` (dispatch templating)
6. `energy.c`
7. `probe.c`, `metropolis.c`
8. `flex.c`, `checkpoint_io.c`
9. `nested.c`
10. `tools/*`
11. **Done, ahead of schedule.** Cleanup: drop `LANGUAGES C` / `CMAKE_C_STANDARD*` once `find src tools -name '*.c'` is empty — `CMakeLists.txt` has been `LANGUAGES CXX` only since Phase 1 finished, no C sources ever remained afterward.

## Phase 2 progress — steps 1 and 2 (params.cpp) done

**Step 1 (leaves).** `error.cpp`, `vector.cpp`, `rotation.cpp`, `random16.cpp`,
`canonicalAA.cpp`: reviewed against all four idiomatic categories, nothing to convert —
Phase 1's cast-only pass already left them idiomatic. `aadict.cpp`: the one real site
(a local `sprintf`-into-fixed-buffer error message) converted to `std::string`.

**Step 2 (params.cpp).** Every `malloc`/`realloc`/`free`'d field in `model_params`,
`FLEX_params` and `simulation_params` that this plan's Phase 2 section named is converted:

- `simulation_params`'s `energy_gradient`, `energy_probe_1_this`, `energy_probe_1_last`
  (36-element `double*`) and `energy_probe_1_calc` (`int*`) → `std::vector`.
- `FLEX_params`'s `output_path`, `outputpdb_filename`, `flex_cmd`, `flex_dir` and
  `model_params`'s `contact_map_file`, `fixed_aalist_file`,
  `external_constrained_aalist_file`(`2`) → `std::string`.
- `FLEX_params`'s `filenames_to_read_in` (`char**`) → `std::vector<std::string>`.
- `model_params`'s `vdw_gamma_gamma_cutoff`, `vdw_gamma_nongamma_cutoff` (702-element
  `double*`) → `std::vector<double>`.
- `model_params`'s `sidechain_properties` (31-element `sidechain_properties_*`) →
  `std::vector<sidechain_properties_>` — the widest blast radius of the pass: it's passed
  by raw pointer into helper functions across 13 files (`aadict.cpp`'s
  `sidechain_vdw_radius`, `sidechain_dihedral`, `hbond_donor`/`acceptor`, `charge`, etc.).
  Those helpers keep their `sidechain_properties_ *` signatures unchanged; every call site
  now passes `.data()` instead of the field directly, safe because `std::vector` guarantees
  contiguous storage.

Left out of this pass (not named in the plan's representative sites, still raw
`char*`/`malloc`): `simulation_params`'s generic strings (`seq`, `sequence`, `infile_name`,
`outfile_name`, `prm`, `checkpoint_filename`) and `MC_lookup_table`/`MC_lookup_table_n`
(dynamically-sized, `NULL`-checked). `params.cpp`'s malloc/realloc/free count dropped from
20/10/26 to 4/0/9 — the remainder is exactly these left-out fields plus the generic
`copy_string` helper they still use.

**A real bug found along the way, not a migration regression in the usual sense:** two
call sites — `vdw.cpp`'s `vdw_cutoff_distances_calculate` and `cdlearn.cpp`'s
per-protein `simulation_params` array setup — `malloc`/`realloc`'d raw `simulation_params`
memory and then called `sim_params_copy` on it (which does `*to = *from`). That's fine for
a plain-old-data struct. It is undefined behaviour the moment the struct owns a
`std::vector`/`std::string` member (`operator=` runs on an unconstructed object), and it
segfaulted `func_adcp_fold_short`/`determinism` immediately once `simulation_params` grew
its first `std::vector` field. Fixed by switching both sites to `new`/`delete`.
`model_params`'s one analogous site (`energy.cpp`'s `energy_probe_1` scratch copy) was
fixed pre-emptively before `model_params` gained its own vector members. **Lesson for
`peptide.c` (Phase 2 step 3, next up):** before converting any field in a struct, grep for
raw `malloc(sizeof(TheStruct))` / `realloc(..., n * sizeof(TheStruct))` on that struct type
across the whole tree and fix those sites to `new`/`delete` first — this bug class is cheap
to have missed and expensive to debug once hit.

Validated per group: full `ctest` suite (13/13) after every conversion; the determinism
test hammered 12× after every group (all clean); the `docking|validation` label gate run
after the `params.cpp` groups specifically, since `realloc`→container conversions are the
exact pattern flagged in the risk register below — one run got a live network fetch of the
3Q47 target and passed `dock_fetch_target`/`dock_3q47_smoke`/`val_3q47_redock` in full.

## Validation gates

- **Every file conversion**: build gcc+clang, run default `ctest` (smoke + functional, ~13 tests, seconds).
- **Any change touching `energy.cpp`, `vdw.cpp`, `metropolis.cpp`, `probe.cpp`, or `main.cpp`'s swap logic**: also configure `-DADCP_DOCKING_TESTS=ON` and run both `docking` and `validation` labels locally. CI now runs both (see below), so this is a fast local pre-check rather than the only gate. `val_3q47_redock` takes **~2 minutes** — measured at 117s, 120s and 140s. An earlier draft of this document claimed ~30 min and used that to justify leaving it out of CI; that number was wrong.
- **`func_adcp_fold_determinism`**: run after every `energy.c`/`vdw.c` change without exception. It's the only automated determinism signal, and template-vs-function-pointer codegen changes can alter floating-point accumulation order in Monte Carlo sums (IEEE 754 non-associativity) — a determinism regression here is a blocking failure to investigate, not acceptable drift.

## FIXED REGRESSION — `energy.cpp`, introduced by `5660a33` (Phase 1 step 6)

**Status: fixed (trigger removed). Mechanism still unexplained — read this before touching
`scoreSideChain` / `scoreSideChainNoClash` again.**

The fix: `tc` is a fixed-size stack array sized at the largest rotamer set declared in
`canonicalAA.h` (81 × 11 × 3 floats, ~10 KB), indexed through the same `View3` wrapper.
That is standard C++ — no VLA, so no `gnu++17` needed — and keeps the buffer on the stack.
Measured 0/15 hangs.

`func_adcp_fold_determinism` hangs intermittently — roughly 1 run in 3, same seed, same
command. Bisected by running the identical fold 12 times per build:

| Build | Hangs |
|---|---|
| Pre-migration C (`5e94b0b`) | 0/12 |
| Before step 6 (`adf422e`) | 0/12 |
| Step 6 (`5660a33`) | 3/12 |
| `90c63a8` (current) | 4/12 |
| Current, with only `energy` reverted to the C file | 0/12 |
| Current, with only `tc` reverted from `std::vector` to a VLA | 0/12 |

So the trigger is precisely one change: in `scoreSideChain` and `scoreSideChainNoClash`,
`float tc[nbRot][nbAtoms][3]` became `std::vector<float> tc_storage` + a `View3` wrapper.
Symptom: `currTargetEnergy` picks up a denormal (~1e-308, a different value each time),
which makes the retry loop at [src/main.cpp:306](src/main.cpp#L306) spin forever, since
every `swapEnergy[i]` compares `>=` against it.

Mechanism: **not established.** Seven hypotheses were tested and every one was ruled out
by measurement:

| Hypothesis | Result |
|---|---|
| `a->SCRot` indexes `tc` out of bounds (only assigned inside `if (score < bestScore)`, and ADCP transmutates residues, so a stale index from an 81-rotamer set could address a 9-rotamer one) | Refuted — an instrumented build logging every `SCRot` outside `[0, nbRot)` recorded **zero** violations on runs that hung |
| A memory error ASan/UBSan would catch | No report in 6 runs, and the hang never reproduced under ASan at all |
| Stack frame size changed by removing the VLA | Refuted — re-adding a `volatile` pad of identical size still hung 4/12 |
| `chain->erg` is `realloc`'d without zeroing, so `totenergy` sums heap residue | Refuted — `memset`ing it still hung 3/12 |
| `View3`'s indexing differs from the C99 VLA parameter | Refuted by review — `p[i*nbAtoms*3 + j*3 + k]` in both, and a C99 VLA parameter also uses the *runtime* `nbAtoms` |
| Strict-aliasing violation exploited at `-O3` | Refuted — `-fno-strict-aliasing` still hung 2/12 |
| Optimization level | Refuted — plain `-O2` still hung 3/12 (so ASan's clean runs came from the sanitizer, not from `-O2`) |

Note the `tc` change made memory *more* initialized, not less — `std::vector<float>`
zero-initializes where the VLA did not — so a naive "the vector reads uninitialized data"
story does not hold. What is established is narrow and empirical: **heap `tc` hangs, stack
`tc` does not**, across every other variable tested. Pinning the mechanism needs MSan or
valgrind; neither is installed on this machine.

**Consequence for the rest of the migration: validate by hammering, not by reasoning.** Run
`func_adcp_fold_determinism`'s command 12–15 times, not once. A single green run proves
nothing here — this regression passed 13/13 plus docking and validation when it was
committed.

Two lessons worth keeping:

- **The container conversion did not catch a latent bug — it converted a silent one into a
  hang.** This is the opposite of the premise in the Context section above, and it is the
  strongest argument yet for building the ASan/UBSan leg *before* converting anything else.
- `5660a33`'s commit message claims "compile-only changes: explicit casts" and "No behavior
  change". That is false: it also rewrote two VLA signatures behind a hand-written `View3`
  template and restructured a `goto` that this plan had scheduled for Phase 2. It passed
  13/13 plus docking and validation at the time, because the failure is intermittent. **A
  green test run is not evidence that a diff is cast-only — read the diff.**

## CI hardening (done after Phase 1, before starting Phase 2)

Phase 2 is made entirely of `malloc`→`std::vector` conversions — the same class of change
that produced the unexplained `energy.cpp` regression. The CI as it stood would not have
caught that, so it was rebuilt first. [.github/workflows/ci.yml](.github/workflows/ci.yml)
now has five jobs:

| Job | What it adds |
|---|---|
| `build-and-test` (gcc, clang) | **Sets `CXX`, not `CC`.** The project is `LANGUAGES CXX`, so `CC` selected nothing; both legs had been building with the default g++, and clang coverage had decayed to zero as the tree converted. Also runs the determinism test with `--repeat until-fail:12`. |
| `sanitizers` | ASan+UBSan with `-fno-sanitize-recover=all` (without it UBSan prints and continues, and the job goes green with real findings). |
| `openmp` | `-DADCP_OPENMP=ON`; `cdlearn` is the only consumer and nothing else exercises it. |
| `mpi` | `-DADCP_MPI=ON`, build only — no test covers the `PARALLEL` path. `adcp_mpi` was never compiled once during the entire migration. |
| `docking` | Now `-L "docking|validation"`, not just `-L docking`. |

The `No COMMON symbols` step is kept but is now **vacuous**: g++ emits no COMMON symbols for
C++ regardless of `-fno-common` (confirmed with `nm`, zero either way). It guarded C
tentative definitions and there is no C left. Delete it if you would rather not carry a
check that cannot fail.

Why the repeat matters: **one green run proves nothing when the failure is intermittent.**
The `energy.cpp` regression hung ~1 run in 3 and still passed 13/13 plus docking and
validation the day it was committed.

## Pre-existing bugs found by ASan (NOT migration regressions)

`sscanf(prm, "Bias=%256[^,]...", contact_map_file)` at
[src/params.cpp:1121](src/params.cpp#L1121) writes into `char contact_map_file[256]`.
`%256[` reads up to 256 characters **and then appends a NUL** — 257 bytes into a
256-byte buffer. ASan reports it as a stack-buffer-overflow on the first run.

Verified present, at the same line, in the pre-migration C build (`5e94b0b`), so it is not
caused by the migration. The same off-by-one appears at **12 sites** across `params.cpp`
and `main.cpp` (`%256[^,]` and `%256s` into 256-byte buffers). Fix is `%255[^,]` / `%255s`
at every site. It needs a ≥256-character path to trigger, so it is latent — and it is
almost certainly *not* the cause of the hang above.

**Fixed** in its own commit, separate from the migration steps: all 12 sites now use
`%255[^,]` / `%255s`.

A second one, found the first time the sanitizer CI leg ran: `dssp2cm`'s `parse_dssp_body`
wrote its string terminator out of bounds. `seq`/`ss` were `calloc`'d with `n_res` bytes but
are indexed by the DSSP residue number, which reaches `n_res-1`, with the terminator written
one past that — so the overflow happens on valid input, not only on the empty smoke fixture.
`k` was also uninitialized, so on empty input the terminator index was indeterminate. Both
present in the pre-migration C. Fixed in `8392e34`; output byte-identical.

Two real bugs within minutes of the sanitizer leg first existing is the argument for having
built it in Phase 1, as this document originally recommended.

## Risk register

| Risk | Where | Mitigation |
|---|---|---|
| Converting a stack VLA to a heap container relocates any out-of-bounds access, turning a benign latent read into garbage data | `energy.cpp`'s `tc` (already bit us, see above); any remaining VLA→container change | Build the ASan/UBSan leg first; hammer `func_adcp_fold_determinism` 12× rather than once, since these failures are intermittent |
| `std::vector` reallocation growth differs from `realloc`, could unmask a latent buffer over/under-read previously hidden in slack space | any `realloc`→`push_back` conversion (`params.c`, `flex.c`) | ASan/UBSan CI leg (Phase 1); full functional+docking+validation run after each |
| Template/inlining changes float accumulation order | `vdw.c`, `energy.c` dispatch conversion | `func_adcp_fold_determinism` after every change, treat any diff as blocking |
| `std::vector`/`std::string` throw `bad_alloc` on OOM vs. C's `NULL`-check-then-`stop()` | any container conversion | Accepted as-is: both are fatal-on-OOM in practice; document as the one sanctioned exception to "no exceptions," don't add `try`/`catch` anywhere else |
| RAII MPI/file guards don't run on the `stop()`→`exit()` abnormal path (`exit()` skips stack unwinding) | `main.c`, `error.c` | Not a regression — current code has the same gap. Don't oversell RAII as fixing abnormal-exit cleanup in commit messages |
| Installed static lib `adcp_core.a` ([src/CMakeLists.txt:33](src/CMakeLists.txt#L33)) — C++ name mangling breaks any out-of-repo consumer linking it directly | install target | No in-repo evidence of external linkage (ADCP is invoked as a subprocess per README); flag for the user to confirm with ADFRsuite maintainers in parallel, don't block on it |

## Relative sizing

Phase 1 (whole tree, 25 files) is the largest by file count, though most diffs are small casts — except `main.c`'s VLA fix, the one unavoidable real change. `params.c` and `peptide.c` idiomatic passes (Phase 2, steps 2-3) are the next largest by blast radius since `peptide.c` owns shared structs every other file reads. `energy.c` is large by line count but narrower in scope. `vdw.c` dispatch templating is small and isolated (4 call sites). `tools/` and the final CMake cleanup are trivial.

## Verification

**Always configure with `-DCMAKE_BUILD_TYPE=Release`, and build into a directory that has no stale objects.** Two false failures burned time during Phase 1 and both are traps for the next person:

- With no `CMAKE_BUILD_TYPE` the build is `-O0`, where a single `func_adcp_fold_determinism` run exceeds 400s against its 180s timeout. The test cannot pass in an unoptimized build; it takes ~54s in Release.
- Renaming `foo.c` → `foo.cpp` leaves the old `foo.c.o` behind in an existing build tree, and CMake keeps linking it. A `build/` tree carried across the renames ended up with 28 objects for a 17-source target — stale `-O0` C objects linked alongside the new C++ ones, producing garbage globals and a hang that looks exactly like a code regression. **`rm -rf` the build directory after every rename step.**

Run `cmake -B build -DADCP_TOOLS=ON -DCMAKE_BUILD_TYPE=Release && cmake --build build && ctest --test-dir build --output-on-failure` after every commit. For checkpoints listed above, additionally: `cmake -B build-dock -DADCP_DOCKING_TESTS=ON && cmake --build build-dock && ctest --test-dir build-dock -L "docking|validation" --output-on-failure` (first run fetches the ~10MB 3Q47 target once, cached thereafter).
