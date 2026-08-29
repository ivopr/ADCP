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

> **"No behavior change" turned out to be false, and stayed false for 60 commits.**
> `5660a33` silently reduced five vector normalisations in `energy.cpp` from double to
> single precision, because `sqrt` on a `float` expression binds to `double sqrt(double)`
> in C but to the `float` overload in C++. It changed every docking run past 200,000
> steps. See "Measured against pristine upstream" below. A rename is not automatically
> behaviour-preserving in this language pair — **verify it, do not assert it.**

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
- **MPI/OpenMP**: no `extern "C"` needed, both headers are C++-safe as-is. Verify empirically via the new CI legs, no code change expected. **This held for OpenMP but not for MPI** — see "Phase 2 progress — MPI unblocked" below: `flex.h` had `#include <mpi.h>` *inside* an `extern "C"` block, which does not compile once `mpi.h` pulls in C++ STL headers. Not caught earlier because nothing had compiled `flex.h` under real `PARALLEL` until then.

Ordering within Phase 2, by risk × blast radius:

1. **Done.** Leaves: `error.c`, `vector.c`, `rotation.c`, `random16.c`, `aadict.c`, `canonicalAA.c`
2. **Done.** `params.c` — do before anything downstream depends on its idiomatic shape
3. `peptide.c` — owns the shared structs. **Split into 3a/3b/3c** (see the step 3
   section below); doing it as one step means changing 13 files at once, which is
   not bisectable against the `energy.cpp` precedent.
   - **3a. Done.** `peptide.cpp`-local cleanup, no header or caller change.
   - **3b. Done**, after adding the NS/checkpoint tests it needed (3b-0). No
     `malloc`/`realloc`/`free` of these three types remains anywhere.
   - **3c. Done, deliberately scoped**: `erg`, `ergt`, `distb` are `std::vector`.
     `aa`/`aat` and the four `triplet` members stay raw pointers **by decision, not
     omission** — see "Why 3c stops where it does" below. The leaks a `Chain`
     destructor would have fixed were fixed by hand instead.
4. **Done.** `main.c` idiomatic pass
5. **Done.** `vdw.c` dispatch templating — measured 5.6% faster, bit-identical
6. **Done, cold paths only.** `energy.c` — the hot allocations, the three
   `ADenergyNoClash` VLAs and the whole `tc`/`View3` region stay as they are, on
   purpose; see the step 6 section below.
7. **Done.** `probe.c` — all 32 diagnostic mask bits swept under ASan, three real
   bugs fixed. `metropolis.c` — dead code only; the per-move path is untouched.
8. **Done.** `flex.c`, `checkpoint_io.c` — dead code and two reproduced
   buffer overflows fixed. The MPI half, once blocked on installing MPI, is now
   unblocked — see "Phase 2 progress — MPI unblocked" below.
9. **Done.** `nested.c` — dead FAST branch removed and a real nested-sampling
   bug fixed. MPI half unblocked the same way.
10. **Done.** `tools/*` — three bugs fixed, all sanitizer-reproduced. This step
    also found that the test suite's ASan gate was being masked; see below.
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

## Phase 2 progress — step 3a (peptide.cpp local cleanup) done

Step 3 as originally written ("`peptide.c` — owns the shared structs") is the
widest change in the migration: converting `Chain`/`Chaint`'s members to
`std::vector` forces ~22 `malloc(sizeof(Chain))` sites, 9 array `realloc` sites
and ~20 `= NULL` initializers across 13 files to change in one commit — the same
heap-container pattern that produced the still-unexplained `energy.cpp` hang, on
the hottest struct in the tree. So it was split; **3a is everything inside
`peptide.cpp` that needs no header change and no caller change.**

Six changes, four commits, each hammered 12× before the next:

- **Deleted dead `repair()`** (197 lines, 3 `malloc`, 3 `free`). `static` and never
  called — the live path is `initialize` → `repair_multichain` → `repair_segment`.
  It was also the only place a later vector conversion would have been a real
  behavior change: `repair` read `diag[i]`/`cosn[i]` at indices it never wrote,
  where `repair_segment` zeroes its block explicitly.
- **Fixed an `xaa_prev` leak.** `build_peptide_from_sequence` `malloc`'d over the
  pointer `allocmem_chain` had just allocated at the identical size, leaking
  `(Nchains+1)*sizeof(triplet)` per call. Safe to drop: `xaa_prev[1]` and
  `xaa_prev[chainid]` are both written before use and index 0 is never read.
- **Fixed a VLA stack overrun** in `repair_multichain`. `chain_starts`/`chain_ends`
  were VLAs (non-standard C++, and `CMAKE_CXX_EXTENSIONS` is `OFF`) sized
  `Nchains`, but the loop writes `chain_starts[next_chain]` *before* the check that
  `next_chain+1 == Nchains` — a malformed PDB smashed the stack and only then
  reported the inconsistency. Now `std::vector<int>` sized `Nchains+1`; every index
  expression unchanged. Phase 1 step 7 fixed `main.cpp`'s VLAs and missed these.
- **`str_without_separator`/`chain_ids` → `std::string`/`std::vector<int>`**, built
  with `push_back`. Drops 2 `malloc`, 2 pointless shrinking `realloc`, 2 `free`
  (one sitting ~150 lines away at the end of the function).
- **The three amino-acid list readers → one RAII helper.**
  `mark_fixed_aa_from_file` and `mark_constrained_aa_from_file` (×2) were
  near-identical `fopen`/`fscanf`/`fclose` blocks whose mid-loop `stop()` skipped
  the close. Now one `mark_aa_from_file()` holding the handle in a `unique_ptr`.
  **Also a real fix:** the range check was only `next >= chain->NAA`, so a negative
  index in the user-supplied list was an out-of-bounds write into `chain->aa`.
  Verified the previous build accepts `-3` silently and the new one rejects it.

**Deliberately not fixed: `aat_init`'s guards** (`src/peptide.cpp`). Both read
`sizeof(chaint)->aat`, which parses as `sizeof((chaint)->aat)` — i.e. `sizeof(AA*)`
== 8, a compile-time constant, not a capacity check. The conditions are therefore
always true, so **every `fulfill()` call reallocs all four `Chaint` buffers and
re-runs the field-copy loop, on the Monte Carlo hot path.** Left alone on purpose
and marked with a `ponytail:` comment: `Chaint` carries no size field to check
against, and changing hot-path allocation frequency in the same step as everything
else would make a determinism regression unbisectable. Fix in **3c**, where
`std::vector::resize` is the real no-op early-out for free.

Validated: 13/13 `ctest` and `func_adcp_fold_determinism` **12/12** after every
commit; ASan+UBSan clean including a 15200-atom fold through the new list-reader
path; `docking|validation` 3/3 with a live 3Q47 fetch (`val_3q47_redock` 114s).
The `fixed=` path has no test fixture, so it was checked by hand: stderr and every
ATOM record byte-identical against the parent commit.

### Notes for 3b and 3c — do not re-derive these

- **Before converting any member, convert the struct's allocation sites first.**
  This is the `params.cpp` lesson again. The single-object sites are `main.cpp:182,
  240, 978, 979, 984`; `flex.cpp:585, 587, 731, 892, 894`; `nested.cpp:94, 899,
  902, 906`; `peptide.cpp` (`pdbin`'s `tempchain`, which never escapes — make it a
  stack `Chain`); `vdw.cpp:1281`; `cdlearn.cpp:54, 634`; `checkpoint_io.cpp:367,
  370, 636, 639`.
- **Highest-risk sites are the growth-loop array reallocs** of live `Chain`s:
  `checkpoint_io.cpp:553`, `checkpoint_io.cpp:581`, `tools/cdlearn.cpp:669`. Each
  call bitwise-relocates every previously-constructed element. Also
  `nested.cpp:1107`, `flex.cpp:767`, `vdw.cpp:1280` (`malloc(sizeof(Chain)*28)`),
  `cdlearn.cpp:673, 686, 721`.
- **`AA` is POD and pointer-free** — its `vector` members are `double[3]` typedefs.
  It needs no conversion, and **3c must not touch it**: 8 whole-`AA` assignments in
  `metropolis.cpp` (362, 535, 693, 1094, 1227, 1364, 1511, 1775) and a `memcpy` in
  `tools/bfactor.cpp:142` all become UB the moment `AA` gains a non-trivial member.
- **`copybetween` is already a deep field-by-field copy**, not `*to = *from`, and
  all 20 callers go through it — including `main.cpp`'s replica-exchange pool,
  which never pointer-swaps. So the swap logic needs no change in 3c; only the
  `malloc`/`free` at `main.cpp:240`/`539`. (Its error message at the `Nchains`
  mismatch prints `to->NAA, from->NAA` — copy-paste bug, harmless, fix in 3c.)
- **`Chain::flex_data` is never initialized** at any creation site except
  `flex.cpp:730`, and `freemem_chain` never touches it. A real constructor in 3c
  fixes this for free — a migration win, not a separate fix.
- `allocmem_chain`'s `NAA*NAA*sizeof(double)` and `aat_init`'s
  `5*NAA*NAA*sizeof(double)` are `int` products that overflow above `NAA ≈ 46k`/`20k`.
  `std::vector` with a `size_t` product removes this.
- **`static char line[83]` in `getaa` is not a mechanical `std::string` swap.** It
  is cross-call parser state by design — the `do {} while (fgets(...))` loop
  inspects the line left over from the previous call, which is how the
  ENDMDL/END/TER sentinels and the residue-boundary return work, and the direct
  `line[30..37]` indexing relies on the array being 83 bytes regardless of the
  actual line length. Leave it unless you are building a `PdbReader` that owns
  `FILE*` + the line buffer.

### Testing note

`tests/run_fold_test.sh` does `rm -rf "$WD"` on a fixed workdir, so **two `ctest`
runs against the same build directory will wipe each other's scratch space** and
produce a spurious determinism failure. Hammer sequentially, one `ctest` at a time.
`--repeat until-fail:12` takes ~10.5 min, which exceeds some command timeouts —
two batches of 6 is equivalent.

## Phase 2 progress — steps 3, 4 and 5 done

**3b-0: the riskiest code had no tests.** Nested sampling runs only under `-n`, and
no test passed `-n`, so `nested.cpp` (1325 lines) and `checkpoint_io.cpp` (991
lines) had zero coverage — including `store_chain`'s realloc-grow loop over the
live NS population, the single highest-risk site in step 3. Two tests now drive it
by feeding a multi-model PDB (`tests/data/output.pdb` already has 1000 `MODEL`
records; the scripts truncate to 20, and `read_in_from_pdb` stores one `cpoint`
per model):

- `func_adcp_ns_smoke` — all 20 snapshots stored, same seed gives bit-identical
  NS output, energies finite and **monotonically decreasing** (the defining NS
  invariant, asserted instead of a stored baseline because these energies are
  post-MC and therefore `rand()`-dependent — the same reason `run_fold_test.sh`
  stores no baseline energy).
- `func_adcp_checkpoint_roundtrip` — checkpoints written, `-R` resumes at the
  checkpointed iteration rather than restarting at 0, two restarts agree. It
  deliberately does **not** compare a restart against an uninterrupted run: the
  checkpoint stores the population but not the RNG state, so those legitimately
  diverge. Verified empirically; do not "fix" that by asserting equality.

**Known limit of both, verified by mutation:** perturbing a stored chain by 0.001
still passes. They catch crashes, miscounts and nondeterminism, not small
deterministic corruption. That mutation *does* shift every value in the
`NS energy sequence:` line the NS test prints, so **every conversion in this area
must also diff that line against the parent commit on the same machine.** It was
byte-identical across all five commits below.

**3b: no `malloc`/`realloc`/`free` of `Chain`, `Chaint` or `Biasmap` remains.**
That is the greppable exit criterion, and it holds. Single objects became
`new T{}`/`delete` (value-initialization also closes the `flex_data`/`ll`/
`Nchains` holes the hand-nulling never covered); `pdbin`'s `tempchain` became a
stack `Chain`. Arrays split by whether they grow: the growth loops became
`std::vector` (`cpoints` — so `Chain **cpoints` is now `std::vector<Chain>&`
across four `checkpoint_io.h` signatures plus `new_amplitude`; and cdlearn's
`all_chains`/`all_chaints`/`all_biasmaps`), while fixed-size allocations became
`new[]`/`delete[]` (`input_chains`) or `resize` (`chaincopies`,
`all_chains_sim`). ~20 consumers taking a plain `Chain*` were left alone and get
`.data()`.

**3c so far: `Chain::erg`, `Chaint::ergt` and `Biasmap::distb` are `std::vector`.**
The largest allocations in the tree (`erg` is NAA², `ergt` is 5·NAA²). `resize`
also made the lengths `size_t`, fixing an `int` overflow above NAA ≈ 46000/20000,
and made `aat_init`'s `ergt` path a genuine no-op early-out.
`allocmem_chain`/`freemem_chain`/`freemem_chaint` were kept as shims, so their ~40
call sites did not change.

### Why 3c stops where it does — a decision, not an omission

`Chain::aa`/`Chaint::aat` and the four `triplet` members (`xaa`, `xaa_prev`,
`xaat`, `xaat_prev`) remain raw pointers. Converting them was costed and declined:

- **`aa`/`aat`** would need `getpdb`'s `AA**` signature reworked (plus
  `tools/bfactor.cpp`), 17 `.data()` sites and **164 pointer-arithmetic sites**
  rewritten, 82 of them in `metropolis.cpp`'s hot Monte Carlo path.
- **The `triplet` members** cannot go in a `std::vector` at all —`triplet` is
  `double[3][3]`, not a valid element type. The `std::array` equivalent indexes
  identically for `xaa[i][j][k]` but does **not** decay to `double(*)[3]`, and
  **117 sites pass `xaa[i]` straight into a `triplet` parameter** (`acidate`,
  `casttriplet`, `transset`, …). It would need a wrapper struct with a conversion
  operator. Note also that `matrix` and `triplet` are the same type, and
  `peptide.cpp` passes `xaa[0]` into a `matrix` parameter — so the wrapper could not
  even be a distinctly-typed one.

Converting **both** is the only thing that would have paid for the work, because
only then can `Chain` have a destructor, and the destructor is what would have
retired the manual `freemem_chain`/`freemem_chaint` protocol and the leaks. That is
~300 mechanical edits in the file whose sibling produced the still-unexplained
`energy.cpp` hang. **The leaks were fixed directly instead** (see below), which is
where the actual safety value was.

Leaving them raw is safe as things stand: `Chain` has no destructor, so its implicit
move copies those pointers exactly as the old bitwise `realloc` move did, and 3b
removed every `realloc`/`memcpy` of a `Chain`. Consequence to know about:
`aat_init`'s `sizeof(chaint)->aat` guard is still broken and still always true (it
parses as `sizeof(AA*)` == 8), so `aat`/`xaat` are still realloc'd on every
`fulfill()`. `ergt` no longer cares — `resize` is a real no-op.

### Steps 4 and 5

**Step 4 (`main.cpp`)** was small: the dead `swapFile`/`swapname` deleted, three
`sprintf`-into-buffer filenames to `std::string`. Deliberately unchanged, so nobody
redoes the analysis: `checkpoint_filename[256]` stays a char buffer because it is an
`sscanf` target, and **no `FILE*` RAII wrapper was added** — `fptr1` is already
`fclose`d, `fptr`/`fptr_pdb` are function-`static` and process-lifetime, and
`sim_params->infile`/`outfile` are owned by `sim_params`.

**Step 5 (`vdw.cpp` dispatch) was measured before being written**, because it is a
performance change on a hot path and the risk register calls float-accumulation
changes blocking. Fixed fold workload, median of 5, spread under 1%:

| build | time |
|---|---|
| function pointer (before) | 11.17 s |
| calls devirtualized only | 10.97 s |
| **bodies templated on the potential** | **10.54 s — 5.6% faster** |

ATOM records byte-identical throughout, so the accumulation order did not change.
Specializing the whole body is what buys most of it; merely removing the indirect
call does not. `probe.cpp`'s `vdw_fn` is left alone — this document lists it as
firing every MC step, but `vdw_contributions` is commented out of the test dispatch
table and is unreachable.

That work surfaced a **latent linkage bug**: `vdw_lj` and `vdw_hard_cutoff` were
defined `inline` in `vdw.cpp` while `include/vdw.h` declares them without it.
`probe.cpp` calls them through a function pointer and needs a real symbol, which
existed only because `vdw.cpp` happened to take their address. Once the callers
became templates, every use inlined, no out-of-line copy was emitted and the link
broke. Fixed by matching the header.

### The MPI build was broken and had never been compiled

`MPI_Comm *MPI_COMM = mpi_comm;` assigns a `void*` to a `T*` implicitly — legal in
C89, illegal in C++ — at **7 sites** in `nested.cpp` and `checkpoint_io.cpp`, left
over from Phase 1. So `-DADCP_MPI=ON` could never have built, consistent with this
document's own note that `adcp_mpi` was never compiled during the migration. Fixed
with explicit casts. That immediately paid for itself: with the PARALLEL files
compiling, the very next check caught `mpi_send_chain`/`mpi_rec_chain` passing the
newly-vectorized `erg` where a `void*` was expected — a break no other gate here
would have found.

MPI is not installed on the dev machine, so the PARALLEL paths are verified with
`g++ -std=c++17 -fsyntax-only -DPARALLEL -Iinclude -I<stub>` against a hand-written
stub `mpi.h` (~20 symbols). **That checks syntax and types, not linking or
behaviour** — CI's `mpi` job remains the real gate, and it should now actually get
past the compiler.

### Pre-existing leaks — fixed directly

Since 3c stops short of a `Chain` destructor, these were fixed by hand rather than
waiting for RAII:

- `flex.cpp`'s `ns_for_flex_processor` never released `chain` or `chaint`
  (`finalize_flex` releases `flex_data`/`input_chains`/`rmsd`, not the chain).
- `nested.cpp`'s `new_amplitude` has a `if (P == 1) { … return; }` early exit that
  skipped the cleanup at the end of the function, once per recalibration.
- `nested.cpp`'s `only_output_checkpoint` exit had its per-element `freemem_chain`
  loop commented out. `cpoints` is a `std::vector<Chain>` now, so its destructor
  releases the *element storage* — but `aa`/`xaa`/`xaa_prev` are still raw and
  still need the explicit loop. Do not delete it thinking the vector covers it.
- `probe.cpp`'s `initialize_displacement` released `chain_init` but not the stack
  `Chaint` or the `displacements` chain, **and** leaked `G_to_delete`/`G2_to_delete`
  — a fifth leak nobody had listed, found by running the path rather than reading it.

**That last one is the lesson.** `initialize_displacement` runs only under test mask
`0x2000000` and no test sets it, so ASan had never seen it. Running it by hand
reported `10400 byte(s) leaked in 200 allocation(s)`; after the fix, nothing. The
other mask-gated `probe.cpp` diagnostics have still never been run under ASan and
are a reasonable place to look for more of the same.

Two of the four are `#ifdef PARALLEL` and are verified only by the stub-`mpi.h`
syntax check, never by running.

## Phase 2 progress — step 6 (energy.cpp, cold paths only) done

`energy.cpp` splits cleanly into cold setup/diagnostic code and hot per-MC-move
code. **Only the cold half was touched.** Nothing in this step puts a heap
container, or any new indirection, on a per-move path — which matters more than it
used to, because of what step 6-0 turned up (see "MECHANISM RESOLVED" below).

**Real bugs fixed:**

- **Three heap overflows in the setup parsers.** `gridmap_initialise`,
  `transpts_initialise` and `ramaprob_initialise` each size a buffer from a count in
  the file header, then write one element per input line without ever comparing the
  two. A file with more data lines than its header declares silently overflows the
  heap. All three now stop and warn; they are cold parse paths, so the check is free.
- **A dead `fopen` retry** in `biasmap_initialise`: on failure the old code re-ran the
  identical `fopen` with identical arguments. It cannot succeed the second time, so
  only the `!= "NULL"` test it was ANDed with ever did anything.

**Cleanups:**

- The four `#ifdef LJ_HBONDED_HARD` blocks in `all_vdw`. `LJ_HBONDED_HARD` is never
  `#define`d anywhere and no CMake target sets it — only `PARALLEL` is — so this was
  dead in every configuration. It also carried a bug no compiler ever saw:
  `((i1>=1) || (i1<chain->NAA))` uses `||` where `&&` is meant, so the bounds check is
  vacuous. Deleted rather than repaired. **`vdw.cpp`'s own `LJ_HBONDED_HARD` block is
  still there** — it guards `exclude_hard`, which step 5 templated.
- Twelve commented-out `malloc`/`free` corpses, including two copies of
  `//free(tc),…` left from the heap-`tc` version.
- The cold VLA in `energy_matrix_calculate` → `std::vector` (docking-only code).
- `biasmap_initialise`'s `FILE*` → `unique_ptr`, collapsing three repeated
  `if (fin) fclose(fin);` blocks. Readability only — the old pairing was correct.
- **The last `goto` is gone.** MIGRATION.md's other scheduled site (old
  `energy.c:2876`) had already been removed by the `distb` conversion.
- `energy.h`'s `#ifdef __cplusplus` fork around the `scoreSideChain*` prototypes:
  the `#else` branch held C99 VLA-parameter declarations and has been dead text since
  the last `.c` file went away.

### The `goto` had one trap worth recording

The outer loop must `break` rather than fall out through its condition. Failing the
condition runs `i++` one extra time, and `i` is read *after* the loop to report how
many contact-map rows were actually read — so getting it wrong silently turns
`Go-type bias: 6x12` into `7x12`. The flag version breaks out of both loops.

### How the uncovered code was validated

`biasmap_initialise`'s contact-map branch has **no test coverage** — every test passes
`Bias=NULL`, which returns before the read loop. So groups (c) and (d) were validated
by hand: a generated 12×12 contact map, plus a deliberately truncated one that drives
the EOF path the `goto` used to serve. Both produce byte-identical stderr and ATOM
records across the change, and both are clean under ASan+UBSan, which had never seen
this code. That harness is worth rebuilding if anyone touches `bias()` — it is ~200
lines that no test reaches.

### Deliberately not done, and why

- **The five hot allocations** — `sbond_energy`'s `cyslist`/`cyspos`/`cysdist` and
  `secondary_radius_of_gyration`'s `s_com`/`weights`. Biggest RAII win in the file
  (8 `free`s and a quadratic realloc loop) but they run per MC move.
  `secondary_radius_of_gyration` additionally has zero coverage — it needs
  `srgy_param`/`hphobic_srgy_param` nonzero.
- **The three `ADenergyNoClash` VLAs**, coupled to `metropolis.cpp`'s matching VLA on
  the same call, and the whole `tc`/`View3` region. (One exception has since been made
  here: the five `sqrt((double)…)` casts that restore the C precision semantics. They are
  a correctness fix, not a conversion — **do not remove them**, see "Measured against
  pristine upstream".)
- **The ten file-scope setup allocations stay raw pointers.** The plan called for
  `std::vector`, but they are already freed at `main.cpp:1099-1107`, and
  `Xpts`/`Ypts`/`Zpts` are read from `metropolis.cpp`'s `transmove` and
  `peptide.cpp` — a move path. The actual defect was the missing bounds checks, and
  those are fixed without touching a type or a header.
- **The 19 remaining `energy.h` declarations with no external callers stay public.**
  Making them `static` changes inlining, and most of them (`bias`, `hbond`,
  `hydrophobic`, `electrostatic`, `sidechain_hbond`, `stress`, `proline`) run per MC
  move. Tidying a header is not worth codegen churn in this file.
- **Latent, documented, not fixed:** casting `±inf` to `int` in `getindex` is UB and
  runs on every fold move — benign today only because of the `< 0` check right after.
  It is in the innermost function in the program.

### Coverage holes in energy.cpp, for whoever comes next

`bias()` and the contact-map reader (~200 lines); the gradient/finite-difference block
(~300 lines, needs `-t 2000`); `energy2cyclic`/`cyclic_energy`;
`secondary_radius_of_gyration`; and the non-trivial `external`/`external2` bodies —
all have zero automated coverage. The grid-map and `scoreSideChain*` code is reached
**only** by the opt-in 3Q47 docking tests.

## Phase 2 progress — step 7 (probe.cpp, metropolis.cpp) done

Two files needing opposite treatment: `probe.cpp` is entirely cold diagnostics,
`metropolis.cpp` is the per-move Monte Carlo engine.

### probe.cpp — the mask-bit sweep found three real bugs

`tests()` maps table entry N to mask bit N, and the default masks are `0x8803`
(MC) / `0x18803` (NS) — so **only 5 of 32 bits had ever executed**, and 27 had
never been under a sanitizer. Ran a fold once per bit under ASan+UBSan on a
20-residue peptide. **Three bits reported; 29 were clean** — the clean result is
worth recording too, because it bounds how much is left to find here.

| bit | diagnostic | finding |
|---|---|---|
| `0x1000000` | `number_of_contacts` | two `calloc(NAA)` never freed, on any path. ASan: 1680 bytes × 20 objects, twice — once per output snapshot, so it accumulated. Now `std::vector<int>`. |
| `0x8000000` | `CA_geometry` | heap-buffer-overflow **read**. The `aa[1..16]` block had no guard at all; the second block reads `aa[24]` while guarded by `NAA > 17`, so it overran on anything under 25 residues. Guards are now `> 16` and `> 24`. |
| `0x40000000` | `vdw_max_gamma` | leaked the `->sequence` that `sim_params_copy` had just allocated. |

**The third one is the interesting one.** `copy_string` deliberately does *not*
free its target — the free is commented out with *"would crash things in
model_params_copy"* — so `build_peptide_from_sequence` overwriting
`my_sim_params->sequence` orphans the previous copy. Fixed locally in
`vdw_cutoff_distances_calculate`, matching the `->seq` line directly above it.
`copy_string` itself was left alone: making it free would touch every caller in
the tree, and the comment says that was already tried. **Any other
`copy_string` into an already-populated field is the same bug** — that is the
generalisable finding, not the single site.

It only bites when the function runs *after* `->sequence` is populated, which is
why `main.cpp`'s startup call has never shown it.

Suspected by reading but **not** reproduced by any sanitizer, so documented rather
than patched: shift UB in `hpattern` (`1 << (15 - j + i)` can go negative or past
31), uninitialized `nca`/`cacb` in `all_torsions`' `G2_`-without-`G__` branch
(MSan territory, not installed), and `hbm` left dangling after `free`.

Also deleted six diagnostics that had a definition, a header declaration and a
*commented-out* table entry — `vdw_contributions`, `cm_ideal_4`, `cm_alpha_8`,
`cm_native_go`, `fasta`, `hbond_pattern`. **Bit numbering is the thing that could
have broken there**; it did not, verified by count and by re-running the three
bits the sweep had reported on.

### metropolis.cpp — ~400 dead lines, and nothing else

`crankshaft_adk()` (262 lines, zero callers) and `flipChain()` (130 lines, zero
callers), plus 11 commented `//free(ADEnergy_Chaint)` fossils sitting next to the
stack VLAs that replaced the heap version, and the never-used `reject` static.

**Nothing on the per-move path was touched** — `allowed()`, `crankshaft()`,
`crankshaftcyclic()`, `move()`, `reModNum()`, `indMoved()`, and every VLA. In
particular `crankshaft()`'s 113-line one-shot table build stays inline: extracting
cold code from a hot function is the textbook improvement and is exactly the
layout edit step 6 proved can flip the hang rate.

### transopt's `mod` parameter is inert — measured, deliberately left inert

```c
int maxStep = 10;  int maxNoImprovStep = 3;
if (mod == 1){ int maxStep = 30; int maxNoImprovStep = 5; }   /* shadows, dies at the brace */
```

`main.cpp` passes `mod == 1` from **six of its seven** call sites, so "hard
minimization" has never differed from the soft mode. Making it real was measured
on the 3Q47 redock rather than assumed:

| | val_3q47_redock | dock_3q47_smoke | RMSD | top targetE |
|---|---|---|---|---|
| shipped (inert) | 130 s | 8.6 s | **0.72 Å** | **-28.95** |

Note: those RMSD/targetE figures were measured while the `sqrt` precision defect
described under "Measured against pristine upstream" was present. Post-fix the same
run gives 0.66 Å / -28.5508. The conclusion (`mod` is inert, the change is not worth
shipping) is unaffected — both arms were measured on the same defective build.
| `maxStep = 30` | 156 s | 19.3 s | 0.80 Å | -28.34 |

**+20% on the validation run, 2.25× on the smoke test, for no quality gain** — the
RMSD and energy differences sit inside the noise of a 16-run stochastic search and
if anything favour the current behaviour. So the change was **not** shipped: ADCP's
published results were produced with this budget, and whether hard minimization
should finally do something is a science decision for the ADFRsuite maintainers.
The shadowing is replaced by a comment carrying these numbers, so the next reader
finds the measurement instead of "fixing" it silently.

### Latent, documented, not fixed

`transmove`'s `while (vecind2 == vecind1)` loop spins forever when `chain->NAA == 2`
(one residue): `rand()%1 + 1` is identically 1. It needs `NAA >= 3` to terminate
and `external_potential_type == 5` to be reached. Fixing it means adding a branch
to a warm function, which is not worth it for a one-residue peptide.

## Phase 2 progress — steps 8 and 9 (partly done; the MPI half is blocked)

### The MPI build has never compiled ANY of the MPI code

`src/CMakeLists.txt` makes `PARALLEL` **`PRIVATE` to the `adcp_mpi` executable**,
which is only `main.cpp`. `adcp_core` — containing `nested.cpp`,
`checkpoint_io.cpp`, `flex.cpp`, `probe.cpp` — is compiled **once, without
`PARALLEL`**, and linked into both `adcp` and `adcp_mpi`. Consequences:

- `adcp_mpi` is an MPI `main` linked against a **serial `nestedsampling()`**. Every
  rank would run the same complete serial simulation.
- `mpi_send_chain` / `mpi_rec_chain` are **never emitted into any object file**.
- **~1150 lines of MPI code are unbuilt by every configuration**, `-DADCP_MPI=ON`
  included.

This is why step 6's fix (7 implicit `void*`→`T*` conversions) only got `main.cpp`
compiling — nothing else was ever being compiled with `PARALLEL` at all. **Every
"validated on MPI" or "the CI mpi job covers it" claim in this document, past or
future, should be read against this.** The CI `mpi` job builds a binary that cannot
work.

**Blocked, not fixed.** Fixing it properly means compiling `adcp_core` twice (an
`OBJECT` library, or a second static lib with `PARALLEL` and `MPI::MPI_CXX`), and
verifying that requires a real MPI. `libopenmpi-dev` is available in apt but not
installed, and installing it needs `sudo`. Until then the stub-`mpi.h`
`-fsyntax-only` check remains the only local gate — it checks syntax and types, not
linking or behaviour.

### A real nested-sampling bug

The serial live-point draw was `copies = 1 + (rand*N) % N`, i.e. `1..N`. But
`find_worst` swaps the heap root out to `chainhash[heaplength]` and decrements, so
with `P == 1` **the discarded worst point sits at index `N`** and the live heap is
`1..N-1`. So with probability `1/N` the new sample was seeded from the point being
discarded — nested sampling's validity rests on drawing from the *surviving* prior
mass. Confirmed by instrumentation: a run drew `copies==20` with `N==20` whose `ll`
equalled `logLstar` exactly, and `logLstar` is by definition the discarded point's
likelihood. The PARALLEL branch already did it right with `% (N-P)`.

**The NS energy sequence baseline changed as a result**, deliberately:

```
before  40.058829 31.225769 10.764044 7.463047 6.393763
after   40.058829 688.128472 81.995422 10.222738 7.086147
```

### The evidence estimate is now checked against a closed-form answer

Everything else in the suite asserts that ADCP does not crash, is reproducible, or lands
in a plausible range. Nothing checked that nested sampling computes the **right number** —
the peptide likelihood has no closed form to check against, so there was nothing to
compare to.

`tests/ns_evidence_test.cpp` (ctest `num_ns_evidence`, label `functional`) closes that.
It links `adcp_core` and drives the real `find_worst()` / `update_NS_parameters()` /
`alpha = exp(-1/N)` bookkeeping against a problem that does have one:

    prior       x ~ Uniform(0,1)
    likelihood  L(x) = exp(-x/tau)
    evidence    Z = tau * (1 - exp(-1/tau))

`L` is monotone in `x`, so `{L > L*}` is exactly `{x < x*}` and a constrained sample is
`Uniform(0, x*)` drawn **exactly**. That is deliberate: removing MCMC quality from the
measurement leaves the bookkeeping as the only thing under test. Measured, N=100,
3000 iterations, four fixed seeds: errors of -0.02 to +0.14 nats against
`log Z = -1.6161986619`, all inside 3 sqrt(H/N), mean inside 2 sigma.

Three properties make it a better gate than anything else here:

- **It asserts a value, not self-consistency.** The first test in the repo that does.
- **It is libc-independent.** It uses its own `std::mt19937`, not `rand()`, so unlike
  every other test it is portable enough to carry an exact expected number — the reason
  `run_fold_test.sh` and `run_ns_test.sh` deliberately store none.
- **It guards this section's fix, and proves the guard is not vacuous.** Both draw
  expressions are evaluated from the same uniform and counted: the shipped
  `% (N-1)` selects the discarded point at index `N` **0** times, the pre-`35fb3fb`
  `% N` would have selected it **122** times in 12000 draws — 1.02%, matching the
  predicted `1/N`. A guard that cannot fail is not a guard, so the test asserts both
  halves. It also asserts the layout invariant the fix rests on: after `find_worst`,
  `chainhash[N].ll == logLstar` and no survivor shares that index.

Mutation-tested, since a test that has only ever passed proves nothing. Scaling
`log_DeltaX` by 1.02 in `update_NS_parameters`, and taking `logLstar` from the heap root
instead of `chainhash[N-P+1]` in `find_worst`, are both caught.

**Known limit:** exact constrained sampling means the seed cannot influence the evidence
here, so the numerical half does not detect a reintroduced seed bug — the draw counters
do. Detecting it numerically would need a finite-step MCMC model, which trades a sharp
deterministic check for a noisy statistical one. Not worth it.

### …which exposed that the NS test's own assertion was wrong

`run_ns_test.sh` asserted the energies decrease **monotonically**, described in step
3b-0 as "the defining invariant of nested sampling". It is not an invariant of this
implementation. Measured on *unmodified* code, varying only the seed:

| seed | sequence | monotonic |
|---|---|---|
| 12345 | 40.06 31.23 10.76 7.46 6.39 | yes |
| 99 | 40.06 13.14 **13.69** 7.46 6.39 | no |
| 4242 | 40.06 13.14 **30.20** 7.46 6.39 | no |
| 555 | 40.06 **313.37** 13.14 10.22 7.46 | no |

3 of 8 arbitrary seeds violate it. It had been passing only because the seed is
pinned at 12345 — so the test was one seed change from failing on correct code, and
it **actively vetoed the sampling fix above** on evidence that was an artefact of its
own strictness. Replaced with a property that does hold on 15 of 15 seeds: the last
reported energy is below the first. **Lesson: an invariant asserted from theory and
confirmed on one seed is not confirmed.**

### Unreachable feature branches deleted

- **`nested.cpp`'s `FAST` branch** — `FAST` is `#define`d nowhere in the tree or in
  any CMake file, so the `#else` of `#ifndef FAST` (241 lines: a second
  `collect_chains`/`return_and_reheap_chains` pair, the whole instruction-set
  machinery) plus two `#ifdef FAST` blocks were unreachable in every configuration.
  That orphaned the `Instructions` type: the surviving functions took an
  `Instructions*` they never referenced and both call sites passed `NULL`.
- **`flex.cpp`'s `setup_bias_hbonds`** — gated on `flex_params.only_bias_hbonds == 1`,
  which is initialised to 0 twice and assigned nowhere else. It was the only writer
  of `FLEX_data::Hbond_aaH`/`Hbond_aaO`, so deleting it also removed 12
  grow-one-element-at-a-time `realloc`s of the `p = realloc(p, n)` form that leaks
  and then dereferences NULL on failure, plus both struct members.

### Two buffer overflows, both reproduced

- **`flex.cpp`, seven `char[256]` path buffers.** `output_path` holds up to 255 chars
  (`"%255[^,]"`), and every site concatenated a prefix or suffix into 256 bytes —
  worst of all `create_directory` (`"mkdir -p " + path`), which then goes to
  `system()`. ASan on a 250-char path: `stack-buffer-overflow, WRITE of size 260`.
  All seven are `std::string` now. Reachable without MPI via `-t 40000 -p FLEX=1,…`,
  which is how it was tested — no suite test covers this code.
- **`checkpoint_io.cpp`'s `read_checkpoint_header`.** `fscanf(..., "%s\n", seq)` with
  no field width, into `malloc(NAA+2)` where `NAA` came from the same file's header.
  A hand-corrupted checkpoint gives `heap-buffer-overflow, WRITE of size 401`. Now
  width-limited and length-checked, plus a `NAA >= 1` guard and a `free` of any
  previous `seq` that the old code leaked.

### Still outstanding here

Not done, and worth naming: `checkpoint_io.cpp`'s loop bound
`cpoints->aa[NAA-1].chainid` is a value read *from the file* used to index
`xaa_prev`, and the consistency check that would catch it runs *after* the read loop.
The `sprintf` into `malloc(1010)` sites from argv-controlled paths. `nested.cpp`'s
`fclose(fptr)` with no NULL check and two unchecked `fopen`s. And everything gated
behind getting a real MPI.

## Phase 2 progress — step 10 (tools/) done, and the ASan gate was blind

### The test suite was hiding sanitizer failures

**CTest's `PASS_REGULAR_EXPRESSION` overrides the exit code.** If the expected
string is printed, the test passes even when the process then aborts. All 11
`smoke_*` tests use one — so **an ASan or UBSan failure in any `tools/` binary was
invisible to the suite.**

Not hypothetical: `dssp2cm` leaks under ASan on `/dev/null`, the exact input
`smoke_dssp2cm` feeds it, and the test reported PASS. **Every "ASan suite clean"
claim in this document has therefore been partly blind for `tools/`.** What still
stands: the `func_*` tests (their shell scripts assert and exit non-zero), and every
`src/` finding that was reproduced by running a binary directly — which is how all
of them were actually found.

Fixed by applying `FAIL_REGULAR_EXPRESSION` for sanitizer output to every test. It is
evaluated independently of the pass pattern, so it needs no restructuring and is a
no-op in the normal build. The `PASS_REGULAR_EXPRESSION`s stay: several smoke tests
pass *because* the binary calls `stop()`/`exit(EXIT_FAILURE)`, and the pattern is
what separates "failed correctly" from "crashed".

**`WILL_FAIL` inverts `FAIL_REGULAR_EXPRESSION` too**, which would have re-masked
`smoke_mergie`. That test had no output assertion at all — it passed on any non-zero
exit, and was in fact passing because a `malloc((size_t)-8)` failed. It now asserts
the usage message.

### Three bugs, all reproduced under a sanitizer

Found by running each tool directly under ASan on its own fixture, not by trusting
ctest:

- **`dssp2cm`** — `parse_dssp_header`'s three `calloc`s for `map`/`seq`/`ss` were
  never freed. The leak `smoke_dssp2cm` was silently passing over.
- **`cdlearn`** — `realloc(..., strlen(retval+1))` on two adjacent lines where
  `strlen(retval)+1` is meant. `strlen(retval+1)` is `strlen(retval)-1`, and the
  `strcpy` on the next line writes `strlen(retval)+1` bytes: **a 2-byte heap
  overflow, twice**, reachable from a mistyped filename. glibc's fortify catches it
  first — `*** buffer overflow detected ***`.
- **`mergie`** — with no arguments `argv[1]` is `NULL` and `files = argc-2 = -1`, so
  `malloc(files * sizeof(FILE*))` requests `(size_t)-8`. ASan:
  `allocation-size-too-big 0xfffffffffffffff8`. Added the missing `argc` check.

### The checkpoint `chainid` overflow (leftover from step 8)

`read_checkpoint_entry` used the last residue's `chainid` — read straight from the
file — as the bound of a loop writing into `xaa_prev`, which holds only `Nchains+1`
triplets. The consistency check that would catch a mismatch runs *after* the whole
read loop.

Demonstrated by generating a real checkpoint and changing **one character** (last
residue's `chainid` 1 → 5, header still `Nchains=1`):
`heap-buffer-overflow, WRITE of size 8 in read_checkpoint_entry`. The write lands
because the energy matrix follows in the file, so `fscanf` keeps finding numbers.

Fixed **at the read, not at the loop** — `chainid` is validated against
`cpoints->Nchains` right after the `fscanf` that produces it. That also covers the
three sites in `peptide.cpp` that index `xaa_prev` with the same field, which a guard
on the loop alone would have left exposed.

### Found by reading, then confirmed and fixed — the tools/ survey closed out

The prior section flagged five suspects across five files with no sanitizer
reproduction; per this document's own rule ("do not patch on suspicion"), all five
were re-investigated with crafted inputs before touching any code. All **nine**
distinct findings across those five files reproduced, and are now fixed. A tenth,
more serious bug was found incidentally while building a repro for one of them.

**Sanitizer environment note, before the findings:** ASan's completeness turned out
to depend on the compiler. `cm.cpp`'s `s[num]` finding (below) did **not** reproduce
under `-DCMAKE_CXX_COMPILER` left at its default (`/usr/bin/c++`, i.e. g++) with a
bare `printf("%s", s)` — GCC's libsanitizer doesn't intercept the internal scan
inside `printf`'s `%s` handling as completely as LLVM's compiler-rt does. Confirmed
by three independent checks: (1) an isolated `malloc(n+1)`/`printf("%s")` pattern
compiled with `clang++ -fsanitize=address` catches it immediately, at `-O0` and
`-O2` alike (at `-O2`, the compiler folds `printf("%s\n", s)` into `puts(s)`, which
clang's ASan still catches); (2) adding an explicit `strlen(s)` call ahead of the
`printf` makes **g++'s** ASan catch it too, at the exact same address; (3) after the
fix, that same explicit `strlen(s)` probe is clean. So the bug was real and CI's
`sanitizers` job (which also builds with the default `c++` / g++) would **not**
have caught it — worth remembering for any future `printf("%s", ...)`-shaped
finding: a clean run there is not proof the string is not over-read.

**1. `cm.cpp`, `write_contacts` — `s[num]` never written before `printf("%s")`.**
`s = malloc(num+1)` but the fill loop only ever touches `s[0..num-1]`; the reserved
terminator byte was read straight off the allocator. Reproduced on the *existing*
`smoke_cm` fixture with no crafted input at all (see the sanitizer note above for why
this needed an explicit `strlen` probe to show up under the CI toolchain). Fixed:
`s[num] = '\0';` before the print.

**2. `cm.cpp`, `parse_input` — bound check ran after the write, and used `>` not
`>=`.** A ≥2731-residue PDB wrote `ca[2730]` (one past the last valid index, 2729)
before the too-big check ever ran. Reproduced with a synthetic 2735-residue PDB
(sequential `ATOM ... CA` records, generated, not committed — 2700+ lines is not
worth carrying in the repo); UBSan: `index 2730 out of bounds for type 'double
[2730][3]'` at the write. Fixed by checking `num >= 2730` immediately after each of
the two `num++` sites, before the corresponding write — which surfaced a second,
self-inflicted bug during the fix: an early `return num;` at the new check must
return `num`, not the old `num + 1`, since index `num` itself was never populated.
Verified: `ctest`-fixture output unchanged; the 2735-residue input now prints 2730
rows and a clean "This file is too big!" notice instead of crashing.

**3. `bfactor.cpp`, `parse_input` — `ava`/`sqa` sized from the first model only.**
`getpdb()` overwrites the global `NAA` on every call; `ava`/`sqa` are allocated once,
before the multi-model loop, from whichever `NAA` the *first* model reported. A later
model with more residues overruns both. Reproduced with a 2-model PDB (model 2 has
2 extra residues cloned from the fixture's last one) — real
`heap-buffer-overflow` in `subtract()` (called from `update()`) under ASan. Fixed by
rejecting the mismatch outright (`stop("bfactor: all models must have the same
number of residues")`) rather than attempting a numerically-dubious resize-and-
continue — averaging B-factors across models of different length isn't a
well-defined operation this tool was ever meant to support. Committed as
`smoke_bfactor_multimodel` (`tests/data/bfactor_multimodel.pdb`).

**4. `dssp2cm.cpp`, `parse_dssp_body` — `k`/`i`/`j` unchecked against `n_res`.**
All three come straight from a DSSP body record's fixed-width columns and index
`seq`/`ss`/`map`, which are sized off the header's residue count. Reproduced with a
3-line synthetic DSSP file whose one body record declares `k=99` against a
header-declared size of 4: `heap-buffer-overflow` in `parse_dssp_body` at the
`ss[k] = t` write. Fixed with a range check right after both `sscanf`s (`i`/`j` may
legitimately be `0`, the "no contact" sentinel — only genuinely out-of-range values
are rejected), skipping the malformed record. **That fix needed a second fix**: the
function's final terminator write (`seq[++k] = '\0'`) used the same `k`, so a
rejected record left a stale out-of-range value sitting there for when the input
ended. Added a separate `last_k`, updated only when a record is accepted.

**5. `dssp2cm.cpp`, `print_contacts` — `char str[11]` overrun on ≥4 contacts one
side.** The buffer reserves exactly 2 extra slots each direction (indices 8-9
forward, 1-2 backward) beyond the fixed ±1/±2 neighbours; nothing bounded the loops
that fill them. Reproduced with a 7-line synthetic DSSP giving one residue four
beta-sheet contacts on its forward side: UBSan, `index 11 out of bounds for type
'char [11]'`. Fixed by capping the fill loops at the buffer's real capacity
(`k < 10` forward, `k >= 0` backward) — extra contacts beyond that are now dropped
instead of overflowing, and as a side effect the terminator at `str[10]` can no
longer be silently overwritten by a 3rd forward contact either.

**6. `dssp2cm.cpp`, `parse_dssp_header` — `n_res`/`n_chains` uninitialised when the
header line is missing.** ASan doesn't catch use-of-uninitialized-value; confirmed
with MSan (`clang++ -fsanitize=memory` — no MSan-instrumented libc++ needed, this
tool has no STL surface) on the *existing* `smoke_dssp2cm` (`/dev/null`) input:
`MemorySanitizer: use-of-uninitialized-value` at the `calloc(n_res * n_res, ...)`
call. Fixed: `n_res = 0, n_chains = 1` at declaration, plus a stderr warning when the
header search comes up empty, so the empty-structure fallback is now deliberate
rather than accidental.

**7. `cdlearn.cpp` — `fscanf(list_file, "%s", ...)` with no field width into
`char[1024]`.** Reproduced with a `-L` list file containing one 2000-byte token:
ASan `stack-buffer-overflow` inside `scanf_common`. Fixed with `%1023s` (matches the
buffer-size-1 idiom already used at the 12 `%256[^,]`→`%255[^,]` sites from
`1e82869`), plus a second, related fix: the subsequent `strcpy`+`strcat` onto
`pdb_filename` (also 1024 bytes, `next_pdb_filename` + `".pdb"`) could still overflow
even with a correctly-bounded `next_pdb_filename` — replaced with a bounds-checked
`snprintf` that `stop()`s with a clear message instead. Committed as
`smoke_cdlearn_longtoken` (`tests/data/cdlearn_longtoken.txt`).

**8. `cdlearn.cpp` — `sprintf(index, "_%d", iter+iter_start)` into `char[10]`,
`iter_start` from `-I` on argv.** This was the one MIGRATION.md previously called the
most expensive to reproduce — it needed a minimal but *real* PDB + `.icm` contact
map + non-empty learn string, none of which exist anywhere in the repo, to reach the
vulnerable line. Built one: reused `tests/data/peptide12.pdb`, a hand-built
all-zeros `.icm` (the contact-map format is `(NAA-1)²` whitespace-separated
numbers/symbols — trivial once read from `biasmap_initialise()`), a one-line list
file, `-l V` (any single valid learn-string character suffices), `-I 100000000`.
The restart-file `sprintf` fires on the very first loop iteration (`iter==0`, since
`iter % 100 == 0`), so no real iteration count was needed. Reproduced first try:
ASan `stack-buffer-overflow`, `WRITE of size 11` into `index[10]`. Fixed: widened to
`char index[16]` (an `int` is at most 11 digits + sign + `_` + NUL = 13) and switched
to `snprintf`. Not committed as a permanent fixture — the full pipeline is more
setup than a single test file justifies; the recipe above is reproducible from
scratch in a couple of minutes if needed again.

**9. `ramachandran.cpp` — `n1[1]`/`n1[2]` uninitialised on the first `angle()`
call.** Only `n1[0]` is reset at the `start:` label; the first residue pair's
`eta = angle(n1, n2)` reads `n1` before line 202 ever writes it. Confirmed with MSan
on the *existing* `smoke_rama_pdb` fixture: `use-of-uninitialized-value` inside
`angle()` (`vector.cpp:308`), called from `ramachandran.cpp:195`. Fixed: reset
`n1[1]`/`n1[2]` to `NaN` alongside `n1[0]` — consistent with every other "no previous
residue" placeholder on that line, and the post-fix output for the first row is
`nan` for `eta` exactly like `phi`/`chi`/`chi2` already were.

**10. Found incidentally, not one of the original five: `cdlearn.cpp` never calls
`ramaprob_initialise()`.** Building finding 8's repro pipeline reached
`energy_matrix_calculate()` for the first time in this tool's history (matches this
document's own repeated note that nothing has ever tested `cdlearn` past its early
`stop()`s) and immediately segfaulted: `ramaprob`/`alaprob`/`glyprob` are only
`malloc`'d inside `ramaprob_initialise()`, which `main.cpp` calls before any energy
calculation but `cdlearn.cpp` never calls at all — a guaranteed NULL-pointer read in
`ramabias()` on **any** real CD-learning iteration, in every build, forever. Fixed by
adding the same `ramaprob_initialise()` call `main.cpp` makes, in the same relative
position (right after `model_param_read`, before `initialize_sidechain_properties`).
This does mean `cdlearn` now genuinely requires `ramaprob.data` in its working
directory like `adcp` always has — `smoke_cdlearn` needed `WORKING_DIRECTORY
${TEST_WORK}` added (that directory already gets `ramaprob.data` copied into it for
the `func_adcp_*` tests) to keep passing.

**All three now closed out:**

**`ramachandran.cpp`'s `goto` is gone.** The `start:` label's body (two statements:
reset every derived angle to `NaN`) was reached two ways — the `goto start` for the
very first residue pair, and normally as the `else` branch's tail on a later chain
break. Both reached the identical code, so the goto was replaced by duplicating those
two statements in place of it and deleting the now-unreferenced label — the standard,
behavior-preserving fix for a "skip straight to this block" goto, not worth a helper
function for two lines used in exactly two places. Verified exactly as this document
said it should be: full stdout of `rama < tests/data/peptide12.pdb` is byte-identical
before and after (`cmp`, not just `smoke_rama_pdb`'s header-line regex).

**An `msan` CI job exists now**, but scoped down from a whole-suite run after two real
findings during setup:
- `tools/cdlearn.cpp` had `#include<omp.h>` unconditionally, for a header nothing in
  the file actually uses (only the `#pragma omp parallel for` compiler pragma appears,
  no `omp_*` runtime call anywhere) — g++ ships its own `omp.h` so this was invisible
  there, but clang has none without `libomp-dev`, which isn't installed. Deleted the
  dead include; this also makes the existing `build-and-test` clang leg's cdlearn
  build not silently depend on that package.
- Plain glibc `fopen()` produces a confirmed, stack-layout-dependent MSan false
  positive on this toolchain — reproduced in complete isolation (a five-line program
  outside this repo, `fopen("/etc/hostname", "r")` via the same RAII pattern
  `biasmap_initialise` uses) — because there is no MSan-instrumented libc available
  here. It fires inside `main.cpp`'s startup `fopen`, which every `func_adcp_*` test
  reaches, and would make them intermittently red for reasons with nothing to do with
  the code under test. `dssp2cm`/`rama`'s smoke tests read only from stdin and never
  call `fopen`, so the job runs the whole tree but scopes `ctest` to just
  `smoke_dssp2cm`/`smoke_rama.*` — confirmed clean across repeated local runs. This is
  a documented, evidence-based scope decision, not an oversight: a broader MSan leg
  would need an MSan-instrumented libc to be reliable, which is a bigger undertaking
  than this gap justifies on its own.

**`cdlearn`'s OpenMP path has real coverage now** (`func_cdlearn_openmp`, two
proteins, `-DADCP_OPENMP=ON` only) — and building it surfaced a second real,
previously-unknown bug on the way: `sim_params_copy` at the per-protein setup loop
(before this fix, [tools/cdlearn.cpp:729-732](tools/cdlearn.cpp#L729)) cloned the
*single*, shared `sim_params` into every `sim_params_sim[j]` slot, but nothing had
ever set that shared struct's `seq`/`sequence`/`NAA` from an actually-read protein —
`main.cpp` always calls `update_sim_params_from_chain()` after every chain read;
`cdlearn.cpp`'s protein-loading loop never did. Every real CD-learning iteration
therefore either crashed (`"sequence is not present in sim_params for MC lookup
table calculation"`, on the very first `move()` call) or — had the shared struct
happened to hold a stale value from some other run — silently used the *wrong*
protein's sequence for every slot but one. Fixed by calling
`update_sim_params_from_chain(&(all_chains[i]), &(sim_params_sim[i]))` for each
protein individually, right after its own `sim_params_copy`. Confirmed: the same
minimal PDB+`.icm`+learn-string pipeline built for finding #8 above now runs a
complete, real CD-learning iteration end to end (`FINISHED!`) instead of stopping
at that error. Per the earlier decision, **this does not fix `rand()`'s
thread-safety** — `move()` (`src/metropolis.cpp`) calls the global, non-thread-safe
`rand()` ~28 times, and every OpenMP thread calls it concurrently. Attempted to
confirm the race directly with ThreadSanitizer; TSan itself failed to initialize in
this sandboxed environment (`FATAL: ThreadSanitizer: unexpected memory mapping`,
a known TSan/container limitation unrelated to this codebase) — the race stands on
the same first-principles reasoning this document already gave it (shared mutable
global state, concurrent unsynchronized access), just not independently confirmed by
a sanitizer. `func_cdlearn_openmp` only asserts the parallel path completes without
crashing or hanging, run 5× locally with no flakiness — it deliberately does not
assert determinism, since there is none to assert.

**Full regression, every configuration this pass could plausibly affect:** default
build 18/18, ASan+UBSan 18/18, the new `msan` job's actual scope (`smoke_dssp2cm`/
`smoke_rama*`) 3/3 clean across repeated runs, `-DADCP_MPI=ON` 19/19,
`-DADCP_OPENMP=ON` 19/19 (18 plus the new `func_cdlearn_openmp`). None of these
changes touch any previously-tested behavior.

## Phase 2 progress — MPI unblocked

Steps 8 and 9 above were "partly done" for one reason: `src/CMakeLists.txt` made
`PARALLEL` `PRIVATE` to `adcp_mpi`'s `main.cpp` only. `adcp_core` — `nested.cpp`,
`checkpoint_io.cpp`, `flex.cpp`, `probe.cpp` — was compiled **once, without
`PARALLEL`**, and that same non-PARALLEL library was linked into both `adcp` and
`adcp_mpi`. So `adcp_mpi` really was an MPI `main` wired to a serial
`nestedsampling()`; `mpi_send_chain`/`mpi_rec_chain` had never been emitted into
any object file, in any configuration, in this project's history — confirmed by
`nm` on the pre-fix binary before starting this work.

**The fix, exactly as this document already prescribed:** a second static library.
`src/CMakeLists.txt` now hoists the source list into `ADCP_CORE_SOURCES` and, only
when `ADCP_MPI` is on, builds `adcp_core_mpi` from those same sources with
`PARALLEL` defined and linked against `MPI::MPI_CXX`; `adcp_mpi` links against
`adcp_core_mpi` instead of the serial `adcp_core`. `adcp`'s build, the install
target, and every non-MPI target are untouched — verified by a clean rebuild of the
default configuration and `ctest` still 16/16.

**Installing `libopenmpi-dev` and building `-DADCP_MPI=ON` for real (not the
syntax-only stub) found one real bug immediately**, the way the sanitizer leg did
in Phase 1: [include/flex.h:16](include/flex.h#L16) had `#include <mpi.h>` *inside*
its `extern "C" { ... }` block. OpenMPI's `mpi.h` pulls in C++ STL headers
(`<map>`, `<utility>`, …), and giving those C linkage doesn't compile — hundreds of
"template with C linkage" errors from libstdc++, none of them at the actual fault
line. Fixed the same way [include/peptide.h:12](include/peptide.h#L12) already
fixes the identical problem for `<vector>`: move the `#include <mpi.h>` outside the
`extern "C"` block, with a comment pointing at the precedent. None of the other
PARALLEL-guarded headers have the same bug — `nested.h`, `checkpoint_io.h`,
`probe.h`, `random16.h` all leave `#include <mpi.h>` to their `.cpp` files, at file
scope, outside any `extern "C"`.

After that one-line fix, the whole tree — including `nested.cpp`'s and
`checkpoint_io.cpp`'s `mpi_send_chain`/`mpi_rec_chain`, `setup_communicators`, and
`flex.cpp`'s `ns_for_flex_processor` — compiles clean against the real headers.
Confirmed with `nm`: `adcp_core_mpi`/`adcp_mpi` now contain `mpi_send_chain`/
`mpi_rec_chain`; `adcp_core`/`adcp` (the serial targets) contain neither.

**Runtime, not just compile-time:** ran `mpirun -np 2 adcp_mpi` on the same fold
workload as `run_fold_test.sh`. Real replica-exchange swaps occurred between the
two ranks (`swap N : rank 0 : ... <=> ...` / matching `rank 1` lines, genuine
`MPI_Send`/`MPI_Recv` traffic), and both ranks printed "successfully finished."
This exercises `main.cpp`'s swap logic, which — worth being precise about — was
already compiled under `PARALLEL` before this fix (`PARALLEL` was always `PRIVATE`
to `adcp_mpi`'s `main.cpp`); what this run actually confirms is that MPI works
end-to-end on this machine and that the newly-added `adcp_core_mpi` link doesn't
break anything `adcp_mpi` was already doing.

A second `mpirun -np 2 adcp_mpi -n` run (mirroring `run_ns_test.sh`'s fixture and
`-r 2x10`) also completed without crashing, and its single, non-duplicated stream
of `totalE` lines matching the documented serial NS sequence is *consistent with*
the collaborative `P > 1` path in `nested.cpp` running rather than two independent
per-rank simulations — but that alone wasn't verified with the same rigor as the MC
swap test. **Getting real signal there needed a deliberate multi-rank test**, and
now there is one.

### `tests/ns_evidence_mpi_test.cpp` — the P>1 companion to the closed-form NS test

`tests/ns_evidence_test.cpp` already checks `find_worst()`/`update_NS_parameters()`
against a closed-form evidence, but hardcodes `P=1` and never touches a network. The
new test drives the actual cross-rank functions instead — `collect_chains()`,
`return_and_reheap_chains()`, and through them `mpi_send_chain()`/`mpi_rec_chain()`
— under a real `mpirun`, against the same problem (`x ~ Uniform(0,1)`,
`L(x) = exp(-x/tau)`, closed-form `Z`).

The trick that makes this cheap: since a constrained sample is drawn exactly
(`Uniform(0, x*)`, no MCMC needed), a `Chain` with nothing but `.ll` set — allocated
via `allocmem_chain(&c, 1, 1)`, skipping `build_peptide_from_sequence`,
`biasmap_initialise`, and `aat_init` entirely — is everything
`collect_chains`/`return_and_reheap_chains`/`mpi_send_chain`/`mpi_rec_chain` touch
for this problem; none of them interpret a Chain's geometry, they just move
whatever is there. The initial population is built identically on every rank from
one shared RNG seed, then each rank keeps only the points `store_chain`'s real
round-robin scheme would have given it (point `i` belongs to processor
`(i-1) % P` at local index `(i-1) / P`) — that reproduces `nestedsampling()`'s
exact chainhash/cpoints shape without needing our own initial-population broadcast.

Registered as `num_ns_evidence_mpi`, `mpirun -np 4`, only built/run when
`ADCP_MPI` is on. **Passes**: two seeds, both within 3σ of the analytic evidence
(2.05σ and 0.75σ), plus a per-draw invariant (`ll > logLstar` after every
constrained sample — violated only if `logLstar` desynced across ranks, which is
exactly the failure mode a `mpi_send_chain`/`mpi_rec_chain` bug would produce).
Manually re-run at `-np 1, 2, 3, 5, 7` (including counts that don't divide `N=100`
evenly) — all pass, all in well under a second. Full `ctest` under
`-DADCP_MPI=ON` is 17/17; the default build is unaffected (16/16, unchanged) since
the new target only exists inside `if(ADCP_MPI)`.

This closes the gap the previous paragraph left open: the NS/FLEX cross-rank
chain-exchange machinery is now verified, not just compiled-and-not-crashing.

**Now done:** [.github/workflows/ci.yml](.github/workflows/ci.yml)'s `mpi` job runs
`ctest --test-dir build --output-on-failure` after the build, so `num_ns_evidence_mpi`
(and the rest of the non-docking suite) runs on every push/PR, not just when a
developer configures `-DADCP_MPI=ON` locally. One thing this needed: OpenMPI's
`mpirun` refuses to start more ranks than it detects cores unless told
`--oversubscribe`, and `num_ns_evidence_mpi` is fixed at `-np 4`
([tests/CMakeLists.txt](tests/CMakeLists.txt)) — a real risk on CI runners with
fewer than 4 vCPUs, unrelated to whether the code works. Added `--oversubscribe`
to that test's `mpirun` invocation; confirmed locally it's a no-op when cores
are plentiful (same pass, same output) and prevents the hard failure when
they're not.

## Validation gates

- **Every file conversion**: build gcc+clang, run default `ctest` (smoke + functional, ~13 tests, seconds).
- **Any change touching `energy.cpp`, `vdw.cpp`, `metropolis.cpp`, `probe.cpp`, or `main.cpp`'s swap logic**: also configure `-DADCP_DOCKING_TESTS=ON` and run both `docking` and `validation` labels locally. CI now runs both (see below), so this is a fast local pre-check rather than the only gate. `val_3q47_redock` takes **~2 minutes** — measured at 117s, 120s and 140s. An earlier draft of this document claimed ~30 min and used that to justify leaving it out of CI; that number was wrong.
- **`func_adcp_fold_determinism`**: run after every `energy.c`/`vdw.c` change without exception, **12×, not once**. It's the only automated determinism signal, and template-vs-function-pointer codegen changes can alter floating-point accumulation order in Monte Carlo sums (IEEE 754 non-associativity) — a determinism regression here is a blocking failure to investigate, not acceptable drift. **And per "MECHANISM RESOLVED" below, the `5660a33` hang is a latent bug that unrelated codegen changes can resurface, so the 12× hammer applies to every change in any file, not only container conversions.**

- **Any commit claiming "no behaviour change", including a pure rename**: prove it against
  the parent with a bit-level probe, at a workload that actually reaches the changed code.
  For docking that means **≥ 250,000 steps** — the swap machinery in `simulate()` is gated
  behind `swapMutateSteps` (200000), so anything shorter exercises none of it. A 10,000-step
  comparison reported a clean match across all 16 seeds while `5660a33`'s precision defect
  was live. The cheap oracle: seed 4242, `-r 1x250000`, redock options, ~8 s per build.

## FIXED REGRESSION — `energy.cpp`, introduced by `5660a33` (Phase 1 step 6)

**Status: symptom gone. Mechanism now IDENTIFIED as indirect, and the underlying bug
is still latent — read the "MECHANISM RESOLVED" subsection below before trusting the
fix or the hypotheses table.**

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

### MECHANISM RESOLVED — `tc` was a red herring; the real bug is still latent

> **The `sqrt` precision defect found later in this same commit is NOT this bug.** That one
> lives in `scoreSideChain`/`scoreSideChainNoClash`, which — as this section itself proves —
> are never called on the folding path. Fixing it changed no fold output: the fold ATOM md5
> is `6a438d0a673006235fccd2b1b7007ba3` before and after. **The hang below remains
> unexplained and still latent.** Two independent defects entered the tree in `5660a33`;
> only one is now understood.

Re-measured from scratch at `5660a33`, in a throwaway worktree, Release build.
Three facts, each measured:

| Measurement | Result |
|---|---|
| `5660a33` as committed (heap `std::vector` `tc`), fold command | **3 hangs / 13 runs** |
| `5660a33` with **only** `tc` changed to the fixed stack array | **0 hangs / 16 runs** |
| Entry counters compiled into `scoreSideChain` **and** `scoreSideChainNoClash`, fold command | **0 calls, in 6 of 6 runs** |

The first two reproduce this document's original bisect. The third breaks it.

**`scoreSideChain` and `scoreSideChainNoClash` are never called on the folding path.**
Every call to either sits inside `if ((int) mod_params->external_r0[0] == 1)` in
`ADenergyNoClash`; `external_r0[]` defaults to `0.0` (`params.cpp:308`, `:434`) and is
only set by `-p external=…`. `tests/run_fold_test.sh` passes only `-p Bias=NULL`. The
docking tests pass `external=5,con,1.0,1.0,Opt=1,…` — they are the only tests that
reach this code. The documented *symptom* is likewise docking-only: the
`currTargetEnergy` spin loop is inside `if (protein_model.opt == 1)` (`main.cpp:188`),
and `opt` defaults to `0`.

So changing the storage class of a buffer **in a function that never executes on this
path** moves the hang rate from ~23% to zero. `tc`'s semantics cannot be the cause.
The effect is indirect — code layout, alignment, inlining and register-allocation
ripple out from the changed function and perturb something else on the fold path.

**This is why all seven hypotheses in the table above were refuted: every one of them
assumed `tc` was executing.** They were testing the wrong function.

Consequences, and they are not comfortable:

- **The stack-`tc` fix is a layout coincidence, not a root-cause fix.** The real
  defect is still in the tree, still latent, and any change that perturbs codegen
  could resurface it — including changes in unrelated files.
- **Therefore the 12× hammer stays mandatory for every change, permanently**, not just
  for container conversions. That is the one mitigation that actually works here, and
  this document's "validate by hammering, not by reasoning" line is now load-bearing
  rather than rhetorical.
- The original claim that the fold test hangs is **correct** — do not "fix" it back.
  What is wrong is the causal attribution to `tc`.
- Finding the real bug needs the hang caught in the act. It did **not** reproduce under
  gdb (0/5) or with entry counters compiled in (0/6), the same timing-sensitivity this
  document already recorded for ASan. MSan or valgrind, still not installed here,
  remain the tools of choice.

Do not delete the stack-`tc` array — it is measurably suppressing the symptom. Just do
not believe it fixed anything.

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
| `mpi` | `-DADCP_MPI=ON`, then `ctest`. Since "Phase 2 progress — MPI unblocked" (below), this actually builds `nested.cpp`/`checkpoint_io.cpp`/`flex.cpp` with `PARALLEL` too, not just `main.cpp` — previously `adcp_mpi` linked those from the always-serial `adcp_core`. `num_ns_evidence_mpi` (`mpirun -np 4`, `--oversubscribe` for CI runners with fewer cores) is this leg's regression guard for the whole PARALLEL path. |
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

## Measured against pristine upstream — the whole migration, end to end

Every claim above compares a commit against its parent. HEAD has also been measured
against `1c1a330` ("added license stuff", 2025-10-28) — the last upstream commit,
before any of this work. Two behavioural changes were found; **one was a defect and is
now fixed, so exactly one intended difference remains.**

1. Nested sampling — `35fb3fb`, deliberate, documented above.
2. **Docking past 200,000 steps — `5660a33`, unintended. Root-caused and FIXED.**
   `sqrt` applied to a `float` expression binds to `double sqrt(double)` in C but to
   the **`float` overload** in C++, so five vector normalisations in `energy.cpp`
   silently dropped from double to single precision. One extra rounding at ~1.2e-7
   relative, amplified by the MC search into a different trajectory once the swap
   machinery engages at `swapMutateSteps`. The five sites now cast to `double`, with a
   comment saying why — **do not remove those casts.** No other math call in `src/`
   takes a `float` argument, so the blast radius was exactly those five lines.

   After the fix, docking is bit-identical to pristine at every budget tested,
   including all 16 seeds at the full 2.5M steps, and the redock returns to
   targetE -28.5508 / RMSD 0.66 Å.

   **Generalise this.** Any `<math.h>` function applied to a `float` expression changes
   meaning between C and C++, silently, with no warning — `-Wdouble-promotion` warns
   about the opposite direction and does not catch it. Every Phase 1 "rename only"
   commit is suspect on the same grounds. This one was found only because a bit-level
   comparison against pristine upstream was run at a step count large enough to reach
   the affected code; the 10,000-step comparison that preceded it reported a clean
   match.

`5660a33` is also the commit blamed for the still-latent fold hang above — but the two are
**separate defects that happened to arrive together**, and it is important not to conflate
them. The `sqrt` bug lives in `scoreSideChain`/`scoreSideChainNoClash`, which never execute
on the folding path; fixing it left the fold ATOM md5 unchanged at
`6a438d0a673006235fccd2b1b7007ba3`. The hang is still unexplained, still latent, and the
12× hammer still applies.

Verification of the fix, all on one machine, all reproducible:

| check | result |
|---|---|
| seed 4242, 250k steps | `-13.6251` — pristine value |
| docking, 16 seeds × 10,000 steps | 0/16 differ from pristine (targetE + ATOM md5) |
| docking, 16 seeds × 2,500,000 steps | 0/16 differ from pristine (targetE + ATOM md5) |
| full redock, top-ranked | targetE -28.5508, RMSD 0.66 Å |
| fold md5 / final energy | `6a438d0a…07ba3` / 5.936278 — unchanged |
| NS energy sequence, checkpoint roundtrip | unchanged |
| `ctest build` | 15/15 |
| `ctest build-asan` (ASan+UBSan) | 15/15 |
| `ctest build-dock -L "docking\|validation"` | 3/3 |

The 0.66 Å figure restores what [README.md](README.md) and
[tests/CMakeLists.txt](tests/CMakeLists.txt) always documented; the 0.72 Å recorded in the
`transopt` section above was measured while the defect was live.

A full-length pristine docking baseline exists after all: the segfaulting binary
completes under `-fsanitize=address -fsanitize-recover=address` with no source
change, because the one-past write lands in an ASan redzone. Proof that it touches
no live data, the three transparency controls, the tables, the `-fcommon` build
recipe, and the finding that the `swapChains` crash threshold is far below the
200,000 steps `7e2192d` reports: [docs/compares/1.md](docs/compares/1.md).

## Risk register

| Risk | Where | Mitigation |
|---|---|---|
| Converting a stack VLA to a heap container relocates any out-of-bounds access, turning a benign latent read into garbage data | `energy.cpp`'s `tc` (already bit us, see above); any remaining VLA→container change | Build the ASan/UBSan leg first; hammer `func_adcp_fold_determinism` 12× rather than once, since these failures are intermittent |
| `std::vector` reallocation growth differs from `realloc`, could unmask a latent buffer over/under-read previously hidden in slack space | any `realloc`→`push_back` conversion (`params.c`, `flex.c`) | ASan/UBSan CI leg (Phase 1); full functional+docking+validation run after each |
| **A `<math.h>` function applied to a `float` expression silently changes precision between C and C++** — `sqrt`, `pow`, `fabs`, `exp`, … resolve to the `float` overload in C++, where C had only the `double` one. No warning by default; `-Wdouble-promotion` warns about the opposite direction and does not catch it | any `.c`→`.cpp` rename touching float math; hit `energy.cpp` in `5660a33` | Grep every converted file for math calls whose arguments are `float`, and cast to `double` explicitly. Verified done for `src/`: after the fix, `sqrt`/`fabs`/`pow`/`exp`/`log`/`sin`/`cos`/`acos`/`atan2` have **no** remaining `float`-argument call site |
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
