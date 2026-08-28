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

## Validation gates

- **Every file conversion**: build gcc+clang, run default `ctest` (smoke + functional, ~13 tests, seconds).
- **Any change touching `energy.cpp`, `vdw.cpp`, `metropolis.cpp`, `probe.cpp`, or `main.cpp`'s swap logic**: also configure `-DADCP_DOCKING_TESTS=ON` and run both `docking` and `validation` labels locally. CI now runs both (see below), so this is a fast local pre-check rather than the only gate. `val_3q47_redock` takes **~2 minutes** — measured at 117s, 120s and 140s. An earlier draft of this document claimed ~30 min and used that to justify leaving it out of CI; that number was wrong.
- **`func_adcp_fold_determinism`**: run after every `energy.c`/`vdw.c` change without exception, **12×, not once**. It's the only automated determinism signal, and template-vs-function-pointer codegen changes can alter floating-point accumulation order in Monte Carlo sums (IEEE 754 non-associativity) — a determinism regression here is a blocking failure to investigate, not acceptable drift. **And per "MECHANISM RESOLVED" below, the `5660a33` hang is a latent bug that unrelated codegen changes can resurface, so the 12× hammer applies to every change in any file, not only container conversions.**

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
