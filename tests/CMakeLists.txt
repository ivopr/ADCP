# =============================================================================
# ADCP Test Suite
# =============================================================================
#
# Labels:
#   smoke      — tool binaries run against a committed PDB fixture (seconds)
#   functional — adcp folding: determinism + energy sanity (~1 min)
#
# Every test asserts something about its output. A test that only checks for
# exit code 0 would pass against a binary that does nothing, so there are none.
#
# Usage:
#   ctest --test-dir build -L smoke
#   ctest --test-dir build --output-on-failure -j4
#   ctest --test-dir build -N                    # list tests
# =============================================================================

set(SRC_BIN   ${CMAKE_BINARY_DIR}/src/adcp)
set(TOOLS_BIN ${CMAKE_BINARY_DIR}/tools)
set(FIXTURE   ${CMAKE_CURRENT_SOURCE_DIR}/data/peptide12.pdb)
set(RAMAPROB  ${CMAKE_SOURCE_DIR}/data/ramaprob.data)

# Scratch area in the BUILD tree — tests never write into the source tree.
set(TEST_WORK ${CMAKE_CURRENT_BINARY_DIR}/work)
file(MAKE_DIRECTORY ${TEST_WORK})
file(COPY ${RAMAPROB} DESTINATION ${TEST_WORK})

# ---------------------------------------------------------------------------
# Sanitizer output must never be masked.
# ---------------------------------------------------------------------------
# CTest's PASS_REGULAR_EXPRESSION *overrides the exit code*: if the expected
# string is printed, the test passes even when the process then aborts. Every
# smoke test below uses one, so an ASan/UBSan failure in a tools/ binary was
# invisible -- dssp2cm leaks on the very input smoke_dssp2cm feeds it, and the
# test reported PASS. FAIL_REGULAR_EXPRESSION is evaluated independently, so
# applying it to every test restores the sanitizer build as a real gate without
# disturbing the pass patterns.
#
# The PASS_REGULAR_EXPRESSIONs stay: several smoke tests pass *because* the
# binary calls stop()/exit(EXIT_FAILURE), and the pattern is what distinguishes
# "failed correctly" from "crashed".
set(SANITIZER_FAIL_RE "AddressSanitizer|LeakSanitizer|UndefinedBehaviorSanitizer|runtime error:")

# ---------------------------------------------------------------------------
# Smoke — tool binaries. Each asserts on real output, not just exit status.
# ---------------------------------------------------------------------------

# lipa: emits a PDB model built from Pauling geometry
add_test(NAME smoke_lipa
  COMMAND sh -c "${TOOLS_BIN}/lipa < /dev/null")
set_tests_properties(smoke_lipa PROPERTIES
  LABELS "smoke" TIMEOUT 10
  PASS_REGULAR_EXPRESSION "MODEL.*ENDMDL"
)

# cdlearn: prints its amino-acid alphabet and count.
# Needs ramaprob.data in the working directory: cdlearn now calls
# ramaprob_initialise() (previously missing entirely -- ramaprob/alaprob/
# glyprob stayed NULL and any real CD-learning iteration segfaulted the
# first time ramabias() dereferenced them), matching what main.cpp already
# does before any peptide is built.
add_test(NAME smoke_cdlearn COMMAND ${TOOLS_BIN}/cdlearn)
set_tests_properties(smoke_cdlearn PROPERTIES
  LABELS "smoke" TIMEOUT 30
  WORKING_DIRECTORY ${TEST_WORK}
  PASS_REGULAR_EXPRESSION "ABCDEFGHIGKLMNGPQRSTGVWGYZ 26"
)

# cdlearn: a list-file token longer than the 1024-byte buffer must be
# rejected, not overflow it. Regression guard for the unbounded
# fscanf(list_file, "%s", ...) fixed in tools/cdlearn.cpp.
add_test(NAME smoke_cdlearn_longtoken
  COMMAND ${TOOLS_BIN}/cdlearn -L ${CMAKE_CURRENT_SOURCE_DIR}/data/cdlearn_longtoken.txt -l V
)
set_tests_properties(smoke_cdlearn_longtoken PROPERTIES
  LABELS "smoke" TIMEOUT 10
  WORKING_DIRECTORY ${TEST_WORK}
  PASS_REGULAR_EXPRESSION "PDB base name from list file is too long"
)

# statistics: reports an empty sample rather than dividing by zero
add_test(NAME smoke_statistics
  COMMAND sh -c "${TOOLS_BIN}/statistics < /dev/null")
set_tests_properties(smoke_statistics PROPERTIES
  LABELS "smoke" TIMEOUT 10
  PASS_REGULAR_EXPRESSION "Count = 0"
)

# mergie: no files to merge -> must print usage and fail, not emit garbage.
# This used to be WILL_FAIL with no output check, so it passed on any non-zero
# exit -- including the (size_t)-8 malloc that argc-2 == -1 used to produce.
# WILL_FAIL also inverts FAIL_REGULAR_EXPRESSION, which would have re-masked
# sanitizer output here, so assert the message instead.
add_test(NAME smoke_mergie COMMAND ${TOOLS_BIN}/mergie)
set_tests_properties(smoke_mergie PROPERTIES
  LABELS "smoke" TIMEOUT 5
  PASS_REGULAR_EXPRESSION "Usage:.*FORMAT FILE"
)

# cm: contact map -> residue letter followed by a 0/1 row
add_test(NAME smoke_cm
  COMMAND sh -c "${TOOLS_BIN}/cm < ${FIXTURE}")
set_tests_properties(smoke_cm PROPERTIES
  LABELS "smoke" TIMEOUT 15
  PASS_REGULAR_EXPRESSION "[01]*[A-Z][01][01]+"
)

# oops: PostScript projection -> bracketed 3D vectors
add_test(NAME smoke_oops
  COMMAND sh -c "${TOOLS_BIN}/oops < ${FIXTURE}")
set_tests_properties(smoke_oops PROPERTIES
  LABELS "smoke" TIMEOUT 15
  PASS_REGULAR_EXPRESSION "\\[ *-?[0-9]+\\.[0-9]+"
)

# rama: with no input, must say so instead of printing a bogus table
add_test(NAME smoke_rama
  COMMAND sh -c "${TOOLS_BIN}/rama < /dev/null")
set_tests_properties(smoke_rama PROPERTIES
  LABELS "smoke" TIMEOUT 10
  PASS_REGULAR_EXPRESSION "No amino acids recognized"
)

# rama: must terminate on a real PDB. The fixture ends in ENDMDL, which used to
# leave getaa() returning its end-of-model sentinel forever instead of EOF.
add_test(NAME smoke_rama_pdb
  COMMAND sh -c "${TOOLS_BIN}/rama < ${FIXTURE}")
set_tests_properties(smoke_rama_pdb PROPERTIES
  LABELS "smoke" TIMEOUT 10
  PASS_REGULAR_EXPRESSION "Id +phi"
)

# dssp2cm: echoes the (empty) sequence and secondary-structure tracks
add_test(NAME smoke_dssp2cm
  COMMAND sh -c "${TOOLS_BIN}/dssp2cm < /dev/null")
set_tests_properties(smoke_dssp2cm PROPERTIES
  LABELS "smoke" TIMEOUT 5
  PASS_REGULAR_EXPRESSION "seq:"
)

# bfactor: reads the fixture and emits averaged coordinates with B-factors
add_test(NAME smoke_bfactor COMMAND ${TOOLS_BIN}/bfactor -f ${FIXTURE})
set_tests_properties(smoke_bfactor PROPERTIES
  LABELS "smoke" TIMEOUT 15
  PASS_REGULAR_EXPRESSION "ATOM.*ALA A"
)

# bfactor: empty input must terminate with an error.
# Regression guard — this used to loop forever because getpdb() returns EOF
# (-1), which is truthy, so `while (getpdb(...))` never exited.
add_test(NAME smoke_bfactor_empty
  COMMAND sh -c "${TOOLS_BIN}/bfactor < /dev/null")
set_tests_properties(smoke_bfactor_empty PROPERTIES
  LABELS "smoke" TIMEOUT 10
  PASS_REGULAR_EXPRESSION "no amino acids read from input"
)

# bfactor: a later model with MORE residues than the first must be rejected,
# not overrun ava/sqa (sized from the first model only). Regression guard for
# the heap-buffer-overflow fixed in tools/bfactor.cpp's parse_input().
add_test(NAME smoke_bfactor_multimodel
  COMMAND sh -c "${TOOLS_BIN}/bfactor < ${CMAKE_CURRENT_SOURCE_DIR}/data/bfactor_multimodel.pdb")
set_tests_properties(smoke_bfactor_multimodel PROPERTIES
  LABELS "smoke" TIMEOUT 10
  PASS_REGULAR_EXPRESSION "all models must have the same number of residues"
)

# ---------------------------------------------------------------------------
# Functional — adcp folding (no docking grids required)
# ---------------------------------------------------------------------------

# Asserts reproducibility under a fixed seed plus a finite, plausible energy.
add_test(NAME func_adcp_fold_determinism
  COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/run_fold_test.sh
          ${SRC_BIN} ${TEST_WORK}/fold ${RAMAPROB}
)
set_tests_properties(func_adcp_fold_determinism PROPERTIES
  LABELS "functional" TIMEOUT 180
  PASS_REGULAR_EXPRESSION "PASS: deterministic"
)

# Cheap coverage of the short path; the determinism test covers the long one.
add_test(NAME func_adcp_fold_short COMMAND ${SRC_BIN}
  -r 1000x500 -t 2 -s 12345 APGVGV -p Bias=NULL -o fold_short.pdb
)
set_tests_properties(func_adcp_fold_short PROPERTIES
  LABELS "functional" TIMEOUT 60
  WORKING_DIRECTORY ${TEST_WORK}
  PASS_REGULAR_EXPRESSION "successfully finished"
)

# ---------------------------------------------------------------------------
# Numerical — nested sampling against a closed-form evidence
# ---------------------------------------------------------------------------
#
# The only test in the suite that checks a computed VALUE against an
# independently known answer rather than against itself. It links adcp_core and
# drives the real find_worst()/update_NS_parameters() bookkeeping on a problem
# with an analytic Z. Uses its own RNG, not rand(), so it is libc-independent.

add_executable(ns_evidence_test ns_evidence_test.cpp)
target_link_libraries(ns_evidence_test PRIVATE adcp_core)

add_test(NAME num_ns_evidence COMMAND ns_evidence_test)
set_tests_properties(num_ns_evidence PROPERTIES
  LABELS "functional" TIMEOUT 60
  PASS_REGULAR_EXPRESSION "PASS: nested sampling recovers the analytic evidence"
)

# P>1 companion: drives the actual cross-rank collect_chains()/
# return_and_reheap_chains()/mpi_send_chain()/mpi_rec_chain() -- unblocked in
# MIGRATION.md's "Phase 2 progress -- MPI unblocked" -- under a real mpirun,
# against the same closed-form evidence. Only buildable/runnable once MPI is
# actually configured; adcp_core_mpi (the PARALLEL-compiled core) only exists
# under ADCP_MPI.
if(ADCP_MPI)
  add_executable(ns_evidence_mpi_test ns_evidence_mpi_test.cpp)
  target_compile_definitions(ns_evidence_mpi_test PRIVATE PARALLEL)
  target_link_libraries(ns_evidence_mpi_test PRIVATE adcp_core_mpi MPI::MPI_CXX)

  set(NS_MPI_RANKS 4)
  # --oversubscribe: OpenMPI's mpirun otherwise refuses to start more ranks
  # than it detects cores, which is a real risk on CI runners with fewer than
  # 4 vCPUs. No-op when there are enough cores (verified locally).
  add_test(NAME num_ns_evidence_mpi
    COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${NS_MPI_RANKS} ${MPIEXEC_PREFLAGS}
            --oversubscribe $<TARGET_FILE:ns_evidence_mpi_test> ${MPIEXEC_POSTFLAGS}
  )
  set_tests_properties(num_ns_evidence_mpi PROPERTIES
    LABELS "functional" TIMEOUT 120
    # Harmless when not running as root; needed when CI/containers do.
    ENVIRONMENT "OMPI_ALLOW_RUN_AS_ROOT=1;OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1"
    PASS_REGULAR_EXPRESSION "PASS: nested sampling \\(MPI, P=${NS_MPI_RANKS}\\) recovers the analytic evidence"
    FAIL_REGULAR_EXPRESSION "${SANITIZER_FAIL_RE}"
  )
endif()

# ---------------------------------------------------------------------------
# Functional — cdlearn's OpenMP path
# ---------------------------------------------------------------------------
#
# cdlearn's "#pragma omp parallel for" loop over proteins (tools/cdlearn.cpp)
# has had zero coverage in every build: ADCP_OPENMP defaults off, and
# smoke_cdlearn (no args) dies at "No list file given" long before reaching
# it. This test drives it for real with two proteins, so the loop actually
# parallelizes across more than one iteration.
#
# What this does NOT do: assert determinism. move() (src/metropolis.cpp)
# calls the global, non-thread-safe rand() ~28 times per call, and every
# OpenMP thread calls it concurrently -- a genuine data race on rand()'s
# internal state (glibc doesn't crash on it in practice, it just produces a
# racy, non-reproducible interleaving). Fixing that would mean threading a
# per-thread RNG state through move() and everything it calls -- code shared
# by adcp's serial path, nested sampling, and MPI -- which is far larger and
# riskier than this coverage gap justifies on its own (this project already
# has one still-unexplained latent determinism bug from a much smaller
# codegen change; see MIGRATION.md's "MECHANISM RESOLVED" section). A
# decision, not an omission: this test only asserts the parallel path runs a
# real iteration to completion without crashing or hanging.
if(ADCP_OPENMP)
  # cdlearn builds its output filename as "out." + the exact -L argument, with
  # no fopen NULL-check -- an absolute -L path makes that concatenation an
  # unwritable path, leaving a NULL FILE* that segfaults on the first write.
  # Passing a plain relative name (matching WORKING_DIRECTORY below) avoids
  # tripping that latent bug rather than fixing it here, out of this item's
  # scope (documenting + covering the OpenMP path, not auditing cdlearn.cpp
  # further).
  file(COPY
    ${CMAKE_CURRENT_SOURCE_DIR}/data/cdlearn_openmp_protA.pdb
    ${CMAKE_CURRENT_SOURCE_DIR}/data/cdlearn_openmp_protA.icm
    ${CMAKE_CURRENT_SOURCE_DIR}/data/cdlearn_openmp_protB.pdb
    ${CMAKE_CURRENT_SOURCE_DIR}/data/cdlearn_openmp_protB.icm
    ${CMAKE_CURRENT_SOURCE_DIR}/data/cdlearn_openmp_list.txt
    DESTINATION ${TEST_WORK}
  )
  add_test(NAME func_cdlearn_openmp
    COMMAND ${TOOLS_BIN}/cdlearn -L cdlearn_openmp_list.txt -l V -i 1
  )
  set_tests_properties(func_cdlearn_openmp PROPERTIES
    LABELS "functional" TIMEOUT 60
    WORKING_DIRECTORY ${TEST_WORK}
    PASS_REGULAR_EXPRESSION "FINISHED!"
    FAIL_REGULAR_EXPRESSION "${SANITIZER_FAIL_RE}"
  )
endif()

# ---------------------------------------------------------------------------
# Functional — nested sampling and checkpointing
# ---------------------------------------------------------------------------
#
# Everything above runs the plain MC folding path. Nested sampling only runs
# under -n, so until these existed nested.cpp and checkpoint_io.cpp — including
# store_chain's realloc-grow loop over the live NS population — had no coverage
# at all. Both drive that loop by feeding a multi-model PDB: one snapshot read
# and stored per MODEL record.

set(SNAPSHOTS ${CMAKE_CURRENT_SOURCE_DIR}/data/output.pdb)

add_test(NAME func_adcp_ns_smoke
  COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/run_ns_test.sh
          ${SRC_BIN} ${TEST_WORK}/ns ${RAMAPROB} ${SNAPSHOTS}
)
set_tests_properties(func_adcp_ns_smoke PROPERTIES
  LABELS "functional" TIMEOUT 300
  PASS_REGULAR_EXPRESSION "PASS: deterministic NS"
)

# Restarts are compared against each other, not against an uninterrupted run:
# the checkpoint stores the population but not the RNG state, so trajectories
# legitimately diverge after the restart point.
add_test(NAME func_adcp_checkpoint_roundtrip
  COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/run_checkpoint_test.sh
          ${SRC_BIN} ${TEST_WORK}/ckpt ${RAMAPROB} ${SNAPSHOTS}
)
set_tests_properties(func_adcp_checkpoint_roundtrip PROPERTIES
  LABELS "functional" TIMEOUT 300
  PASS_REGULAR_EXPRESSION "PASS: .* checkpoints written"
)

# ---------------------------------------------------------------------------
# Docking — the AutoDock grid path (opt-in: needs the published 3Q47 target)
# ---------------------------------------------------------------------------
#
# Everything above runs the folding path, where the receptor term is exactly
# 'extE 0.000000'. These tests are the only ones that exercise docking, which
# is how a stack overflow that crashed every production docking run survived.
#
# Enable with -DADCP_DOCKING_TESTS=ON. The target (~10 MB) is fetched once into
# the build tree and cached; tests SKIP rather than fail when it is unavailable.

option(ADCP_DOCKING_TESTS "Enable docking tests (downloads the 3Q47 target)" OFF)

# Applies to every test defined above; see SANITIZER_FAIL_RE.
set_tests_properties(
  smoke_lipa smoke_cdlearn smoke_cdlearn_longtoken smoke_statistics smoke_mergie
  smoke_cm smoke_oops
  smoke_rama smoke_rama_pdb smoke_dssp2cm smoke_bfactor smoke_bfactor_empty
  smoke_bfactor_multimodel
  func_adcp_fold_determinism func_adcp_fold_short
  func_adcp_ns_smoke func_adcp_checkpoint_roundtrip num_ns_evidence
  PROPERTIES FAIL_REGULAR_EXPRESSION "${SANITIZER_FAIL_RE}"
)

if(ADCP_DOCKING_TESTS)
  set(TARGET_CACHE ${CMAKE_CURRENT_BINARY_DIR}/target)
  set(STAGE ${TARGET_CACHE}/stage)

  add_test(NAME dock_fetch_target
    COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/fetch_target.sh ${TARGET_CACHE}
  )
  set_tests_properties(dock_fetch_target PROPERTIES
    LABELS "docking" TIMEOUT 900
    SKIP_RETURN_CODE 77
    FIXTURES_SETUP dock_target
  )

  # Regression guard for the swapChains overflow. Runs past swapMutateSteps
  # (200000) on purpose: a shorter run passes against the broken binary.
  add_test(NAME dock_3q47_smoke
    COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/run_dock_test.sh
            ${SRC_BIN} ${STAGE} ${TEST_WORK}/dock
  )
  set_tests_properties(dock_3q47_smoke PROPERTIES
    LABELS "docking" TIMEOUT 300
    SKIP_RETURN_CODE 77
    FIXTURES_REQUIRED dock_target
    PASS_REGULAR_EXPRESSION "PASS: docking ran deterministically"
  )

  # Full redocking of the native 3Q47 peptide, ranked as runADCP.py ranks.
  # Measured: 16 runs x 2.5M steps redocks to 0.66 A in ~96 s on 14 cores.
  # Threshold is the conventional 2.5 A docking success criterion.
  add_test(NAME val_3q47_redock
    COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/run_redock_validation.sh
            ${SRC_BIN} ${STAGE} ${TEST_WORK}/redock 16 2500000 2.5
  )
  set_tests_properties(val_3q47_redock PROPERTIES
    LABELS "validation" TIMEOUT 1800
    SKIP_RETURN_CODE 77
    FIXTURES_REQUIRED dock_target
    PASS_REGULAR_EXPRESSION "PASS: redocked to"
  )

  set_tests_properties(dock_fetch_target dock_3q47_smoke val_3q47_redock
    PROPERTIES FAIL_REGULAR_EXPRESSION "${SANITIZER_FAIL_RE}"
  )
endif()
