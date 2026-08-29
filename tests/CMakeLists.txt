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

# cdlearn: prints its amino-acid alphabet and count
add_test(NAME smoke_cdlearn COMMAND ${TOOLS_BIN}/cdlearn)
set_tests_properties(smoke_cdlearn PROPERTIES
  LABELS "smoke" TIMEOUT 30
  PASS_REGULAR_EXPRESSION "ABCDEFGHIGKLMNGPQRSTGVWGYZ 26"
)

# statistics: reports an empty sample rather than dividing by zero
add_test(NAME smoke_statistics
  COMMAND sh -c "${TOOLS_BIN}/statistics < /dev/null")
set_tests_properties(smoke_statistics PROPERTIES
  LABELS "smoke" TIMEOUT 10
  PASS_REGULAR_EXPRESSION "Count = 0"
)

# mergie: no MPI output to merge -> must fail rather than emit garbage
add_test(NAME smoke_mergie COMMAND ${TOOLS_BIN}/mergie)
set_tests_properties(smoke_mergie PROPERTIES
  LABELS "smoke" TIMEOUT 5 WILL_FAIL TRUE
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
  smoke_lipa smoke_cdlearn smoke_statistics smoke_mergie smoke_cm smoke_oops
  smoke_rama smoke_rama_pdb smoke_dssp2cm smoke_bfactor smoke_bfactor_empty
  func_adcp_fold_determinism func_adcp_fold_short
  func_adcp_ns_smoke func_adcp_checkpoint_roundtrip
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
