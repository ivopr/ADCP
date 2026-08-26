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
# Smoke — tool binaries. Each asserts on real output, not just exit status.
# ---------------------------------------------------------------------------

# lipa: emits a PDB model built from Pauling geometry
add_test(NAME smoke_lipa COMMAND ${TOOLS_BIN}/lipa)
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
add_test(NAME smoke_statistics COMMAND ${TOOLS_BIN}/statistics)
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
add_test(NAME smoke_rama COMMAND ${TOOLS_BIN}/rama)
set_tests_properties(smoke_rama PROPERTIES
  LABELS "smoke" TIMEOUT 10
  PASS_REGULAR_EXPRESSION "No amino acids recognized"
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
