# =============================================================================
# ADCP Test Suite — GROMACS-style 3-tier testing via CTest
# =============================================================================
#
# Labels:
#   smoke      — tool binaries (seconds)
#   functional — adcp folding (seconds)
#   validation — docking with grids (hours, opt-in)
#
# Usage:
#   ctest --test-dir build -L smoke                    # smoke only
#   ctest --test-dir build -L "smoke|functional"       # smoke + functional
#   ctest --test-dir build                             # all available
#   ctest --test-dir build --output-on-failure -j4     # parallel, verbose
#   ctest --test-dir build -N                          # list tests
#
# =============================================================================

set(SRC_BIN ${CMAKE_BINARY_DIR}/src/adcp)
set(TOOLS_BIN ${CMAKE_BINARY_DIR}/tools)
set(TEST_DATA ${CMAKE_CURRENT_SOURCE_DIR}/data)

# Ensure test data directory has ramaprob.data
file(MAKE_DIRECTORY ${TEST_DATA})
execute_process(
  COMMAND ${CMAKE_COMMAND} -E create_symlink
    ${CMAKE_SOURCE_DIR}/data/ramaprob.data
    ${TEST_DATA}/ramaprob.data
  RESULT_VARIABLE _sym_rc
)
if(NOT _sym_rc EQUAL 0)
  execute_process(
    COMMAND ${CMAKE_COMMAND} -E copy_if_different
      ${CMAKE_SOURCE_DIR}/data/ramaprob.data
      ${TEST_DATA}/ramaprob.data
  )
endif()

# ---------------------------------------------------------------------------
# Nível 1 — Smoke tests: tools
# ---------------------------------------------------------------------------

add_test(NAME smoke_lipa COMMAND ${TOOLS_BIN}/lipa)
set_tests_properties(smoke_lipa PROPERTIES
  LABELS "smoke"
  TIMEOUT 10
  PASS_REGULAR_EXPRESSION "MODEL.*ENDMDL"
)

add_test(NAME smoke_rama COMMAND ${TOOLS_BIN}/rama)
set_tests_properties(smoke_rama PROPERTIES
  LABELS "smoke"
  TIMEOUT 10
)

add_test(NAME smoke_cm COMMAND ${TOOLS_BIN}/cm)
set_tests_properties(smoke_cm PROPERTIES
  LABELS "smoke"
  TIMEOUT 10
)

add_test(NAME smoke_bfactor
  COMMAND sh -c "echo '' | ${TOOLS_BIN}/bfactor"
)
# bfactor: disabled — getpdb() infinite loop on stdin (known bug)
set_tests_properties(smoke_bfactor PROPERTIES
  LABELS "smoke"
  TIMEOUT 5
  DISABLED TRUE
)

add_test(NAME smoke_statistics COMMAND ${TOOLS_BIN}/statistics)
set_tests_properties(smoke_statistics PROPERTIES
  LABELS "smoke"
  TIMEOUT 10
  PASS_REGULAR_EXPRESSION "Count = 0"
)

add_test(NAME smoke_dssp2cm
  COMMAND sh -c "echo '' | ${TOOLS_BIN}/dssp2cm"
)
set_tests_properties(smoke_dssp2cm PROPERTIES
  LABELS "smoke"
  TIMEOUT 5
)

add_test(NAME smoke_oops COMMAND ${TOOLS_BIN}/oops)
set_tests_properties(smoke_oops PROPERTIES
  LABELS "smoke"
  TIMEOUT 10
)

add_test(NAME smoke_mergie COMMAND ${TOOLS_BIN}/mergie)
set_tests_properties(smoke_mergie PROPERTIES
  LABELS "smoke"
  TIMEOUT 5
  WILL_FAIL TRUE
)

add_test(NAME smoke_cdlearn COMMAND ${TOOLS_BIN}/cdlearn)
set_tests_properties(smoke_cdlearn PROPERTIES
  LABELS "smoke"
  TIMEOUT 30
  PASS_REGULAR_EXPRESSION "ABCDEFGHIGKLMNGPQRSTGVWGYZ 26"
)

# ---------------------------------------------------------------------------
# Nível 2 — Functional tests: adcp folding (no docking grids needed)
# ---------------------------------------------------------------------------

add_test(NAME func_adcp_fold_short COMMAND ${SRC_BIN}
  -r 1000x500 -t 2 APGVGV -p Bias=NULL -o output.pdb
)
set_tests_properties(func_adcp_fold_short PROPERTIES
  LABELS "functional"
  TIMEOUT 30
  WORKING_DIRECTORY ${TEST_DATA}
)

add_test(NAME func_adcp_fold_long COMMAND ${SRC_BIN}
  -r 5000x1000 -t 2 APGVGVAPGVGV -p Bias=NULL -o output.pdb
)
set_tests_properties(func_adcp_fold_long PROPERTIES
  LABELS "functional"
  TIMEOUT 60
  WORKING_DIRECTORY ${TEST_DATA}
)

# ---------------------------------------------------------------------------
# Nível 3 — Validation: docking with grid maps (opt-in)
# ---------------------------------------------------------------------------

option(ADCP_VALIDATION "Enable physical validation tests (requires grid maps)" OFF)

if(ADCP_VALIDATION)
  add_test(NAME val_adcp_docking_prepare
    COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/validation/prepare.sh
      ${SRC_BIN}
      ${CMAKE_CURRENT_SOURCE_DIR}/validation/systems.json
      ${CMAKE_CURRENT_BINARY_DIR}/validation/run
  )
  set_tests_properties(val_adcp_docking_prepare PROPERTIES
    LABELS "validation"
    TIMEOUT 30
  )

  add_test(NAME val_adcp_docking_analyze
    COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/validation/analyze.sh
      ${CMAKE_CURRENT_BINARY_DIR}/validation/run
      ${CMAKE_CURRENT_BINARY_DIR}/validation/results
  )
  set_tests_properties(val_adcp_docking_analyze PROPERTIES
    LABELS "validation"
    TIMEOUT 300
  )
endif()
