#!/bin/sh
# =============================================================================
# ADCP Physical Validation — Analyze
# Analyzes results from previously-run docking validation tests.
# Usage: ./analyze.sh <run_dir> <results_dir>
# =============================================================================
set -eu

RUN_DIR="${1:?Usage: $0 <run_dir> <results_dir>}"
RESULTS_DIR="${2:?}"

mkdir -p "${RESULTS_DIR}"

REPORT="${RESULTS_DIR}/report.md"
CSV="${RESULTS_DIR}/stats.csv"

echo "=== ADCP Validation Analyze ==="
echo "Run dir:     ${RUN_DIR}"
echo "Results dir: ${RESULTS_DIR}"
echo

# Header
echo "# ADCP Validation Report" > "${REPORT}"
echo "" >> "${REPORT}"
echo "Generated: $(date '+%Y-%m-%d %H:%M:%S')" >> "${REPORT}"
echo "" >> "${REPORT}"
echo "| # | Test | Status | Exit | Log size | PDB size |" >> "${REPORT}"
echo "|---|------|--------|------|----------|----------|" >> "${REPORT}"

echo "test,status,exit_code,log_bytes,pdb_bytes" > "${CSV}"

TOTAL=0
PASSED=0
FAILED=0

for DIR in "${RUN_DIR}"/*/; do
  NAME=$(basename "${DIR}")
  TOTAL=$((TOTAL + 1))
  STATUS="FAIL"
  EXIT_CODE="?"

  LOG_FILE="${DIR}${NAME}.log"
  PDB_FILE="${DIR}${NAME}_output.pdb"
  LOG_SIZE=0
  PDB_SIZE=0

  if [ -f "${LOG_FILE}" ]; then
    LOG_SIZE=$(wc -c < "${LOG_FILE}" | tr -d ' ')
    if grep -q "successfully finished" "${LOG_FILE}" 2>/dev/null; then
      STATUS="PASS"
      PASSED=$((PASSED + 1))
    fi
  fi

  if [ -f "${PDB_FILE}" ]; then
    PDB_SIZE=$(wc -c < "${PDB_FILE}" | tr -d ' ')
    if grep -q "ATOM\|HETATM" "${PDB_FILE}" 2>/dev/null; then
      if [ "${STATUS}" = "FAIL" ]; then
        STATUS="PASS"
        PASSED=$((PASSED + 1))
      fi
    fi
  fi

  if [ "${STATUS}" = "FAIL" ]; then
    FAILED=$((FAILED + 1))
  fi

  echo "| ${TOTAL} | ${NAME} | ${STATUS} | ${EXIT_CODE} | ${LOG_SIZE} | ${PDB_SIZE} |" >> "${REPORT}"
  echo "${NAME},${STATUS},${EXIT_CODE},${LOG_SIZE},${PDB_SIZE}" >> "${CSV}"

  echo "  [$(echo ${STATUS} | tr 'A-Z' 'a-z')] ${NAME}"
done

# Summary
echo "" >> "${REPORT}"
echo "## Summary" >> "${REPORT}"
echo "" >> "${REPORT}"
echo "| Metric | Value |" >> "${REPORT}"
echo "|--------|-------|" >> "${REPORT}"
echo "| Total  | ${TOTAL} |" >> "${REPORT}"
echo "| Passed | ${PASSED} |" >> "${REPORT}"
echo "| Failed | ${FAILED} |" >> "${REPORT}"

if [ "${TOTAL}" -gt 0 ]; then
  PASS_PCT=$(echo "scale=1; ${PASSED} * 100 / ${TOTAL}" | bc 2>/dev/null || echo "?")
  echo "| Pass % | ${PASS_PCT}% |" >> "${REPORT}"
fi

echo "" >> "${REPORT}"
echo "---" >> "${REPORT}"
echo "Report: ${REPORT}" >> "${REPORT}"
echo "CSV:    ${CSV}" >> "${REPORT}"

echo
echo "Total: ${TOTAL} | Passed: ${PASSED} | Failed: ${FAILED}"
echo "Report: ${REPORT}"
echo "CSV:    ${CSV}"
