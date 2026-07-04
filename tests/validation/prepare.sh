#!/bin/sh
# =============================================================================
# ADCP Physical Validation — Prepare
# Generates input files and run scripts for docking validation tests.
# Usage: ./prepare.sh <adcp_binary> <systems.json> <run_dir>
# =============================================================================
set -eu

ADCP_BIN="${1:?Usage: $0 <adcp_binary> <systems.json> <run_dir>}"
SYSTEMS_JSON="${2:?}"
RUN_DIR="${3:?}"

mkdir -p "${RUN_DIR}"

echo "=== ADCP Validation Prepare ==="
echo "Binary:      ${ADCP_BIN}"
echo "Systems:     ${SYSTEMS_JSON}"
echo "Run dir:     ${RUN_DIR}"
echo

if ! command -v jq > /dev/null 2>&1; then
  echo "WARNING: jq not found — JSON parsing limited. Using defaults."
  HAS_JQ=0
else
  HAS_JQ=1
fi

# Write run script header
RUN_SCRIPT="${RUN_DIR}/run_all.sh"
cat > "${RUN_SCRIPT}" << 'HEADER'
#!/bin/sh
set -eu
RUN_DIR="$(cd "$(dirname "$0")" && pwd)"
cd "${RUN_DIR}"
HEADER

COUNT=0

if [ "${HAS_JQ}" -eq 1 ]; then
  TEST_COUNT=$(jq '.systems | length' "${SYSTEMS_JSON}")

  i=0
  while [ "${i}" -lt "${TEST_COUNT}" ]; do
    NAME=$(jq -r ".systems[${i}].name" "${SYSTEMS_JSON}")
    DIR=$(jq -r ".systems[${i}].dir" "${SYSTEMS_JSON}")
    ARGS=$(jq -r ".systems[${i}].adcp_args" "${SYSTEMS_JSON}")
    i=$((i + 1))

    TEST_DIR="${RUN_DIR}/${NAME}"
    mkdir -p "${TEST_DIR}"

    echo "# ${NAME}: ${ARGS}" >> "${RUN_SCRIPT}"
    echo "echo 'Running ${NAME}...'" >> "${RUN_SCRIPT}"
    echo "${ADCP_BIN} ${ARGS} -o ${NAME}_output.pdb > ${NAME}.log 2>&1 || echo 'FAIL: ${NAME} exited ' \$?" >> "${RUN_SCRIPT}"
    echo "" >> "${RUN_SCRIPT}"

    echo "  [+] ${NAME}"
    COUNT=$((COUNT + 1))
  done
else
  # Fallback: use a built-in default test
  TEST_DIR="${RUN_DIR}/ala4_basic"
  mkdir -p "${TEST_DIR}"
  echo "# ala4_basic: default docking test" >> "${RUN_SCRIPT}"
  echo "${ADCP_BIN} -r 10000x3000 -t 2 AAAA -p \"Bias=NULL,external=5,con8,2,1.0,Opt=1,0.5,0.5,-0.5\" -o ala4_output.pdb > ala4.log 2>&1 || echo 'FAIL'" >> "${RUN_SCRIPT}"
  echo "  [+] ala4_basic (default)"
  COUNT=1
fi

chmod +x "${RUN_SCRIPT}"

echo
echo "Prepared ${COUNT} test systems in ${RUN_DIR}"
echo "To run: ${RUN_SCRIPT}"
echo "To analyze: ./analyze.sh ${RUN_DIR} <results_dir>"
