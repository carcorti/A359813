#!/usr/bin/env bash
# oracle_test.sh — verifica a(1)..a(16) contro il b-file OEIS con MR a 12 basi.
#
# Uso: ./oracle_test.sh [THREADS]     (default 16)
#
# Tempo atteso sul Ryzen 9 7940HS: ~45 min a 16 thread, ~5h 30min a 1 thread.

set -euo pipefail

BIN=${BIN:-./A359813}
THREADS=${1:-16}
LOG=${LOG:-oracle_${THREADS}t.log}

make release

echo "Oracle suite start"
echo "  bin    : ${BIN}"
echo "  threads: ${THREADS}"
echo "  mode   : Montgomery + paranoid12 (12 basi deterministiche)"
echo "  log    : ${LOG}"

"${BIN}" --oracle --threads "${THREADS}" --backend mont --paranoid12 --log "${LOG}"

echo
echo "Oracle suite completed successfully on ${THREADS} thread(s)."
echo "Log written to: ${LOG}"
