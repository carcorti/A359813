#!/usr/bin/env bash
# campaign_a18.sh — lancio campagna a(18) con Montgomery + 12 basi.
#
# Uso: ./campaign_a18.sh [THREADS] [BACKEND] [CKPT_PATH] [LOG_PATH]
# Default: 14 thread (margine termico), backend mont,
#          checkpoint A359813_state_a18.bin, log a18.log
#
# Tempo atteso sul Ryzen 9 7940HS: ~12h 50min (mont, 14 thread).
#
# Nota: a(18) richiede che ORACLE_A[17] = 40386401580 sia nel sorgente
# (gia' presente in questa versione).

set -euo pipefail

THREADS=${1:-14}
BACKEND=${2:-mont}
CKPT=${3:-A359813_state_a18.bin}
LOG=${4:-a18.log}

make release

echo "Launching a(18) campaign"
echo "  threads   : ${THREADS}"
echo "  backend   : ${BACKEND}"
echo "  bases     : 12 (deterministic)"
echo "  checkpoint: ${CKPT}"
echo "  log       : ${LOG}"
./A359813 --run 18 --threads "${THREADS}" --resume --backend "${BACKEND}" --paranoid12 --checkpoint "${CKPT}" --log "${LOG}"
