#!/usr/bin/env bash
# campaign_a17.sh — lancio campagna a(17) con Montgomery + 12 basi.
#
# Uso: ./campaign_a17.sh [THREADS] [BACKEND] [CKPT_PATH] [LOG_PATH]
# Default: 16 thread, backend mont, checkpoint A359813_state_main.bin, log a17_main.log
#
# Tempo atteso sul Ryzen 9 7940HS: ~2h 20min (mont) / ~4h 10min (u128).

set -euo pipefail

THREADS=${1:-16}
BACKEND=${2:-mont}
CKPT=${3:-A359813_state_main.bin}
LOG=${4:-a17_main.log}

make release

echo "Launching a(17) campaign"
echo "  threads   : ${THREADS}"
echo "  backend   : ${BACKEND}"
echo "  bases     : 12 (deterministic)"
echo "  checkpoint: ${CKPT}"
echo "  log       : ${LOG}"
./A359813 --run 17 --threads "${THREADS}" --resume --backend "${BACKEND}" --paranoid12 --checkpoint "${CKPT}" --log "${LOG}"
