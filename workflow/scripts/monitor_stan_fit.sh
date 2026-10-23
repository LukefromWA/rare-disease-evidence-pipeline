#!/usr/bin/env bash
# Monitor a CmdStan console log while a model runs.
# Usage: monitor_stan_fit.sh /path/to/console.log
set -euo pipefail

LOG="${1:?Usage: monitor_stan_fit.sh <cmdstan_console_log>}"

if [ ! -f "$LOG" ]; then
  echo "Log not found: $LOG" >&2
  exit 1
fi

echo "Monitoring: $LOG"
echo "Press Ctrl+C to stop watching. This does not stop Stan."

tail -n 40 -f "$LOG" | sed -u \
  -e 's/Chain [0-9]*/\x1b[1;36m&\x1b[0m/g' \
  -e 's/Iteration:.*/\x1b[1;32m&\x1b[0m/g' \
  -e 's/[Ee]rror/\x1b[1;31m&\x1b[0m/g' \
  -e 's/[Ww]arning/\x1b[1;33m&\x1b[0m/g' \
  -e 's/divergen/\x1b[1;33m&\x1b[0m/gI' \
  -e 's/treedepth/\x1b[1;33m&\x1b[0m/gI' \
  -e 's/finished/\x1b[1;32m&\x1b[0m/gI'
