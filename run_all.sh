#!/usr/bin/env bash
# run_all.sh -- execute the publication pipeline end to end, inside the
# qualified container (see Dockerfile):
#   1. scripts/*.R  in numeric order  -- computation; outputs to processed/
#   2. figures/*.R  in name order     -- rendering;   outputs to output/figures/
#   3. tables/*.R   in name order     -- rendering;   outputs to output/tables/
# Each script runs as a fresh R process from the repository root; its console
# output is written to meta/logs/<script>.log. The run stops at the first
# failure.
#
# Usage:
#   ./run_all.sh              run everything
#   ./run_all.sh 310          start from the first script whose name begins
#                             with the argument (e.g. 310, figure_2, table_s4)
#   ./run_all.sh --list       print the execution order and exit
set -euo pipefail
cd "$(dirname "$0")"

# MSigDB gene sets: msigdbr 26.1.0 reads a pinned Zenodo release
# (msigdb.2026.1.zip, MD5 512ba99c6827141a9d471972b812d4ac) from the R user
# cache. Keep that cache inside the repository (raw/cache/R/msigdbr/) so the
# gene-set data is a versioned raw input and no run depends on the network.
export R_USER_CACHE_DIR="$PWD/raw/cache"

order=()
for f in scripts/*.R figures/*.R tables/*.R; do
  [ -e "$f" ] && order+=("$f")
done

if [ "${1:-}" = "--list" ]; then
  printf '%s\n' "${order[@]}"
  exit 0
fi

start="${1:-}"
skipping=0
[ -n "$start" ] && skipping=1

mkdir -p meta/logs
for f in "${order[@]}"; do
  name="$(basename "$f" .R)"
  if [ "$skipping" -eq 1 ]; then
    case "$name" in
      "$start"*) skipping=0 ;;
      *) continue ;;
    esac
  fi
  echo "== $f"
  if ! Rscript "$f" > "meta/logs/$name.log" 2>&1; then
    echo "FAILED: $f (see meta/logs/$name.log)" >&2
    exit 1
  fi
done
echo "== done: ${#order[@]} scripts listed, run finished"
