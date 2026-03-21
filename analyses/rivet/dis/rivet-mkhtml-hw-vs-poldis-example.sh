#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 2 || $# -gt 3 ]]; then
  echo "usage: $0 <herwig-yoda> <poldis-ref-yoda> [output-dir]" >&2
  exit 1
fi

herwig_yoda=$1
reference_yoda=$2
output_dir=${3:-plots_mc_vs_poldis}
reflabel=${REFLABEL:-POLDIS NNLO}
ratiolabel=${RATIOLABEL:-MC/POLDIS}

exec rivet-mkhtml \
  "$herwig_yoda" \
  "$reference_yoda" \
  --reflabel="$reflabel" \
  --ratiolabel "$ratiolabel" \
  -o "$output_dir"
