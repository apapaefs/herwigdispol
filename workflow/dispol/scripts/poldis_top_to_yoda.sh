#!/usr/bin/env bash
set -euo pipefail

script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)

exec python3 "$script_dir/poldis_top_to_yoda.py" \
  --unpol "${1:-dijets_unpol.top}" \
  --pol "${2:-dijets_pol.top}" \
  --out "${3:-POLDIS_MC_DIS_BREIT_ref.yoda.gz}" \
  --analysis MC_DIS_BREIT \
  --order NNLO
