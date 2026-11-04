#!/bin/bash
# Regenerate every manuscript and SM figure into figures/out/ (or $BQSTO_FIGOUT).
#   PY=/path/to/python bash figures/make_all.sh      (default: python on PATH)
# Data: $DATA_DIR (default data/sim; fetch with data/download.sh) and data/derived/, see data/README.md.
set -e
cd "$(dirname "$0")"
PY=${PY:-python}
for f in fig1 fig2 fig3 fig4 fig5 figS1 figS3; do echo "== $f"; $PY $f.py; done
echo "all figures in ${BQSTO_FIGOUT:-$(pwd)/out}"
