#!/bin/bash
#PJM -L rscgrp=b-batch-mig
#PJM -L gpu=1
#PJM -L elapse=01:00:00
#PJM -j
# One run: the BQ reference case of Figs. 2-4 (D = 320 nm, J = 6.0e11 A/m^2).
#   pjsub simulations/micromagnetic/jobs/single.sh      (submit from the repository root)
cd "$PJM_O_WORKDIR"
source simulations/micromagnetic/jobs/env_example.sh
python simulations/micromagnetic/run.py --D 320 --J 6.0 --bq --precision float64 --outdir data/sim/J_sweep/bq/6.0
python analysis/postprocess_fft.py data/sim/J_sweep/bq/6.0
