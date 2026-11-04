#!/bin/bash
#PJM -L rscgrp=b-batch-mig
#PJM -L gpu=1
#PJM -L elapse=01:00:00
#PJM -j
# Current sweep, D = 320 nm, one bulk sub-job per current density J = PJM_BULKNUM / 10 (in 1e11 A/m^2).
# The paper uses J = 1.0 ... 20.0 in steps of 0.1 (sub-jobs 10-200); TAG = bq or no_bq.
#   pjsub --bulk --sparam 10-200 -x TAG=bq    simulations/micromagnetic/jobs/j_sweep.sh
#   pjsub --bulk --sparam 10-200 -x TAG=no_bq simulations/micromagnetic/jobs/j_sweep.sh
cd "$PJM_O_WORKDIR"
source simulations/micromagnetic/jobs/env_example.sh
J=$(python -c "print(f'{${PJM_BULKNUM}/10:.1f}')")
OUT=${DATA_DIR:-data/sim}/J_sweep/${TAG}/${J}
python simulations/micromagnetic/run.py --D 320 --J "$J" --${TAG/_/-} --precision float32 --outdir "$OUT"
python analysis/postprocess_fft.py "$OUT"
