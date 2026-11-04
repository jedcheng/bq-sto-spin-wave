#!/bin/bash
#PJM -L rscgrp=b-batch-mig
#PJM -L gpu=1
#PJM -L elapse=02:00:00
#PJM -j
# Diameter sweep, one bulk sub-job per diameter D = 10 * PJM_BULKNUM nm (paper: D = 100 ... 700 nm, sub-jobs 10-70).
# With BQ at J = 6.0, without BQ at J = 5.5 (x 1e11 A/m^2).
#   pjsub --bulk --sparam 10-70 -x TAG=bq    simulations/micromagnetic/jobs/d_sweep.sh
#   pjsub --bulk --sparam 10-70 -x TAG=no_bq simulations/micromagnetic/jobs/d_sweep.sh
cd "$PJM_O_WORKDIR"
source simulations/micromagnetic/jobs/env_example.sh
D=$(( PJM_BULKNUM * 10 ))
if [ "$TAG" = bq ]; then J=6.0; DIR=bq; else J=5.5; DIR=no_bq_5.5; fi     # folder names as on the data share
OUT=${DATA_DIR:-data/sim}/D_sweep/${DIR}/${D}
python simulations/micromagnetic/run.py --D "$D" --J "$J" --${TAG/_/-} --precision float32 --outdir "$OUT"
python analysis/postprocess_fft.py "$OUT"
