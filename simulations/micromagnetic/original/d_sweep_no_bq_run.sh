#!/bin/bash
#PJM -L rscgrp=b-batch
#PJM -L gpu=1
#PJM -L elapse=10:00:00
#PJM -j

module load cuda/12.2.2 cudnn/8.9.7 nccl/2.22.3 pytorch-cuda/2.3.1-12.2.2


source $HOME/venv/magnumnp_v2_gpu/bin/activate



python3.11 run.py --d ${PJM_BULKNUM} --J 5.5


# python3.11 run.py --d 32 --J 6.0
