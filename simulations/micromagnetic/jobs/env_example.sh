# Environment used for the published runs on a PJM cluster with H100 GPUs (site-specific; adapt to your system).
# Any CUDA build of PyTorch 2.3 with the latest magnum.np should work; see simulations/README.md.
module load cuda/12.2.2 cudnn/8.9.7 nccl/2.22.3 pytorch-cuda/2.3.1-12.2.2   # provides torch 2.3.1+cu121
source "$HOME/venv/magnumnp/bin/activate"                                     # venv with magnum.np
