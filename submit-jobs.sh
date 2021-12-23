#!/bin/bash
#SBATCH --account=rrg-aspuru
#SBATCH --gres=gpu:1              # Number of GPU(s) per node
#SBATCH --cpus-per-task=6         # CPU cores/threads
#SBATCH --mem=12000M              # memory per node
#SBATCH --time=0-01:00            # time (DD-HH:MM)

module load python/3.6 scipy-stack
module load StdEnv/2020 gcc/9.3.0
module load rdkit/2021.03.3

source  ~/env/reinvent/bin/activate

time python data_structs.py
time python train_prior.py --num-epochs 50 # --verbose
time python main.py --num-steps 10 --batch-size 500

deactivate
