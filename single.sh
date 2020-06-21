#!/bin/bash
#SBATCH --job-name single-run
#SBATCH --time 02:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
mkdir -p data
#SBATCH -o ./data/out.txt
#SBATCH -e ./data/err.txt
#SBATCH --ntasks-per-node=1
source activate nest-log
python $PWD/run_network.py 
cp $PWD/* $PWD/data/

