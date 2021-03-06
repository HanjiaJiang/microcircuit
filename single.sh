#!/bin/bash
#SBATCH --job-name hj-sg
#SBATCH --time 02:00:00
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=0
mkdir -p data
#SBATCH -o ./data/out-sg.txt
#SBATCH -e ./data/err-sg.txt

source activate nest-log
python $PWD/run_network.py ${TMPDIR}
cp $PWD/* $PWD/data/
