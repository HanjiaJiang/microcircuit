#!/bin/bash
#SBATCH -o ./out-ltc.txt
#SBATCH -e ./err-ltc.txt
#SBATCH --job-name hj-ltc
#SBATCH --mem=4G
#SBATCH --time=0
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 1
#SBATCH --partition=hambach,blaustein,hamstein
#SBATCH --mail-type=FAIL # notifications for job done & fail
#SBATCH --mail-user=h.jiang@fz-juelich.de

source activate nest-log
mkdir -p out

snakemake --unlock\

snakemake --jobs 20\
          --cluster-config cluster.json\
          --cluster "sbatch -n {cluster.n} \
                            -o out/ltc.{jobid}.out \
                            --cpus-per-task {cluster.nCPUs} \
                            --mem {cluster.mem} \
                            --time {cluster.time} \
                            --partition=hambach,blaustein,hamstein \
                            --mail-type=FAIL \
                            --mail-user=h.jiang@fz-juelich.de"\
          --configfile config.yml\
          --use-conda\
          --rerun-incomplete\
          --jobname "ltc.{jobid}"
