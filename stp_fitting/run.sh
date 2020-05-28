#!/bin/bash
for d in $(find $PWD -maxdepth 1 -type d)
    do 
        cd $d
        sbatch snake-stp.sh
done
