#!/bin/bash

mkdir -p ~/hambach/copy

cp -r ~/git/microcircuit/microcircuit/ \
      ~/git/microcircuit/conn_probs/ \
      ~/git/microcircuit/run_network.py \
      ~/git/microcircuit/snake-*.sh \
      ~/git/microcircuit/cluster.json \
      ~/git/microcircuit/config.yml \
      ~/git/analysis/scan/snake-*.sh \
      ~/git/analysis/scan/analysis_gs.py \
      ./Snakefile \
      ~/hambach/copy
