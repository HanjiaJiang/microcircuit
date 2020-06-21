#!/bin/bash

mkdir -p ~/hambach/copy

cp -r microcircuit/ scans/ *.sh cluster.json config.yml run_network.py create_snakes.py Snakefile Snakefile_template \
      ~/hambach/copy/
