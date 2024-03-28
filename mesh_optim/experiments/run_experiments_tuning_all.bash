#!/bin/bash

python3 run_experiments_tuning1_sample_dist.py -G -E
python3 run_experiments_tuning2_mwt.py
python3 run_experiments_tuning3_mwt2.py
python3 run_experiments_tuning4_weights.py -G -E