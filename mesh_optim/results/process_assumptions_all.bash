#!/bin/bash

DATA_DIR="../experiments/outputs/"

python3 process_assumptions1_time.py -i $DATA_DIR/assumptions
python3 process_assumptions2_weights.py -i $DATA_DIR/assumptions