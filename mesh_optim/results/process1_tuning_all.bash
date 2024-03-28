#!/bin/bash

DATA_DIR="../experiments/outputs/"

python3 process1.py -i $DATA_DIR/tuning1_sample_dist -o csv/tuning1_sample_dist.csv
python3 process1.py --only_info_files -i $DATA_DIR/tuning2_mwt -o csv/tuning2_mwt.csv
for file in $(ls $DATA_DIR/tuning3_mwt2); do
  python3 process1.py --only_info_files -i $DATA_DIR/tuning3_mwt2/$file -o csv/tuning3_mwt2_$file.csv
done
python3 process1.py -i $DATA_DIR/tuning4_weights -o csv/tuning4_weights.csv