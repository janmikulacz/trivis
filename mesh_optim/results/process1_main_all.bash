#!/bin/bash

DATA_DIR="../experiments/outputs/"

python3 process1.py -i $DATA_DIR/main1_tea -o csv/main1_tea.csv

for file in $(find $DATA_DIR/main1_tea -name "*CDT*"); do
  mkdir -p $DATA_DIR/main2_dtea${file#$DATA_DIR/main1_tea}
  cp -r ${file}/* $DATA_DIR/main2_dtea${file#$DATA_DIR/main1_tea};
done
for file in $(find $DATA_DIR/main1_tea -name "*MinLT*"); do
  mkdir -p $DATA_DIR/main2_dtea${file#$DATA_DIR/main1_tea}
  cp -r ${file}/* $DATA_DIR/main2_dtea${file#$DATA_DIR/main1_tea};
done
for file in $(find $DATA_DIR/main1_tea -name "*MinVT*"); do
  mkdir -p $DATA_DIR/main2_dtea${file#$DATA_DIR/main1_tea}
  cp -r ${file}/* $DATA_DIR/main2_dtea${file#$DATA_DIR/main1_tea};
done

python3 process1.py -i $DATA_DIR/main2_dtea -o csv/main2_dtea.csv