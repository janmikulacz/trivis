#!/bin/bash

# Set the test parameters.
executable="../build-Release/performance/eval_method"
out_dir="./outputs"
results_dir="./results/csv"

# Ground truth algorithm.
algorithm_gt="cgal_triangular_expansion"

# Algorithms for evaluation.
algorithms=(
  "cgal_triangular_expansion"
  "cgal_triangular_expansion_inexact_constructions"
  "cgal_rotational_sweep"
  "cgal_rotational_sweep_inexact_constructions"
  "visilibity"
  "trivis"
)

# Load the points set names from point_sets.txt.
points_sets=()
while IFS= read -r line; do
  points_sets+=("$line")
done < points_sets.txt

# Load the maps from maps.txt.
maps=()
while IFS= read -r line; do
  maps+=("$line")
done < maps.txt


for algorithm in "${algorithms[@]}"; do
  append=0
  for map in "${maps[@]}"; do
    for points_set in "${points_sets[@]}"; do
      points_name="${map}_points_${points_set}"
      args=""
      args="${args} --map-name ${map}"
      args="${args} --points-type ${points_set}"
      args="${args} --points-name ${points_name}"
      args="${args} --algorithm-name ${algorithm}"
      args="${args} --algorithm-data-file ${out_dir}/${points_name}_${algorithm}.txt"
      args="${args} --algorithm-gt-name ${algorithm_gt}"
      args="${args} --algorithm-gt-data-file ${out_dir}/${points_name}_${algorithm_gt}.txt"
      args="${args} --out-file ${results_dir}/${algorithm}.csv"
      args="${args} --out-file-append ${append}"
      # Run the test.
      cmd="${executable} ${args}"
      echo "Running: ${cmd}"
      ${cmd}
      append=1
    done # points_set
  done # map
done # algorithm
