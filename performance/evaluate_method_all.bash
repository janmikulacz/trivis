#!/bin/bash

# Set the test parameters.
executable="../build-Release/performance/evaluate_method"
out_dir="./outputs"
results_dir="./results/csv"

# Ground truth algorithm:
algorithm_gt="cgal_triangular_expansion"

# Algorithms to be evaluated on complex maps (i.e., with holes):
algorithms_complex=(
  "trivis"
  "trivis_inexact"
  "visilibity"
  "cgal_triangular_expansion"
  "cgal_triangular_expansion_inexact_constructions"
  "cgal_rotational_sweep"
  "cgal_rotational_sweep_inexact_constructions"
)

# Algorithms to be evaluated on simple maps (i.e., without holes):
algorithms_simple=(
  "trivis"
  "trivis_inexact"
  "visilibity"
  "cgal_triangular_expansion"
  "cgal_triangular_expansion_inexact_constructions"
  "cgal_rotational_sweep"
  "cgal_rotational_sweep_inexact_constructions"
  "cgal_simple_polygon"
  "cgal_simple_polygon_inexact_constructions"
)

# Load the points types from points_types.txt.
points_types=()
while IFS= read -r line; do
  points_types+=("$line")
done < points_types.txt

# Load the maps from maps.txt and maps-simple.txt.
maps_simple=()
while IFS= read -r line; do
  maps_simple+=("$line")
done < maps-simple.txt
while IFS= read -r line; do
  maps_complex+=("$line")
done < maps.txt

for map_type in "simple" "complex"; do
  if [ "${map_type}" == "simple" ]; then
    algorithms=("${algorithms_simple[@]}")
    maps=("${maps_simple[@]}")
  else
    algorithms=("${algorithms_complex[@]}")
    maps=("${maps_complex[@]}")
  fi
  for algorithm in "${algorithms[@]}"; do
    if [ "${map_type}" == "simple" ]; then
      append=0
    else
      append=1
    fi
    for map in "${maps[@]}"; do
      for points_type in "${points_types[@]}"; do
        points_name="${map}_points_${points_type}"
        args=""
        args="${args} --map-type ${map_type}"
        args="${args} --map-name ${map}"
        args="${args} --points-type ${points_type}"
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
      done # points_type
    done # map
  done # algorithm
done # map_type



