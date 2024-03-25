#!/bin/bash

maps=()
while IFS= read -r line; do
  maps+=("$line")
done < maps.txt

i=0
for map in "${maps[@]}"
do
  i=$((i+1))
  echo "Generating points for the ${i}-th map of ${#maps[@]}: ${map}"
  ../build-Release/performance/gen_points --map-name "${map}"
done