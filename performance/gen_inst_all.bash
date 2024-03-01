#!/bin/bash

maps=()
while IFS= read -r line; do
  maps+=("$line")
done < maps.txt
while IFS= read -r line; do
  maps+=("$line")
done < maps-simple.txt


i=0
for map in "${maps[@]}"
do
  i=$((i+1))
  echo "Generating instances for the ${i}-th map of ${#maps[@]}: ${map}"
  ../build-Release/performance/gen_inst --map-name "${map}" --max-points 10000
done