#!/bin/bash

# Set the test parameters.
build_dir="../build-Release"
out_dir="./outputs"
max_points=1000

visilibity_epsilon="1e-12"
test_frameworks=("trivis" "cgal" "visilibity")
cgal_algorithms=("triangular_expansion" "rotational_sweep")
cgal_versions=("default" "inexact_constructions")
declare -A timeouts
timeouts["trivis"]=20
timeouts["cgal_triangular_expansion"]=20
timeouts["cgal_triangular_expansion_inexact_constructions"]=20
timeouts["cgal_rotational_sweep"]=60
timeouts["cgal_rotational_sweep_inexact_constructions"]=60
timeouts["visilibity"]=3600 # doesn't need timeout

# Load the points set names from point_sets.txt.
point_sets=()
while IFS= read -r line; do
  point_sets+=("$line")
done < point_sets.txt

# Load the maps from maps.txt.
maps=()
while IFS= read -r line; do
  maps+=("$line")
done < maps.txt

for test_framework in "${test_frameworks[@]}"; do

  executable="${build_dir}/performance/test_${test_framework}"

  # Set the algorithms to test.
  algorithms=("default")
  if [ "${test_framework}" == "cgal" ]; then
    algorithms=("${cgal_algorithms[@]}")
  fi

  for algorithm in "${algorithms[@]}"; do

    algorithm_versions=("default")
    if [ "${test_framework}" == "cgal" ]; then
      algorithm_versions=("${cgal_versions[@]}")
    fi

    for algorithm_version in "${algorithm_versions[@]}"; do

      method_name="${test_framework}"
      if [ "${algorithm}" != "default" ]; then
        method_name="${method_name}_${algorithm}"
      fi
      if [ "${algorithm_version}" != "default" ]; then
        method_name="${method_name}_${algorithm_version}"
      fi

      timeout="${timeouts[${method_name}]}"

      for point_set in "${point_sets[@]}"; do

        for map in "${maps[@]}"; do

          # Prepare the arguments.
          points_name="${map}_points_${point_set}"
          output_file="${out_dir}/${points_name}_${method_name}.txt"
          args=""
          args="${args} --map-name ${map}"
          args="${args} --points-name ${points_name}"
          args="${args} --max-points ${max_points}"
          if [ "${algorithm}" != "default" ]; then
            args="${args} --algorithm ${algorithm}"
          fi
          if [ "${algorithm_version}" == "inexact_constructions" ]; then
            args="${args} --inexact-constructions 1"
          fi
          args="${args} --out-file ${output_file}"
          if [ "${test_framework}" == "visilibity" ]; then
            args="${args} --epsilon ${visilibity_epsilon}"
          fi

          # Run the test.
          cmd="timeout ${timeout} ${executable} ${args}"
          echo "Running: ${cmd}"
          ${cmd}
          case $? in
            124)
              message="timeout"
              ;;
            0)
              message="done"
              ;;
            *)
              message="crashed"
              ;;
          esac

          # Handle crashes.
          if [ -f "${output_file}" ]; then
            crashed=1 # Assume the program crashed until proven otherwise.
            cmd_after_crash="${cmd} --resume-after-crash 1"
            while [ "${crashed}" -eq 1 ]; do
              # Read the last line of the output file.
              last_line=$(tail -n 1 "${output_file}")
              # The first word of the last line is "DONE" or the id of the last point if the program crashed.
              last_id=$(echo "${last_line}" | awk '{print $1}')
              if [ "${last_id}" == "DONE" ]; then
                crashed=0
              else
                crashed=1
              fi
              if [ "${crashed}" -eq 1 ]; then
                line_count=$(wc -l < "${output_file}")
                if [ "${line_count}" -le 3 ]; then
                  crash_id=0
                else
                  # Check that last id is an integer number:
                  if ! [[ "${last_id}" =~ ^[0-9]+$ ]]; then
                    echo "ERROR: The last id is not an integer number (${last_id})."
                    exit 1
                  fi
                  # Record the crash to the last line of the output file.
                  crash_id=$(( last_id + 1 ))
                fi
                crash_line="${crash_id} ${message}"
                # Get the last byte of the file and if it is not a newline, add a newline.
                last_byte=$(tail -c 1 "${output_file}" | od -An -t x1)
                if [ "${last_byte}" != " 0a" ]; then
                  echo "" >> "${output_file}"
                fi
                echo "${crash_line}" >> "${output_file}"
                # Resume the program from the last point.
                echo "Crashed at point ${crash_id}. Message: ${message}. Resuming..."
                ${cmd_after_crash}
                case $? in
                  124)
                    message="timeout"
                    ;;
                  0)
                    message="done"
                    ;;
                  *)
                    message="crashed"
                    ;;
                esac
              fi
              # Repeat the loop until the output file is complete.
            done
          fi
        done # map
      done # point_set
    done # algorithm_version
  done # algorithm
done # test_framework
