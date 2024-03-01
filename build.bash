#!/bin/bash

echo "Running build.bash"

echo "=========================================================================================="

if [[ -z "${TRIVIS_ROOT}" ]]; then
	echo "Env. var. TRIVIS_ROOT unset! Using \$(pwd) as root_dir."
	root_dir="$(pwd)"
else
	root_dir="${TRIVIS_ROOT}"
fi

if [[ -z "${TRIVIS_BUILD_TYPE}" ]]; then
	echo "Env. var. TRIVIS_BUILD_TYPE unset! Using Release as build_type."
	build_type="Release"
else
	build_type="${TRIVIS_BUILD_TYPE}"
fi

if [[ -z "${TRIVIS_BUILD_DIR}" ]]; then
	echo "Env. var. TRIVIS_BUILD_DIR unset! Using \$(pwd)/build-\${build_type} as build_dir."
	build_dir="$(pwd)/build-${build_type}"
else
	build_dir="${TRIVIS_BUILD_DIR}"
fi

echo "=========================================================================================="

echo "root_dir: ${root_dir}"
echo "build_type: ${build_type}"
echo "build_dir: ${build_dir}"

echo "=========================================================================================="

mkdir -p "${build_dir}"
cd "${build_dir}" || exit 1
cmake -DCMAKE_BUILD_TYPE="${build_type}" "${root_dir}"
make 

echo "=========================================================================================="

echo "DONE"

exit 0
