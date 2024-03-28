set(MeshOptim_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/include)
add_subdirectory(src/mesh_optim)

add_executable(evaluate_mesh evaluate_mesh.cc)
target_compile_definitions(evaluate_mesh PRIVATE
        DEFAULT_MAP_DIR="${PROJECT_SOURCE_DIR}/data/maps"
        DEFAULT_OUT_DIR="${CMAKE_CURRENT_SOURCE_DIR}/experiments/outputs"
)
target_link_libraries(evaluate_mesh PUBLIC MeshOptim TrivisPlus)

add_executable(generate_meshes generate_meshes.cc)
target_compile_definitions(generate_meshes PRIVATE
        DEFAULT_MAP_DIR="${PROJECT_SOURCE_DIR}/data/maps"
        DEFAULT_OUT_DIR="${CMAKE_CURRENT_SOURCE_DIR}/experiments/outputs"
        DEFAULT_WEIGHTS_DIR="${CMAKE_CURRENT_SOURCE_DIR}/experiments/weights"
)
target_link_libraries(generate_meshes PUBLIC MeshOptim TrivisPlus)