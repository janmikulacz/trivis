set(MeshOptim_SRC
        mesh_weights.cc
        min_weight_triangulation.cc
        triangular_meshes.cc
)
add_library(MeshOptim ${MeshOptim_SRC})
target_link_libraries(MeshOptim Trivis)
target_include_directories(MeshOptim PUBLIC ${MeshOptim_INCLUDE_DIRS})
