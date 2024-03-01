set(Trivis_SRC
        bucketing/bucketing.cc
        bucketing/grid.cc
        geom/intersections.cc
        geom/poly_map.cc
        geom/robust_geometry.cc
        mesh/tri_mesh.cc
        mesh/constrained_delaunay.cc
        utils/clipper_geom.cc
        utils/clipper_utils.cc
        utils/random_points.cc
        trivis.cc
        visibility_regions.cc
)
add_library(Trivis ${Trivis_SRC})
target_link_libraries(Trivis Clipper2 RobustPredicate Triangle)
target_include_directories(Trivis PUBLIC ${Trivis_INCLUDE_DIRS})

if (TRIVIS_BUILD_TRIVIS_INEXACT)

    add_library(TrivisInexact ${Trivis_SRC})
    target_compile_definitions(TrivisInexact PUBLIC -DTRIVIS_USE_INEXACT_PREDICATES)
    target_link_libraries(TrivisInexact Clipper2 Triangle)
    target_include_directories(TrivisInexact PUBLIC ${Trivis_INCLUDE_DIRS})

endif ()
