add_subdirectory(lib)

set(COMPILE_DEFINITIONS
        DEFAULT_MAP_DIR="${PROJECT_SOURCE_DIR}/data/maps"
        DEFAULT_MESH_DIR="${PROJECT_SOURCE_DIR}/data/meshes"
        DEFAULT_POINT_DIR="${PROJECT_SOURCE_DIR}/data/points"
)

add_executable(gen_points gen_points.cc)
target_compile_definitions(gen_points PRIVATE ${COMPILE_DEFINITIONS})
target_link_libraries(gen_points PUBLIC TrivisPlus)

add_executable(test_visilibity test_visilibity.cc)
target_compile_definitions(test_visilibity PRIVATE ${COMPILE_DEFINITIONS})
target_link_libraries(test_visilibity PUBLIC TrivisPlus VisiLibity1)

add_executable(test_cgal test_cgal.cc)
target_compile_definitions(test_cgal PRIVATE ${COMPILE_DEFINITIONS})
target_link_libraries(test_cgal PUBLIC TrivisPlus CGAL::CGAL)

add_executable(test_trivis test_trivis.cc)
target_compile_definitions(test_trivis PRIVATE ${COMPILE_DEFINITIONS})
target_link_libraries(test_trivis PUBLIC TrivisPlus)

if (TRIVIS_BUILD_TRIVIS_INEXACT)

    add_executable(test_trivis_inexact test_trivis.cc)
    target_compile_definitions(test_trivis_inexact PRIVATE ${COMPILE_DEFINITIONS})
    target_link_libraries(test_trivis_inexact PUBLIC TrivisPlusInexact)

endif()

add_executable(eval_method eval_method.cc)
target_compile_definitions(eval_method PRIVATE ${COMPILE_DEFINITIONS})
target_link_libraries(eval_method PUBLIC TrivisPlus)
