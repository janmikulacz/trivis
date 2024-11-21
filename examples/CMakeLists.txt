list(APPEND examples vis_2point)
list(APPEND examples vis_region)
list(APPEND examples vis_vertices)
list(APPEND examples vis_graph)
list(APPEND examples sandbox)
foreach (example ${examples})
    add_executable(${example} ${example}.cc)
    target_compile_definitions(${example} PRIVATE
            DEFAULT_MAP_DIR="${PROJECT_SOURCE_DIR}/data/maps"
            DEFAULT_OUT_DIR="${CMAKE_CURRENT_SOURCE_DIR}/outputs"
    )
    target_link_libraries(${example} PUBLIC TrivisPlus)
endforeach ()
if (TRIVIS_BUILD_TRIVIS_PATHFINDER)

    add_executable(pathfinder pathfinder.cc)
    target_compile_definitions(pathfinder PRIVATE
            DEFAULT_MAP_DIR="${PROJECT_SOURCE_DIR}/data/maps"
            DEFAULT_OUT_DIR="${CMAKE_CURRENT_SOURCE_DIR}/outputs"
    )
    target_link_libraries(pathfinder PUBLIC TrivisPlus TrivisPathfinder)

endif ()