# list(APPEND examples vis_graph)
# list(APPEND examples vis_map_vertices)
# list(APPEND examples vis_region)
list(APPEND examples vis_2point)
foreach (eg ${examples})
    add_executable(${eg} ${eg}.cc)
    target_compile_definitions(${eg} PRIVATE
            DEFAULT_MAP_DIR="${PROJECT_SOURCE_DIR}/data/maps"
            DEFAULT_OUT_DIR="${CMAKE_CURRENT_SOURCE_DIR}/outputs"
    )
    target_link_libraries(${eg} PUBLIC TrivisPlus)
endforeach ()