add_library(Trivis STATIC
    src/trivis/clipper_int.cc
    src/trivis/clipper_int_id.cc)

target_include_directories(Trivis PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>)
target_compile_features(Trivis PUBLIC cxx_std_17)

# clipper_int.cc — plain Clipper2, no USINGZ
target_link_libraries(Trivis PRIVATE Clipper2)

# clipper_int_id.cc — needs USINGZ set per-file so it doesn't bleed into clipper.cc.
# Link Clipper2Z via $<TARGET_FILE:...> to get its compiled symbols without
# consuming its PUBLIC USINGZ compile definition across all Trivis sources.
set_source_files_properties(src/trivis/clipper_int_id.cc PROPERTIES
    COMPILE_DEFINITIONS "USINGZ")
target_link_libraries(Trivis PRIVATE $<TARGET_FILE:Clipper2Z>)
add_dependencies(Trivis Clipper2Z)
