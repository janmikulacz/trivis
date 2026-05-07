add_library(Trivis STATIC src/trivis/clipper.cc)
target_include_directories(Trivis PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>)
target_link_libraries(Trivis PRIVATE Clipper2)
target_compile_features(Trivis PUBLIC cxx_std_17)

add_library(TrivisZ STATIC src/trivis/clipper_id.cc)
target_include_directories(TrivisZ PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>)
target_link_libraries(TrivisZ PRIVATE Clipper2Z)
target_compile_features(TrivisZ PUBLIC cxx_std_17)
