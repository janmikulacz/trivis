pybind11_add_module(trivis_clipper src/bindings.cc)

target_compile_features(trivis_clipper PRIVATE cxx_std_17)
target_link_libraries(trivis_clipper PRIVATE Clipper2::Clipper2)

# Place the built .so next to the app so Python can import it without setting PYTHONPATH.
set_target_properties(trivis_clipper PROPERTIES
    LIBRARY_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/app"
)
