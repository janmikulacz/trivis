cmake_minimum_required(VERSION 3.10)
project(Clipper2 VERSION 1.0.6 LANGUAGES C CXX)

set(CMAKE_POSITION_INDEPENDENT_CODE ON)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set_property(GLOBAL PROPERTY USE_FOLDERS ON)

option(CLIPPER2_UTILS "Build utilities" ON)
option(CLIPPER2_EXAMPLES "Build examples" ON)
option(CLIPPER2_TESTS "Build tests" ON)
option(USE_EXTERNAL_GTEST "Use system-wide installed GoogleTest" OFF)
option(BUILD_SHARED_LIBS "Build shared libs" OFF)

set(Clipper2_SRC
        clipper.engine.cpp
        clipper.offset.cpp
        clipper.rectclip.cpp
)

add_library(Clipper2 ${Clipper2_SRC})
target_include_directories(Clipper2 PUBLIC ${Clipper2_INCLUDE_DIRS})
