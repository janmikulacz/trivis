cmake_minimum_required(VERSION 3.20)
project(trivis_map_editor VERSION 0.1.0 LANGUAGES CXX)

find_package(pybind11 REQUIRED)

# ── Clipper2: prefer system installation, fall back to GitHub ─────────────────
find_package(Clipper2 CONFIG QUIET)

if(NOT Clipper2_FOUND)
    message(STATUS "Clipper2 not found via find_package — fetching from GitHub")
    include(FetchContent)
    FetchContent_Declare(
        Clipper2
        GIT_REPOSITORY https://github.com/AngusJohnson/Clipper2.git
        GIT_TAG        Clipper2_2.0.1
        SOURCE_SUBDIR  CPP
        GIT_SHALLOW    TRUE
    )
    set(CLIPPER2_UTILS     OFF CACHE BOOL "" FORCE)
    set(CLIPPER2_EXAMPLES  OFF CACHE BOOL "" FORCE)
    set(CLIPPER2_TESTS     OFF CACHE BOOL "" FORCE)
    FetchContent_MakeAvailable(Clipper2)

    # FetchContent creates Clipper2 without the :: namespace; add an alias so
    # clipper_bindings/CMakeLists.txt can always use Clipper2::Clipper2.
    if(TARGET Clipper2 AND NOT TARGET Clipper2::Clipper2)
        add_library(Clipper2::Clipper2 ALIAS Clipper2)
    endif()
    if(TARGET Clipper2Z AND NOT TARGET Clipper2::Clipper2Z)
        add_library(Clipper2::Clipper2Z ALIAS Clipper2Z)
    endif()
endif()

add_subdirectory(trivis)
add_subdirectory(clipper_bindings)
