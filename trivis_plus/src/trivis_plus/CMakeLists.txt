
unset(ENV{BOOST_ROOT}) # if BOOST_ROOT is set to a different version than found by find_package, it may cause problems
set(Boost_USE_STATIC_LIBS OFF)
set(Boost_USE_MULTITHREADED ON)
set(Boost_USE_STATIC_RUNTIME OFF)
add_definitions(-DBOOST_LOG_DYN_LINK)
find_package(Boost 1.82 REQUIRED COMPONENTS program_options system filesystem log_setup log)

set(TrivisPlus_SRC
        data_loading/load_map.cc
        data_loading/load_mesh.cc
        data_loading/load_visibility_region.cc
        drawing/random_colors.cc
        drawing/drawing.cc
        drawing/fancy_drawing.cc
        utils/log.cc
        utils/ptree_json.cc
)
add_library(TrivisPlus ${TrivisPlus_SRC})
target_link_libraries(TrivisPlus PUBLIC Trivis ${Boost_LIBRARIES} cairo)
target_include_directories(TrivisPlus PUBLIC ${TrivisPlus_INCLUDE_DIRS})

if (TRIVIS_BUILD_TRIVIS_INEXACT)

    add_library(TrivisPlusInexact ${TrivisPlus_SRC})
    target_link_libraries(TrivisPlusInexact PUBLIC TrivisInexact ${Boost_LIBRARIES} cairo)
    target_include_directories(TrivisPlusInexact PUBLIC ${TrivisPlus_INCLUDE_DIRS})

endif()
