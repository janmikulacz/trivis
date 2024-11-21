unset(ENV{BOOST_ROOT}) # if BOOST_ROOT is set to a different version than found by find_package, it may cause problems
find_package(Boost 1.82 REQUIRED)

set(TrivisPathfinder_SRC
        astar.cc
        boost_graph.cc
        dijkstra.cc
        trivis_pathfinder.cc
        utils.cc
        visibility_planner.cc
)

add_library(TrivisPathfinder ${TrivisPathfinder_SRC})
target_link_libraries(TrivisPathfinder PUBLIC Trivis ${Boost_LIBRARIES})
target_include_directories(TrivisPathfinder PUBLIC ${TrivisPathfinder_INCLUDE_DIRS})

if (TRIVIS_BUILD_TRIVIS_INEXACT)

    add_library(TrivisPathfinderInexact ${TrivisPathfinder_SRC})
    target_link_libraries(TrivisPathfinderInexact PUBLIC TrivisInexact ${Boost_LIBRARIES})
    target_include_directories(TrivisPathfinderInexact PUBLIC ${TrivisPathfinder_INCLUDE_DIRS})

endif ()
