cmake_minimum_required(VERSION 3.27) # can be downgraded if needed, but without guarantees
project(trivis VERSION 0.4.0)
set(CMAKE_CXX_STANDARD 17)

# Set the default build options.
set(TRIVIS_BUILD_TRIVIS_INEXACT false)
set(TRIVIS_BUILD_TRIVIS_PLUS false)
set(TRIVIS_BUILD_EXAMPLES false)
set(TRIVIS_BUILD_PERFORMANCE false)
set(TRIVIS_BUILD_MESH_OPTIM false)

# Override the default build options with the environment variables.
if (DEFINED ENV{TRIVIS_BUILD_TRIVIS_INEXACT})
    set(TRIVIS_BUILD_TRIVIS_INEXACT $ENV{TRIVIS_BUILD_TRIVIS_INEXACT})
endif ()
if (DEFINED ENV{TRIVIS_BUILD_TRIVIS_PLUS})
    set(TRIVIS_BUILD_TRIVIS_PLUS $ENV{TRIVIS_BUILD_TRIVIS_PLUS})
endif ()
if (DEFINED ENV{TRIVIS_BUILD_EXAMPLES})
    set(TRIVIS_BUILD_EXAMPLES $ENV{TRIVIS_BUILD_EXAMPLES})
endif ()
if (DEFINED ENV{TRIVIS_BUILD_PERFORMANCE})
    set(TRIVIS_BUILD_PERFORMANCE $ENV{TRIVIS_BUILD_PERFORMANCE})
endif ()
if (DEFINED ENV{TRIVIS_BUILD_MESH_OPTIM})
    set(TRIVIS_BUILD_MESH_OPTIM $ENV{TRIVIS_BUILD_MESH_OPTIM})
endif ()

# Make corrections to the build options based on the mutual dependencies.
if (TRIVIS_BUILD_EXAMPLES OR TRIVIS_BUILD_PERFORMANCE OR TRIVIS_BUILD_MESH_OPTIM)
    set(TRIVIS_BUILD_TRIVIS_PLUS true)
endif ()

# Build TriVis by default.
add_subdirectory(trivis)

if (TRIVIS_BUILD_TRIVIS_PLUS)

    # Build TriVisPlus if requested.
    add_subdirectory(trivis_plus)

    if (TRIVIS_BUILD_EXAMPLES)

        # Build the examples if requested.
        add_subdirectory(examples)

    endif ()
    if (TRIVIS_BUILD_PERFORMANCE)

        # Build the performance tests if requested.
        add_subdirectory(performance)

    endif ()
    if (TRIVIS_BUILD_MESH_OPTIM)

        # Build the mesh optimization code if requested.
        add_subdirectory(mesh_optim)

    endif ()
endif ()
