set(VisiLibity1_SRC
        visilibity.cpp
)
add_library(VisiLibity1 ${VisiLibity1_SRC})
target_include_directories(VisiLibity1 PUBLIC ${VisiLibity1_INCLUDE_DIRS})