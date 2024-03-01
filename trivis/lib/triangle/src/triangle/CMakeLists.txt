add_library(Triangle triangle.c)
target_include_directories(Triangle PUBLIC ${Triangle_INCLUDE_DIRS})
target_compile_options(Triangle PRIVATE -std=c17)
target_compile_options(Triangle PUBLIC -DTRILIBRARY -DANSI_DECLARATORS -DNO_TIMER -DNEVER_EXIT)