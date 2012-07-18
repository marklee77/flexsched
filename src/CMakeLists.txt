set(CMAKE_C_FLAGS "-O3 -DNO_QSORT_R" CACHE TYPE INTERNAL FORCE)

find_package(GLPK REQUIRED)

include_directories(${GLPK_INCLUDE_DIR})

add_executable(flexsched main.c util.c greedy.c linutil.c linearprog.c vecutil.c 
    vector.c hetero.c)

target_link_libraries (flexsched ${GLPK_LIBRARY})

