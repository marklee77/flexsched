#set(CMAKE_C_FLAGS "-O3 -DNO_QSORT_R" CACHE TYPE INTERNAL FORCE)
#include_directories(${CMAKE_SOURCE_DIR}/src/include)

find_package(GLPK REQUIRED)

include_directories(${GLPK_INCLUDE_DIR})

add_library(schedulers util.c greedy/greedy.c 
    linear/linutil.c linear/linearprog.c 
    vector/vecutil.c vector/vector.c vector/hetero.c)

target_link_libraries(schedulers ${GLPK_LIBRARY})
