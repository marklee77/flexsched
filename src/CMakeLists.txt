#set(CMAKE_C_FLAGS "-O3 -DNO_QSORT_R" CACHE TYPE INTERNAL FORCE)
set(CMAKE_C_FLAGS "-g -DNO_QSORT_R" CACHE TYPE INTERNAL FORCE)
include_directories(${CMAKE_SOURCE_DIR}/src/include)

add_subdirectory(lib)
add_subdirectory(tools)
