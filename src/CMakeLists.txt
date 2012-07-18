set(CMAKE_C_FLAGS "-O3 -DNO_QSORT_R" CACHE TYPE INTERNAL FORCE)
set(CMAKE_C_LINK_FLAGS "-lglpk" CACHE TYPE INTERNAL FORCE)

include_directories(${CMAKE_SOURCE_DIR}/src /usr/include/glpk)

add_executable(flexched main.c util.c greedy.c linutil.c linearprog.c vecutil.c 
    vector.c hetero.c)

#CFLAGS = -I/opt/ilog/cplex112/include -O3 -DNO_QSORT_R -DCPLEX
#LDFLAGS = -L/opt/ilog/cplex112/bin/x86-64_debian4.0_4.1 -lcplex112
#CFLAGS = -I$(HOME)/progs/base/include -O3 -DNO_QSORT_R
#LDFLAGS = -L$(HOME)/progs/base/lib -lglpk

